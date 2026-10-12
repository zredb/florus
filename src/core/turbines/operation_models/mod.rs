/// Turbine operation models
///
/// Modular turbine operation models organized by functionality:
/// - base.rs: Base types and trait definitions
/// - helpers.rs: Shared utility functions
/// - simple.rs: SimpleTurbine model
/// - cosine_loss.rs: CosineLossTurbine model
/// - simple_derating.rs: SimpleDeratingTurbine model
/// - awc.rs: AWCTurbine model (placeholder)
/// - peak_shaving.rs: PeakShavingTurbine model (placeholder)
/// - mixed.rs: MixedOperationTurbine model
/// - unified_momentum.rs: UnifiedMomentumTurbine model (full implementation)
/// - controller_dependent.rs: ControllerDependentTurbine model (framework for custom controls)
pub mod awc;
pub mod controller_dependent;
pub mod cosine_loss;
pub mod helpers;
pub mod mixed;
pub mod peak_shaving;
pub mod simple;
pub mod simple_derating;
pub mod unified_momentum;

// Re-export operation models
pub use awc::AWCTurbine;
pub use controller_dependent::ControllerDependentTurbine;
pub use cosine_loss::CosineLossTurbine;
pub use mixed::MixedOperationTurbine;
pub use peak_shaving::PeakShavingTurbine;
pub use simple::SimpleTurbine;
pub use simple_derating::SimpleDeratingTurbine;
pub use unified_momentum::UnifiedMomentumTurbine;

use crate::core::turbines::operation_models::helpers::axial_induction_from_ct;
use crate::types::Float;
use crate::types::{Array2, Array4};
use crate::Array1;

/// Parameters required for turbine power/thrust calculations
#[derive(Debug, Clone)]
pub struct TurbineParameters {
    pub power_table: LookupTable,
    pub thrust_table: LookupTable,
    pub ref_air_density: Float,
    pub cosine_loss_exponent_yaw: Float,
    pub cosine_loss_exponent_tilt: Float,
    pub ref_tilt: Float,
    pub peak_shaving_fraction: Float,
    pub peak_shaving_ti_threshold: Float,
    /// Helix (AWC) tuning parameters (Python power_thrust_table helix_*)
    pub helix_a: Float,
    pub helix_power_b: Float,
    pub helix_power_c: Float,
    pub helix_thrust_b: Float,
    pub helix_thrust_c: Float,
    /// Parameters for the controller-dependent (TUM) model; present only when
    /// the turbine YAML defines controller_dependent_turbine_parameters and a
    /// loadable cp_ct_data_file (Python power_thrust_table["cp_ct_data"]).
    pub controller_dependent: Option<ControllerDependentParams>,
}

/// Fully resolved parameters for the controller-dependent (TUM) operation
/// model (controller_dependent_operation_model.py).
#[derive(Debug, Clone)]
pub struct ControllerDependentParams {
    pub sigma: Float,
    pub beta: Float,
    pub cd: Float,
    pub cl_alfa: Float,
    /// R = rotor_diameter / 2
    pub rotor_radius: Float,
    pub rated_rpm: Float,
    pub generator_efficiency: Float,
    /// rated power in kW (turbine YAML units)
    pub rated_power_kw: Float,
    /// Reference air density (Python control_trajectory uses
    /// power_thrust_table["ref_air_density"] for the aero power, line 83/280)
    pub ref_air_density: Float,
    /// Cp over [tsr, pitch]
    pub cp_lut: ndarray::Array2<Float>,
    /// Ct over [tsr, pitch]
    pub ct_lut: ndarray::Array2<Float>,
    pub pitch_lut: ndarray::Array1<Float>,
    pub tsr_lut: ndarray::Array1<Float>,
}

/// Input context for turbine power/thrust calculations
#[derive(Debug, Clone)]
pub struct TurbineContext<'a> {
    pub velocities: &'a Array4,
    pub air_density: &'a crate::Array1,
    pub yaw_angles: Option<&'a Array2>,
    pub tilt_angles: Option<&'a Array2>,
    pub power_setpoints: Option<&'a Array2>,
    pub turbulence_intensities: Option<&'a Array4>,
    pub awc_modes: Option<&'a ndarray::Array2<String>>,
    pub awc_amplitudes: Option<&'a Array2>,
    pub cubature_weights: Option<&'a Array2>,
    pub correct_cp_ct_for_tilt: Option<&'a ndarray::Array2<bool>>,
    pub average_method: crate::core::rotor_velocity::AveragingMethod,
}

impl<'a> TurbineContext<'a> {
    /// Number of turbines derived from yaw_angles or tilt_angles shape.
    /// velcoity array second dim may be grid points (planar grid), not turbine count.
    pub fn n_turbines(&self) -> usize {
        if let Some(yaw) = self.yaw_angles {
            yaw.shape()[1]
        } else if let Some(tilt) = self.tilt_angles {
            tilt.shape()[1]
        } else {
            self.velocities.shape()[1]
        }
    }
}

/// Base trait for all turbine operation models
pub trait OperationModel: Send + Sync + std::fmt::Debug {
    fn model_name(&self) -> &'static str;
    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2>;
    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2>;
    fn axial_induction(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let ct = self.thrust_coefficient(params, ctx)?;
        Ok(axial_induction_from_ct(&ct))
    }

    fn power_loss_factor(&self, yaw_rad: Float, exponent: Float) -> Float {
        crate::utilities::cosd(yaw_rad.to_degrees()).powf(exponent)
    }
    fn clone_box(&self) -> Box<dyn OperationModel>;
}

// 为 trait object 实现 Clone
impl Clone for Box<dyn OperationModel> {
    fn clone(&self) -> Self {
        self.as_ref().clone_box()
    }
}

pub const POWER_SETPOINT_DEFAULT: Float = 1e12;
pub const POWER_SETPOINT_DISABLED: Float = 0.001;

#[derive(Debug, Clone)]
pub struct LookupTable {
    pub keys: Array1,
    pub values: Array1,
    /// Value returned outside the table's wind-speed range, matching Python's
    /// `interp1d(..., fill_value=..., bounds_error=False)` (Bug #28):
    /// 0.0 for power tables, 0.0001 for thrust tables.
    pub fill_value: Float,
}

pub type PowerTable = LookupTable;
pub type ThrustTable = LookupTable;

impl LookupTable {
    pub fn interpolate(&self, wind_speed: f64) -> f64 {
        let ws = &self.keys;
        let n = ws.len();

        if n == 0 {
            return 0.0;
        }

        if wind_speed <= ws[0] {
            return self.fill_value;
        }

        if wind_speed >= ws[n - 1] {
            return self.fill_value;
        }

        // 二分查找插值区间 [lo, hi]，满足 ws[lo] <= wind_speed <= ws[hi]
        let mut lo = 0;
        let mut hi = n - 1;

        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if wind_speed < ws[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
        }

        let x0 = ws[lo];
        let x1 = ws[hi];
        let y0 = self.values[lo];
        let y1 = self.values[hi];

        y0 + (y1 - y0) * (wind_speed - x0) / (x1 - x0)
    }
}

/// 检查功率曲线是否存在非物理的“波动” (wiggles)。
///
/// # 参数
/// * `power`: 功率数据序列 (通常对应递增的风速)。
/// * `tolerance`: 用于忽略微小数值噪声的阈值 (默认建议 0.001)。
///
/// # 返回
/// * `true`: 曲线平滑，符合物理规律 (单调上升后下降或持平)。
/// * `false`: 曲线存在异常震荡。
pub fn check_smooth_power_curve(power: &[f64], tolerance: f64) -> bool {
    // 边界情况：数据点太少无法计算差分
    if power.len() < 2 {
        return true;
    }

    // 1. 确定预期的方向变化次数 (expected_changes)
    let max_power = power.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let last_val = power.last().copied().unwrap_or(0.0);
    // 如果最后一个点 < 95% 最大功率，认为包含切出 (Cut-out)，预期变化为 2 (升->平->降)
    // 否则认为只包含上升和额定，预期变化为 1 (升->平)
    let expected_changes = if last_val < 0.95 * max_power {
        2.0
    } else {
        1.0
    };

    // 2. 计算每一步的变化方向 (dirs)
    // 逻辑: diff = power[i+1] - power[i]
    // 如果 |diff| > tolerance: sign(diff) (1.0 或 -1.0)
    // 否则: 0.0
    let dirs: Vec<f64> = power
        .windows(2)
        .map(|w| {
            let diff = w[1] - w[0];
            if diff.abs() > tolerance {
                diff.signum() // 返回 1.0, -1.0 或 0.0 (但这里 diff!=0)
            } else {
                0.0
            }
        })
        .collect();

    // 如果 dirs 为空 (原数组长度为1，虽前面已判断，但以防万一)
    if dirs.is_empty() {
        return true;
    }

    // 3. 统计方向改变的次数 (dir_changes)
    // 逻辑: sum(abs(diff(dirs)))
    let dir_changes: f64 = dirs.windows(2).map(|w| (w[1] - w[0]).abs()).sum();

    // 4. 判断是否平滑
    // 注意：由于浮点数运算，这里直接比较即可，因为 dir_changes 通常是整数 (0, 1, 2, 4...)
    dir_changes <= expected_changes
}
