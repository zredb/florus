pub mod awc;
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
pub mod base;
pub mod controller_dependent;
pub mod cosine_loss;
pub mod helpers;
pub mod mixed;
pub mod peak_shaving;
pub mod simple;
pub mod simple_derating;
pub mod unified_momentum;

// Re-export main types
pub use base::{
    OperationModel, TurbineContext, TurbineParameters, POWER_SETPOINT_DEFAULT,
    POWER_SETPOINT_DISABLED,
};

// Re-export operation models
pub use awc::AWCTurbine;
pub use controller_dependent::ControllerDependentTurbine;
pub use cosine_loss::CosineLossTurbine;
pub use mixed::MixedOperationTurbine;
pub use peak_shaving::PeakShavingTurbine;
pub use simple::SimpleTurbine;
pub use simple_derating::SimpleDeratingTurbine;
pub use unified_momentum::UnifiedMomentumTurbine;

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
