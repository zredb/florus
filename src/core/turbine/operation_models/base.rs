//! Base turbine operation models
//!
//! Common types and traits for turbine operation models

use crate::types::Float;
use crate::core::turbine::turbine_type::LookupTable;
use crate::types::{Array2, Array4};
use crate::core::turbine::operation_models::helpers::axial_induction_from_ct;

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
}

impl Default for TurbineParameters {
    fn default() -> Self {
        Self {
            power_table: LookupTable {
                wind_speeds: ndarray::Array1::zeros(0),
                values: ndarray::Array1::zeros(0),
            },
            thrust_table: LookupTable {
                wind_speeds: ndarray::Array1::zeros(0),
                values: ndarray::Array1::zeros(0),
            },
            ref_air_density: 1.225,
            cosine_loss_exponent_yaw: 2.0,
            cosine_loss_exponent_tilt: 2.0,
            ref_tilt: 5.0,
            peak_shaving_fraction: 0.2,
            peak_shaving_ti_threshold: 0.15,
        }
    }
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
    pub awc_amplitudes: Option<&'a Array2>,
    pub cubature_weights: Option<&'a Array2>,
    pub correct_cp_ct_for_tilt: Option<&'a ndarray::Array2<bool>>,
    pub average_method: crate::core::rotor_velocity::AveragingMethod,
}

/// Base trait for all turbine operation models
pub trait OperationModel: Send + Sync {
    fn model_name(&self) -> &'static str;
    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2>;
    fn thrust_coefficient(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2>;
    fn axial_induction(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let ct = self.thrust_coefficient(params, ctx)?;
        Ok(axial_induction_from_ct(&ct))
    }
}

pub const POWER_SETPOINT_DEFAULT: Float = 1e12;
pub const POWER_SETPOINT_DISABLED: Float = 0.001;
