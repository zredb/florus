/// None wake deflection model
///
/// No wake deflection - wakes propagate straight downstream

use crate::types::{Float, Array1, Array2};
use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{GridBase, FlowField};
use std::collections::HashMap;
use ndarray::Array;

/// No wake deflection
#[derive(Debug, Clone)]
pub struct NoneVelocityDeflection {
    pub base: BaseModel,
}

impl NoneVelocityDeflection {
    pub fn new() -> Self {
        Self {
            base: BaseModel::new(HashMap::new(), "wind_vector"),
        }
    }
}

impl Default for NoneVelocityDeflection {
    fn default() -> Self {
        Self::new()
    }
}

impl DeflectionModel for NoneVelocityDeflection {
    fn prepare_function(
        &self,
        _grid: &dyn GridBase,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        Ok(HashMap::new())
    }

    fn function(
        &self,
        x: Array2,
        _y: Array2,
        _yaw_angle: Float,
        _turbulence_intensity: Float,
        _thrust_coefficient: Float,
        _rotor_diameter: Float,
        _model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array2> {
        // Return zeros with same shape as x
        let shape = x.shape();
        Ok(Array::zeros((shape[0], shape[1])))
    }
}
