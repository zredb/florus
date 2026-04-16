//! None wake velocity model
//!
//! No wake deficit - free stream conditions maintained

use crate::types::{Float, Array2, Array4};
use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{Grid, FlowField};
use std::collections::HashMap;
use ndarray::Array;

/// No wake deficit - velocity is unchanged
#[derive(Debug, Clone)]
pub struct NoneVelocity {
    pub base: BaseModel,
}

impl NoneVelocity {
    pub fn new() -> Self {
        Self {
            base: BaseModel::new(HashMap::new(), "wind_vector"),
        }
    }
}

impl Default for NoneVelocity {
    fn default() -> Self {
        Self::new()
    }
}

impl VelocityModel for NoneVelocity {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

    fn function(
        &self,
        _x: Array4,
        _y: Array4,
        _z: Array4,
        _axial_induction: Float,
        _deflection_field: Array2,
        _yaw_angle: Float,
        _turbulence_intensity: Float,
        _thrust_coefficient: Float,
        _hub_height: Float,
        _rotor_diameter: Float,
        _turbine_index: usize,
        _model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4> {
        let shape = _x.shape();
        Ok(Array::zeros((shape[0], shape[1], shape[2], shape[3])))
    }
}
