//! None wake velocity model
//!
//! No wake deficit - free stream conditions maintained

use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, Grid};
use crate::types::{Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

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

    fn clone_box(&self) -> Box<dyn VelocityModel> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn function(
        &self,
        _x: Array4,
        _y: Array4,
        _z: Array4,
        _axial_induction: Array2,
        _deflection_field: Array4,
        _yaw_angle: Array2,
        _turbulence_intensity: Array2,
        _thrust_coefficient: Array2,
        _hub_height: Float,
        _rotor_diameter: Float,
        _turbine_index: usize,
        _model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4> {
        let shape = _x.shape();
        Ok(Array::zeros((shape[0], shape[1], shape[2], shape[3])))
    }
}
