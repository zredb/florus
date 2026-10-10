use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{FlowField, Grid};
/// None wake deflection model
///
/// No wake deflection - wakes propagate straight downstream
use crate::types::{Array1, Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

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
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        Ok(HashMap::new())
    }

    fn clone_box(&self) -> Box<dyn DeflectionModel> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn function(
        &self,
        _x_i: &Array2,
        _y_i: &Array2,
        _yaw_i: &Array1,
        _turbulence_intensity_i: &Array4,
        _thrust_coefficient_i: &Array1,
        _rotor_diameter: Float,
        grid_x: &Array4,
        _model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array4> {
        // Return zeros with same shape as the grid
        let shape = grid_x.shape();
        Ok(Array::zeros((shape[0], shape[1], shape[2], shape[3])))
    }
}
