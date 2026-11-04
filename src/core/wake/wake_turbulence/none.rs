//! None wake turbulence model
//!
//! No additional wake turbulence - ambient turbulence only

use crate::core::wake::{BaseModel, TurbulenceModel};
use crate::core::{FlowField, Grid};
use crate::types::{Array1, Array2, Array4};
use std::collections::HashMap;

/// No additional wake turbulence
#[derive(Debug, Clone)]
pub struct NoneTurbulence {
    pub base: BaseModel,
}

impl NoneTurbulence {
    pub fn new() -> Self {
        Self {
            base: BaseModel::new(HashMap::new(), "scalar"),
        }
    }
}

impl Default for NoneTurbulence {
    fn default() -> Self {
        Self::new()
    }
}

impl TurbulenceModel for NoneTurbulence {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        Ok(HashMap::new())
    }

    fn clone_box(&self) -> Box<dyn TurbulenceModel> {
        Box::new(self.clone())
    }

    fn function(
        &self,
        _ambient_turbulence_intensity: &Array4,
        x: &Array4,
        _x_i: &Array4,
        _rotor_diameters: &Array2,
        _axial_induction: &Array4,
    ) -> anyhow::Result<Array4> {
        Ok(Array4::zeros(x.dim()))
    }
}
