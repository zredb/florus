//! None wake turbulence model
//!
//! No additional wake turbulence - ambient turbulence only

use crate::types::{Float, Array1};
use crate::core::wake::{BaseModel, TurbulenceModel};
use crate::core::{Grid, FlowField};
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

    fn function(
        &self,
        ambient_turbulence_intensity: Float,
        _rotor_diameter_eff: Float,
        downstream_distance: Array1,
        _turbine_type_parameters: &HashMap<String, Float>,
        _model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array1> {
        Ok(Array1::from_elem(downstream_distance.len(), ambient_turbulence_intensity))
    }
}
