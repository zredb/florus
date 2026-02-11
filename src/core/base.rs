use serde::Deserialize;
use serde::Serialize;

use crate::core::flow_field::FlowField;
use crate::core::grid::GridBase;
/// Base wake model traits and structures
///
/// Corresponds to wake/geometry.py and base_classes in Python implementation
use crate::types::{Array1, Array2, Array4, Float};
use std::collections::HashMap;

/// Base structure for wake model parameters
#[derive(Debug, Clone)]
pub struct BaseModel {
    /// Wake model parameters
    pub parameters: HashMap<String, Float>,
    /// Dimension of the model (1D or 3D)
    pub model_dimension: String,
}

impl BaseModel {
    pub fn new(parameters: HashMap<String, Float>, model_dimension: &str) -> Self {
        Self {
            parameters,
            model_dimension: model_dimension.to_string(),
        }
    }
}

/// Velocity deficit wake model trait
pub trait VelocityModel {
    /// Prepare function - returns model-specific arguments
    fn prepare_function(
        &self,
        grid: &dyn GridBase,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>>;

    /// Main velocity deficit function
    /// turbine_index: which turbine is the wake source (for position lookup)
    fn function(
        &self,
        x: Array4,
        y: Array4,
        z: Array4,
        axial_induction: Float,
        deflection_field: Array2,
        yaw_angle: Float,
        turbulence_intensity: Float,
        thrust_coefficient: Float,
        hub_height: Float,
        rotor_diameter: Float,
        turbine_index: usize,
        model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4>;
}

/// Wake deflection model trait
pub trait DeflectionModel {
    /// Prepare function - returns model-specific arguments
    fn prepare_function(
        &self,
        grid: &dyn GridBase,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>>;

    /// Main deflection function
    fn function(
        &self,
        x: Array2,
        y: Array2,
        yaw_angle: Float,
        turbulence_intensity: Float,
        thrust_coefficient: Float,
        rotor_diameter: Float,
        model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array2>;
}

/// Wake turbulence model trait
pub trait TurbulenceModel {
    /// Prepare function - returns model-specific arguments
    fn prepare_function(
        &self,
        grid: &dyn GridBase,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>>;

    /// Main turbulence function
    fn function(
        &self,
        ambient_turbulence_intensity: Float,
        rotor_diameter_eff: Float,
        downstream_distance: Array1,
        turbine_type_parameters: &HashMap<String, Float>,
        model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array1>;
}

/// Wake combination model trait
pub trait CombinationModel {
    /// Prepare function - returns model-specific arguments
    fn prepare_function(
        &self,
        grid: &dyn GridBase,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>>;

    /// Main combination function
    fn function(&self, wake_field: &Array4, velocity_field: &Array4) -> anyhow::Result<Array4>;
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum InterpMethod {
    Nearest,
    Linear,
    Cubic,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_base_model_creation() {
        let mut params = HashMap::new();
        params.insert("ka".to_string(), 0.1);
        params.insert("kb".to_string(), 0.05);

        let model = BaseModel::new(params, "wind_vector");

        assert_eq!(model.model_dimension, "wind_vector");
        assert_eq!(model.parameters.get("ka"), Some(&0.1));
        assert_eq!(model.parameters.get("kb"), Some(&0.05));
    }

    #[test]
    fn test_base_model_parameters() {
        let mut params = HashMap::new();
        params.insert("test_param".to_string(), 42.0);

        let model = BaseModel::new(params.clone(), "scalar");

        // Verify parameter access
        assert_eq!(model.parameters.get("test_param"), Some(&42.0));

        // Verify non-existent parameter returns None
        assert_eq!(model.parameters.get("nonexistent"), None);
    }

    #[test]
    fn test_base_model_clone() {
        let mut params = HashMap::new();
        params.insert("ka".to_string(), 0.1);

        let model1 = BaseModel::new(params, "wind_vector");
        let model2 = model1.clone();

        assert_eq!(model1.model_dimension, model2.model_dimension);
        assert_eq!(model1.parameters.get("ka"), model2.parameters.get("ka"));
    }
}
