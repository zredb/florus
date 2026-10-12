use serde::Deserialize;
use serde::Serialize;

use crate::core::flow_field::FlowField;
use crate::core::grid::Grid;
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
        grid: &dyn Grid,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>>;

    /// Main velocity deficit function
    /// turbine_index: which turbine is the wake source (for position lookup)
    /// axial_induction, yaw_angle, thrust_coefficient are per-findex arrays of
    /// shape (n_findex, n_turbines), matching the Python implementation which
    /// broadcasts these as (n_findex, 1, 1, 1).
    /// turbulence_intensity is the full 4D TI field
    /// (n_findex, n_turbines, grid_y, grid_z); models read the source's
    /// per-point slice via `[[fi, turbine_index, iy, iz]]`, matching Python
    /// where the solver passes the source slice `turbulence_intensity_i`
    /// (n_findex, 1, grid_y, grid_z) broadcast over target turbines.
    /// deflection_field is the wake-source deflection evaluated at every grid
    /// point, shape (n_findex, n_turbines, grid_y, grid_z) — matching the
    /// Python `deflection_field_i` argument.
    /// hub_heights / rotor_diameters are the wake source's values per findex,
    /// shape (n_findex, 1), matching Python's (n_findex, 1, 1, 1) broadcast.
    #[allow(clippy::too_many_arguments)]
    fn function(
        &self,
        x: Array4,
        y: Array4,
        z: Array4,
        axial_induction: Array2,
        deflection_field: Array4,
        yaw_angle: Array2,
        turbulence_intensity: Array4,
        thrust_coefficient: Array2,
        hub_heights: &Array2,
        rotor_diameters: &Array2,
        turbine_index: usize,
        model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4>;

    /// Clone into a boxed trait object
    fn clone_box(&self) -> Box<dyn VelocityModel>;

    /// Downcast support for specialized solver access to the concrete model
    fn as_any(&self) -> &dyn std::any::Any;
}

/// Wake deflection model trait
pub trait DeflectionModel {
    /// Prepare function - returns model-specific arguments
    fn prepare_function(
        &self,
        grid: &dyn Grid,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>>;

    /// Main deflection function
    ///
    /// Computes the deflection of wake source turbine i at every grid point
    /// of the flow domain, mirroring the Python call
    /// `deflection_model.function(x_i, y_i, yaw_i, ti_i, ct_i, D, x=grid.x_sorted, ...)`.
    ///
    /// * `x_i`, `y_i`: source turbine mean position, shape `[n_findex, 1]`
    /// * `yaw_i`: per-findex effective yaw [deg], shape `[n_findex]`
    /// * `turbulence_intensity_i`: per-point TI at the source, the source's
    ///   slice of the evolving 4D TI field, shape `[n_findex, 1, grid_y, grid_z]`
    ///   (Python: `turbulence_intensity_i = turbine_turbulence_intensity[:, i:i+1]`)
    /// * `thrust_coefficient_i`: per-findex Ct of the source, shape `[n_findex]`
    /// * `rotor_diameters`: source rotor diameter per findex, shape
    ///   `[n_findex, 1]` (Python passes `(n_findex, 1, 1, 1)`)
    /// * `grid_x`: full grid x-coordinates, shape `[n_findex, n_turbines, grid_y, grid_z]`
    ///
    /// Returns: deflection field, shape `[n_findex, n_turbines, grid_y, grid_z]`
    #[allow(clippy::too_many_arguments)]
    fn function(
        &self,
        x_i: &Array2,
        y_i: &Array2,
        yaw_i: &Array1,
        turbulence_intensity_i: &Array4,
        thrust_coefficient_i: &Array1,
        rotor_diameters: &Array2,
        grid_x: &Array4,
        model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array4>;

    /// Clone into a boxed trait object
    fn clone_box(&self) -> Box<dyn DeflectionModel>;

    /// Downcast support for specialized solver access to the concrete model
    fn as_any(&self) -> &dyn std::any::Any;
}

/// Wake turbulence model trait
pub trait TurbulenceModel {
    /// Prepare function - returns model-specific arguments
    fn prepare_function(
        &self,
        grid: &dyn Grid,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>>;

    /// Main turbulence function
    ///
    /// Computes the wake-added turbulence intensity on the full grid (4D) for a
    /// single wake source turbine `i`, used by the WAT clamp in `sequential_solver`.
    ///
    /// * `ambient_turbulence_intensity`: shape `[findex, 1, 1, 1]` (inflow TI)
    /// * `x`: full grid x, shape `[findex, n_turbines, grid_y, grid_z]`
    /// * `x_i`: source turbine mean x, shape `[findex, 1, 1, 1]`
    /// * `rotor_diameters`: source rotor diameter per findex, shape
    ///   `[n_findex, 1]` (Python passes `(n_findex, 1, 1, 1)`)
    /// * `axial_induction`: source turbine axial induction, shape `[findex, 1, 1, 1]`
    fn function(
        &self,
        ambient_turbulence_intensity: &Array4,
        x: &Array4,
        x_i: &Array4,
        rotor_diameters: &Array2,
        axial_induction: &Array4,
    ) -> anyhow::Result<Array4>;

    /// Clone into a boxed trait object
    fn clone_box(&self) -> Box<dyn TurbulenceModel>;
}

/// Wake combination model trait
pub trait CombinationModel {
    /// Prepare function - returns model-specific arguments
    fn prepare_function(
        &self,
        grid: &dyn Grid,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>>;

    /// Main combination function
    fn function(&self, wake_field: &Array4, velocity_field: &Array4) -> anyhow::Result<Array4>;

    /// Clone into a boxed trait object
    fn clone_box(&self) -> Box<dyn CombinationModel>;
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
