/// FLS - Freestream Linear Superposition
///
/// Combines wake fields by linear superposition
use crate::types::Array4;
use crate::core::wake::CombinationModel;
use crate::core::{Grid, FlowField};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct FLSCopied;

pub use FLSCopied as FLS;

impl CombinationModel for FLS {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

    fn function(&self, wake_field: &Array4, velocity_field: &Array4) -> anyhow::Result<Array4> {
        // Linear superposition: wake + velocity
        // This is the standard approach for combining wake deficits
        Ok(wake_field.clone() + velocity_field.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
use crate::types::Array1;
use crate::types::Array2;
use ndarray::Array;
    use approx::assert_relative_eq;

    #[test]
    fn test_fls_creation() {
        let fls = FLS;
        // Debug representation is "FLSCopied" due to the pub use alias
        assert_eq!(format!("{:?}", fls), "FLSCopied");
    }

    #[test]
    fn test_fls_function() {
        let fls = FLS;
        
        let wake_field = Array::zeros((1, 2, 3, 3));
        let velocity_field = Array::zeros((1, 2, 3, 3));
        
        let result = fls.function(&wake_field, &velocity_field);
        assert!(result.is_ok());
        
        let combined = result.unwrap();
        assert_eq!(combined.shape().len(), 4);
    }

    #[test]
    fn test_fls_linear_superposition() {
        let fls = FLS;
        
        // Create arrays with specific values
        let mut wake_field = Array::zeros((1, 1, 2, 2));
        wake_field[[0, 0, 0, 0]] = 0.3;
        wake_field[[0, 0, 1, 1]] = 0.5;
        
        let mut velocity_field = Array::zeros((1, 1, 2, 2));
        velocity_field[[0, 0, 0, 0]] = 0.1;
        velocity_field[[0, 0, 1, 1]] = 0.2;
        
        let result = fls.function(&wake_field, &velocity_field).unwrap();
        
        // FLS should add the values
        assert_relative_eq!(result[[0, 0, 0, 0]], 0.4);
        assert_relative_eq!(result[[0, 0, 1, 1]], 0.7);
    }

    #[test]
    fn test_fls_prepare_function() {
        let fls = FLS;
        let empty_4d = ndarray::Array4::zeros((0, 0, 0, 0));
        let flow_field = crate::core::FlowField {
            wind_speeds: Array1::from_vec(vec![8.0]),
            wind_directions: Array1::from_vec(vec![270.0]),
            turbulence_intensities: Array1::from_vec(vec![0.06]),
            air_density: 1.225,
            wind_shear: 0.12,
            wind_veer: 0.0,
            reference_wind_height: 90.0,
            n_findex: 1,
            u_initial_sorted: empty_4d.clone(),
            v_initial_sorted: empty_4d.clone(),
            w_initial_sorted: empty_4d.clone(),
            u_sorted: empty_4d.clone(),
            v_sorted: empty_4d.clone(),
            w_sorted: empty_4d.clone(),
            u: empty_4d.clone(),
            v: empty_4d.clone(),
            w: empty_4d.clone(),
            turbulence_intensity_field: empty_4d.clone(),
            turbulence_intensity_field_sorted: empty_4d,
        };
        let result = fls.prepare_function(&FakeGrid, &flow_field);
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }
    
    // Fake implementations for testing
    struct FakeGrid;
    impl crate::core::Grid for FakeGrid {
        fn n_turbines(&self) -> usize { 1 }
        fn n_findex(&self) -> usize { 1 }
        fn x_sorted(&self) -> &Array4 { panic!() }
        fn y_sorted(&self) -> &Array4 { panic!() }
        fn z_sorted(&self) -> &Array4 { panic!() }
        fn x_sorted_inertial_frame(&self) -> &Array4 { panic!() }
        fn y_sorted_inertial_frame(&self) -> &Array4 { panic!() }
        fn z_sorted_inertial_frame(&self) -> &Array4 { panic!() }
        fn cubature_weights(&self) -> Option<&Array2> { None }
        fn average_method(&self) -> crate::core::AveragingMethod {
            crate::core::AveragingMethod::CubicMean
        }
        fn sorted_indices(&self) -> &Array2 {
            static INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
            INDICES.get_or_init(|| Array2::zeros((1, 1)))
        }
        fn sorted_coord_indices(&self) -> &Array2 {
            static COORD_INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
            COORD_INDICES.get_or_init(|| Array2::zeros((1, 1)))
        }
        fn resolution(&self) -> usize { 1 }
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }
    }
}
