use crate::core::wake::CombinationModel;
use crate::core::{FlowField, Grid};
/// SOSFS - Sum of Squares Free Stream
///
/// Square root of sum of squares of wake field and velocity field
use crate::types::Array4;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct SOSFSCopied;

impl CombinationModel for SOSFSCopied {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

    fn function(&self, wake_field: &Array4, velocity_field: &Array4) -> anyhow::Result<Array4> {
        let shape = wake_field.shape();
        let mut result = Array4::zeros((shape[0], shape[1], shape[2], shape[3]));

        for f in 0..shape[0] {
            for t in 0..shape[1] {
                for y in 0..shape[2] {
                    for z in 0..shape[3] {
                        let deficit = wake_field[[f, t, y, z]];
                        let velocity_deficit = velocity_field[[f, t, y, z]];
                        result[[f, t, y, z]] = (deficit.powi(2) + velocity_deficit.powi(2)).sqrt();
                    }
                }
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;
    use crate::types::Array2;
    use approx::assert_relative_eq;
    use ndarray::Array;

    #[test]
    fn test_sosfs_creation() {
        let sosfs = SOSFS;
        assert_eq!(format!("{:?}", sosfs), "SOSFSCopied");
    }

    #[test]
    fn test_sosfs_function() {
        let sosfs = SOSFS;

        let wake_field = Array::zeros((1, 2, 3, 3));
        let velocity_field = Array::zeros((1, 2, 3, 3));

        let result = sosfs.function(&wake_field, &velocity_field);
        assert!(result.is_ok());

        let combined = result.unwrap();
        assert_eq!(combined.shape().len(), 4);
    }

    #[test]
    fn test_sosfs_sqrt_sum_squares() {
        let sosfs = SOSFS;

        let mut wake_field = Array::zeros((1, 1, 2, 2));
        wake_field[[0, 0, 0, 0]] = 3.0;
        wake_field[[0, 0, 1, 1]] = 0.0;

        let mut velocity_field = Array::zeros((1, 1, 2, 2));
        velocity_field[[0, 0, 0, 0]] = 4.0;
        velocity_field[[0, 0, 1, 1]] = 5.0;

        let result = sosfs.function(&wake_field, &velocity_field).unwrap();

        assert_relative_eq!(result[[0, 0, 0, 0]], (3.0_f64.powi(2) + 4.0_f64.powi(2)).sqrt(), epsilon = 1e-6);
        assert_relative_eq!(result[[0, 0, 1, 1]], 5.0);
    }

    #[test]
    fn test_sosfs_prepare_function() {
        let sosfs = SOSFS;
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
        let result = sosfs.prepare_function(&FakeGrid, &flow_field);
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    // Fake implementations for testing
    struct FakeGrid;
    impl crate::core::Grid for FakeGrid {
        fn n_turbines(&self) -> usize {
            1
        }
        fn n_findex(&self) -> usize {
            1
        }
        fn x_sorted(&self) -> &Array4 {
            panic!()
        }
        fn y_sorted(&self) -> &Array4 {
            panic!()
        }
        fn z_sorted(&self) -> &Array4 {
            panic!()
        }
        fn x_sorted_inertial_frame(&self) -> &Array4 {
            panic!()
        }
        fn y_sorted_inertial_frame(&self) -> &Array4 {
            panic!()
        }
        fn z_sorted_inertial_frame(&self) -> &Array4 {
            panic!()
        }
        fn cubature_weights(&self) -> Option<&Array2> {
            None
        }
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
        fn resolution(&self) -> usize {
            1
        }
        fn as_any(&self) -> &dyn std::any::Any {
            self
        }
    }
}

pub use SOSFSCopied as SOSFS;
