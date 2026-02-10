use crate::core::wake::CombinationModel;
use crate::core::{FlowField, GridBase};
/// MAX model - Maximum wake velocity deficit
///
/// Takes maximum wake velocity deficit to add to base flow field
use crate::types::Array4;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct MAX;

impl CombinationModel for MAX {
    fn prepare_function(
        &self,
        _grid: &dyn GridBase,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

    fn function(&self, wake_field: &Array4, velocity_field: &Array4) -> anyhow::Result<Array4> {
        // Take maximum of wake_field and velocity_field element-wise
        let shape = wake_field.shape();
        let mut result = Array4::zeros((shape[0], shape[1], shape[2], shape[3]));

        for f in 0..shape[0] {
            for t in 0..shape[1] {
                for y in 0..shape[2] {
                    for z in 0..shape[3] {
                        let w_val = wake_field[[f, t, y, z]];
                        let v_val = velocity_field[[f, t, y, z]];
                        result[[f, t, y, z]] = w_val.max(v_val);
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
    fn test_max_creation() {
        let max_model = MAX;
        assert_eq!(format!("{:?}", max_model), "MAX");
    }

    #[test]
    fn test_max_function() {
        let max_model = MAX;

        let wake_field = Array::zeros((1, 2, 3, 3));
        let velocity_field = Array::zeros((1, 2, 3, 3));

        let result = max_model.function(&wake_field, &velocity_field);
        assert!(result.is_ok());

        let combined = result.unwrap();
        assert_eq!(combined.shape().len(), 4);
    }

    #[test]
    fn test_max_takes_maximum() {
        let max_model = MAX;

        let mut wake_field = Array::zeros((1, 1, 2, 2));
        wake_field[[0, 0, 0, 0]] = 0.3;
        wake_field[[0, 0, 1, 1]] = 0.2;

        let mut velocity_field = Array::zeros((1, 1, 2, 2));
        velocity_field[[0, 0, 0, 0]] = 0.5; // Greater than wake
        velocity_field[[0, 0, 1, 1]] = 0.1; // Less than wake

        let result = max_model.function(&wake_field, &velocity_field).unwrap();

        // MAX should take the maximum value
        assert_relative_eq!(result[[0, 0, 0, 0]], 0.5);
        assert_relative_eq!(result[[0, 0, 1, 1]], 0.2);
    }

    #[test]
    fn test_max_prepare_function() {
        let max_model = MAX;
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
        let result = max_model.prepare_function(&FakeGrid, &flow_field);
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    // Fake implementations for testing
    struct FakeGrid;
    impl crate::core::GridBase for FakeGrid {
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
        fn resolution(&self) -> usize {
            1
        }
    }
}
