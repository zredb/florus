//! Wake-Induced Mixing Model
//!
//! Model used to generalize wake-added turbulence in the Empirical Gaussian wake model.
//! It computes the contribution of each turbine to a "wake-induced mixing" term
//! that is used in the velocity deficit and deflection models.
//!
//! Based on the Python implementation from floris-4.6

use crate::types::{Float, Array4};
use crate::core::wake::BaseModel;
use std::collections::HashMap;
use ndarray::Array;

/// Wake-Induced Mixing model parameters
///
/// Computes wake-induced mixing as a function of axial induction
/// and downstream distance.
#[derive(Debug, Clone)]
pub struct WakeInducedMixing {
    pub base: BaseModel,
    pub atmospheric_ti_gain: Float,
}

impl WakeInducedMixing {
    pub fn new(atmospheric_ti_gain: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("atmospheric_ti_gain".to_string(), atmospheric_ti_gain);

        Self {
            base: BaseModel::new(params, "wake_induced_mixing"),
            atmospheric_ti_gain,
        }
    }

    pub fn default() -> Self {
        Self::new(0.0)
    }
}

impl Default for WakeInducedMixing {
    fn default() -> Self {
        Self::new(0.0)
    }
}

impl WakeInducedMixing {
    pub fn prepare_function(&self) -> HashMap<String, Array4> {
        HashMap::new()
    }

    pub fn function(
        &self,
        axial_induction_i: &Array4,
        downstream_distance_d_i: &Array4,
    ) -> Array4 {
        let shape = axial_induction_i.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let n_y = shape[2];
        let n_z = shape[3];

        let mut wake_induced_mixing = Array::zeros((n_findex, n_turbines, n_y, n_z));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let axial_ind = axial_induction_i[[fi, ti, 0, 0]];
                let downstream_dist = downstream_distance_d_i[[fi, ti, 0, 0]];

                if downstream_dist > 0.0 {
                    let mixing = axial_ind / downstream_dist.powi(2);
                    for iy in 0..n_y {
                        for iz in 0..n_z {
                            wake_induced_mixing[[fi, ti, iy, iz]] = mixing;
                        }
                    }
                }
            }
        }

        wake_induced_mixing
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_wake_induced_mixing_creation() {
        let mixing = WakeInducedMixing::new(0.0);
        assert_eq!(mixing.atmospheric_ti_gain, 0.0);
    }

    #[test]
    fn test_wake_induced_mixing_default() {
        let mixing = WakeInducedMixing::default();
        assert_eq!(mixing.atmospheric_ti_gain, 0.0);
    }

    #[test]
    fn test_prepare_function() {
        let mixing = WakeInducedMixing::new(0.0);
        let result = mixing.prepare_function();
        assert!(result.is_empty());
    }

    #[test]
    fn test_function_basic() {
        let mixing = WakeInducedMixing::new(0.0);
        let axial_induction = Array4::from_elem((1, 2, 3, 4), 0.33);
        let downstream_dist = Array4::from_elem((1, 2, 3, 4), 5.0);

        let result = mixing.function(&axial_induction, &downstream_dist);

        let expected = 0.33 / 25.0;

        for fi in 0..1 {
            for ti in 0..2 {
                for iy in 0..3 {
                    for iz in 0..4 {
                        assert_relative_eq!(result[[fi, ti, iy, iz]], expected, epsilon = 1e-6);
                    }
                }
            }
        }
    }

    #[test]
    fn test_function_zero_distance() {
        let mixing = WakeInducedMixing::new(0.0);
        let axial_induction = Array4::from_elem((1, 1, 1, 1), 0.33);
        let downstream_dist = Array4::from_elem((1, 1, 1, 1), 0.0);

        let result = mixing.function(&axial_induction, &downstream_dist);

        assert_relative_eq!(result[[0, 0, 0, 0]], 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_function_inverse_square() {
        let mixing = WakeInducedMixing::new(0.0);
        let axial_induction = Array4::from_elem((1, 1, 1, 1), 1.0);

        let dist_1 = Array4::from_elem((1, 1, 1, 1), 1.0);
        let dist_2 = Array4::from_elem((1, 1, 1, 1), 2.0);
        let dist_4 = Array4::from_elem((1, 1, 1, 1), 4.0);

        let result_1 = mixing.function(&axial_induction, &dist_1);
        let result_2 = mixing.function(&axial_induction, &dist_2);
        let result_4 = mixing.function(&axial_induction, &dist_4);

        assert_relative_eq!(result_2[[0, 0, 0, 0]], result_1[[0, 0, 0, 0]] / 4.0, epsilon = 1e-6);
        assert_relative_eq!(result_4[[0, 0, 0, 0]], result_1[[0, 0, 0, 0]] / 16.0, epsilon = 1e-6);
    }

    #[test]
    fn test_wake_induced_mixing_with_turbine_types() {
        let mixing = WakeInducedMixing::new(0.0);
        assert_eq!(mixing.base.parameters.len(), 1);
        assert_eq!(mixing.base.parameters.get("atmospheric_ti_gain"), Some(&0.0));
    }
}
