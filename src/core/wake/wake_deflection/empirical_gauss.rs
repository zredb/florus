use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{FlowField, Grid};
/// Empirical Gauss wake deflection model
///
/// Extended Gauss model with empirical corrections for different conditions
use crate::types::{Array1, Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Empirical Gauss wake deflection model parameters
///
/// Simplified formulation independent from the velocity deficit model
/// (port of Python EmpiricalGaussVelocityDeflection).
#[derive(Debug, Clone)]
pub struct EmpiricalGaussVelocityDeflection {
    pub base: BaseModel,
    pub horizontal_deflection_gain_d: Float,
    pub vertical_deflection_gain_d: Float,
    pub deflection_rate: Float,
    pub mixing_gain_deflection: Float,
    pub yaw_added_mixing_gain: Float,
}

impl EmpiricalGaussVelocityDeflection {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        horizontal_deflection_gain_d: Float,
        vertical_deflection_gain_d: Float,
        deflection_rate: Float,
        mixing_gain_deflection: Float,
        yaw_added_mixing_gain: Float,
    ) -> Self {
        let mut params = HashMap::new();
        params.insert(
            "horizontal_deflection_gain_d".to_string(),
            horizontal_deflection_gain_d,
        );
        params.insert("deflection_rate".to_string(), deflection_rate);

        Self {
            base: BaseModel::new(params, "wind_vector"),
            horizontal_deflection_gain_d,
            vertical_deflection_gain_d,
            deflection_rate,
            mixing_gain_deflection,
            yaw_added_mixing_gain,
        }
    }

    pub fn default() -> Self {
        Self::new(3.0, -1.0, 22.0, 0.0, 0.0)
    }
}

impl Default for EmpiricalGaussVelocityDeflection {
    fn default() -> Self {
        Self::default()
    }
}

impl DeflectionModel for EmpiricalGaussVelocityDeflection {
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
        x_i: &Array2,
        _y_i: &Array2,
        yaw_i: &Array1,
        _turbulence_intensity_i: &Array1,
        thrust_coefficient_i: &Array1,
        rotor_diameter: Float,
        grid_x: &Array4,
        _model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array4> {
        // Generic-solver path: no tilt, no wake-induced mixing through this
        // interface (the dedicated empirical_gauss_solver supplies both).
        let shape = grid_x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let grid_y_dim = shape[2];
        let grid_z_dim = shape[3];

        // Deflection is evaluated at every grid point from its streamwise
        // distance to the source (Python: x_normalized = (x - x_i) * mask / D)
        let n_pts = n_turbines * grid_y_dim * grid_z_dim;
        let mut x_downstream = Array2::zeros((n_findex, n_pts));
        for fi in 0..n_findex {
            let mut idx = 0;
            for _ti in 0..n_turbines {
                for _iy in 0..grid_y_dim {
                    for _iz in 0..grid_z_dim {
                        x_downstream[[fi, idx]] = grid_x[[fi, _ti, _iy, _iz]] - x_i[[fi, 0]];
                        idx += 1;
                    }
                }
            }
        }
        let tilt_arr = Array1::zeros(n_findex);
        let mixing_arr = Array1::zeros(n_findex);
        let (deflection_y, _deflection_z) = self.deflection_yz(
            &x_downstream,
            yaw_i,
            &tilt_arr,
            &mixing_arr,
            thrust_coefficient_i,
            rotor_diameter,
        );
        let deflection_y =
            deflection_y.into_shape_with_order((n_findex, n_turbines, grid_y_dim, grid_z_dim))?;
        Ok(deflection_y)
    }
}

impl EmpiricalGaussVelocityDeflection {
    /// Emgauss deflection for per-findex downstream distances.
    ///
    /// Port of Python EmpiricalGaussVelocityDeflection.function. `x_downstream`
    /// holds the downstream distance from the wake source (m); the Python
    /// downstream mask (`x > x_i + 0.1`) zeroes the normalized distance, which
    /// makes the log term vanish. Returns (deflection_y, deflection_z), each
    /// shaped (n_findex, n).
    #[allow(clippy::too_many_arguments)]
    pub fn deflection_yz(
        &self,
        x_downstream: &Array2,
        yaw_angle: &Array1,
        tilt_angle: &Array1,
        mixing_i: &Array1,
        thrust_coefficient: &Array1,
        rotor_diameter: Float,
    ) -> (Array2, Array2) {
        let shape = x_downstream.shape();
        let n_findex = shape[0];
        let n = shape[1];

        let deflection_gain_y = self.horizontal_deflection_gain_d * rotor_diameter;
        // vertical_deflection_gain_D == -1 means "same as horizontal"
        let deflection_gain_z = if self.vertical_deflection_gain_d == -1.0 {
            deflection_gain_y
        } else {
            self.vertical_deflection_gain_d * rotor_diameter
        };

        let mut deflection_y = Array::zeros((n_findex, n));
        let mut deflection_z = Array::zeros((n_findex, n));

        for fi in 0..n_findex {
            // CW yaw for consistency with other models
            let yaw_r = -yaw_angle[fi].to_radians();
            let tilt_r = tilt_angle[fi].to_radians();
            let denom = 1.0 + self.mixing_gain_deflection * mixing_i[fi];
            let a_y = deflection_gain_y * thrust_coefficient[fi] * yaw_r / denom;
            let a_z = deflection_gain_z * thrust_coefficient[fi] * tilt_r / denom;

            for j in 0..n {
                let x_dist = x_downstream[[fi, j]];
                let x_normalized = if x_dist > 0.1 {
                    x_dist / rotor_diameter
                } else {
                    0.0
                };
                let log_term = ((x_normalized - self.deflection_rate)
                    / (x_normalized + self.deflection_rate)
                    + 2.0)
                    .ln();
                deflection_y[[fi, j]] = a_y * log_term;
                deflection_z[[fi, j]] = a_z * log_term;
            }
        }

        (deflection_y, deflection_z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_empirical_gauss_deflection_creation() {
        let empirical = EmpiricalGaussVelocityDeflection::default();
        assert_eq!(empirical.horizontal_deflection_gain_d, 3.0);
        assert_eq!(empirical.vertical_deflection_gain_d, -1.0);
        assert_eq!(empirical.deflection_rate, 22.0);
        assert_eq!(empirical.mixing_gain_deflection, 0.0);
        assert_eq!(empirical.yaw_added_mixing_gain, 0.0);
    }

    #[test]
    fn test_empirical_gauss_deflection_yz_unyawed() {
        let empirical = EmpiricalGaussVelocityDeflection::default();
        let x_dist = Array2::from_elem((1, 1), 5.0 * 126.0);
        let yaw = Array1::from_vec(vec![0.0]);
        let tilt = Array1::from_vec(vec![5.0]);
        let mixing = Array1::zeros(1);
        let ct = Array1::from_vec(vec![0.76]);

        let (defl_y, defl_z) = empirical.deflection_yz(&x_dist, &yaw, &tilt, &mixing, &ct, 126.0);

        assert_relative_eq!(defl_y[[0, 0]], 0.0, epsilon = 1e-12);
        // tilt_r = 5 deg in rad; log_term at 5D = ln((5-22)/(5+22) + 2)
        let x_n: Float = (5.0 - 22.0) / (5.0 + 22.0) + 2.0;
        let log_term_5d = x_n.ln();
        let expected_z: Float = 3.0 * 126.0 * 0.76 * tilt[0].to_radians() * log_term_5d;
        assert_relative_eq!(defl_z[[0, 0]], expected_z, epsilon = 1e-9);
    }
}
