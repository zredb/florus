use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{FlowField, Grid};
/// Gauss wake deflection model
///
/// Based on Bastankhah and Porte-Agel (2016) - wake deflection for Gauss velocity model
use crate::types::{Array1, Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Gauss wake deflection model parameters
#[derive(Debug, Clone)]
pub struct GaussVelocityDeflection {
    pub base: BaseModel,
    pub ad: Float,
    pub bd: Float,
    pub alpha: Float,
    pub beta: Float,
    pub ka: Float,
    pub kb: Float,
    pub dm: Float,
    pub eps_gain: Float,
    pub use_secondary_steering: bool,
}

impl GaussVelocityDeflection {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        kd: Float,
        ad: Float,
        bd: Float,
        alpha: Float,
        beta: Float,
        ka: Float,
        kb: Float,
        dm: Float,
    ) -> Self {
        let mut params = HashMap::new();
        params.insert("kd".to_string(), kd);
        params.insert("ad".to_string(), ad);
        params.insert("bd".to_string(), bd);
        params.insert("alpha".to_string(), alpha);
        params.insert("beta".to_string(), beta);
        params.insert("ka".to_string(), ka);
        params.insert("kb".to_string(), kb);
        params.insert("dm".to_string(), dm);
        params.insert("eps_gain".to_string(), 0.2);
        params.insert("use_secondary_steering".to_string(), 1.0); // Convert bool to float for BaseModel compatibility

        Self {
            base: BaseModel::new(params, "wind_vector"),
            ad,
            bd,
            alpha,
            beta,
            ka,
            kb,
            dm,
            eps_gain: 0.2,
            use_secondary_steering: true,
        }
    }
}

impl DeflectionModel for GaussVelocityDeflection {
    fn prepare_function(
        &self,
        grid: &dyn Grid,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        let mut params = HashMap::new();

        // Add wind_veer to parameters
        let wind_veer_value = flow_field.wind_veer;
        params.insert(
            "wind_veer".to_string(),
            Array1::from_elem(grid.x_sorted().shape()[0], wind_veer_value),
        );

        Ok(params)
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
        turbulence_intensity_i: &Array4,
        thrust_coefficient_i: &Array1,
        rotor_diameters: &Array2,
        grid_x: &Array4,
        model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array4> {
        // Get wind_veer from model_args if available
        let wind_veer = if let Some(wind_veer_arr) = model_args.get("wind_veer") {
            if !wind_veer_arr.is_empty() {
                wind_veer_arr[[0]] // Using single index for Array1
            } else {
                0.0 // Default value
            }
        } else {
            0.0 // Default value
        };

        let shape = grid_x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let grid_y_dim = shape[2];
        let grid_z_dim = shape[3];

        let mut deflection = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));

        for fi in 0..n_findex {
            // Opposite sign convention in this model
            let yaw = -yaw_i[fi];
            let ct = thrust_coefficient_i[fi];
            let x_src = x_i[[fi, 0]];
            let rotor_diameter = rotor_diameters[[fi, 0]];

            // Initial velocity deficits (Python uses the freestream-normalized
            // freestream_velocity = 1; it cancels in u_r/(u + u0))
            let u_initial = 1.0; // normalized
            let tilt = 0.0; // assuming no tilt

            let ct_yaw_tilt = ct * crate::utilities::cosd(tilt) * crate::utilities::cosd(yaw);
            let sqrt_ct_yaw_tilt = (1.0 - ct_yaw_tilt).max(0.0).sqrt();
            let u_r = u_initial * ct_yaw_tilt / (2.0 * (1.0 - sqrt_ct_yaw_tilt));
            // Python: u0 = freestream_velocity * np.sqrt(1 - ct_i) (no yaw factor)
            let u0 = u_initial * (1.0 - ct).max(0.0).sqrt();

            let c0 = 1.0 - u0 / u_initial;
            let m0 = c0 * (2.0 - c0);
            let e0 = c0 * c0 - 3.0 * (1.0 / 12.0_f64).exp() * c0 + 3.0 * (1.0 / 3.0_f64).exp();

            // Initial Gaussian wake expansion
            let sigma_z0 = rotor_diameter * 0.5 * (u_r / (u_initial + u0)).sqrt();
            let sigma_y0 =
                sigma_z0 * crate::utilities::cosd(yaw) * crate::utilities::cosd(wind_veer);

            // Yaw parameters (skew angle and distance from centerline)
            let theta_c0 = self.dm * (0.3 * yaw.to_radians() / crate::utilities::cosd(yaw));
            let theta_c0 = theta_c0 * (1.0 - sqrt_ct_yaw_tilt);

            // Start of the near wake (xR == x_i in this model)
            let xr = x_src;

            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        // Per-point TI (Python: turbulence_intensity_i is the
                        // evolving 4D field slice (n_findex, 1, gy, gz),
                        // broadcast over target turbines by index)
                        let ti_pt = turbulence_intensity_i[[fi, 0, iy, iz]];

                        // Length of near wake (Python x0 uses sqrt(1 - ct) in the
                        // denominator and (1 + sqrt(1 - ct*cosd(yaw))) in the numerator)
                        let x0 = rotor_diameter
                            * (crate::utilities::cosd(yaw) * (1.0 + sqrt_ct_yaw_tilt))
                            / (2.0_f64.sqrt()
                                * (4.0 * self.alpha * ti_pt
                                    + 2.0 * self.beta * (1.0 - (1.0 - ct).max(0.0).sqrt())))
                            + x_src;

                        // Wake expansion parameters (per point)
                        let ky = self.ka * ti_pt + self.kb;
                        let kz = self.ka * ti_pt + self.kb;

                        let delta0 = theta_c0.tan() * (x0 - x_src); // initial wake deflection

                        let x = grid_x[[fi, ti, iy, iz]];

                        // Deflection in the near wake:
                        // delta_near = ((x - xR) / (x0 - xR)) * delta0 + (ad + bd*(x - x_i))
                        // masked to xR <= x <= x0
                        let mut delta_near_wake = 0.0;
                        if x >= xr && x <= x0 && (x0 - xr).abs() > 1e-10 {
                            delta_near_wake =
                                ((x - xr) / (x0 - xr)) * delta0 + (self.ad + self.bd * (x - x_src));
                        }

                        // Deflection in the far wake
                        let mut delta_far_wake = 0.0;
                        if x > x0 {
                            let sigma_y = ky * (x - x0) + sigma_y0;
                            let sigma_z = kz * (x - x0) + sigma_z0;

                            let m0_sqrt = m0.max(0.0).sqrt();
                            let middle_term_arg =
                                ((sigma_y * sigma_z) / (sigma_y0 * sigma_z0)).sqrt();
                            let ln_delta_num = (1.6 + m0_sqrt) * (1.6 * middle_term_arg - m0_sqrt);
                            let ln_delta_den = (1.6 - m0_sqrt) * (1.6 * middle_term_arg + m0_sqrt);

                            let mut middle_term = 0.0;
                            if ln_delta_den.abs() > 1e-10
                                && ln_delta_num > 0.0
                                && ln_delta_den > 0.0
                            {
                                let log_val = (ln_delta_num / ln_delta_den).ln();
                                middle_term = theta_c0 * e0 / 5.2
                                    * ((sigma_y0 * sigma_z0) / (ky * kz * m0)).sqrt()
                                    * log_val;
                            }

                            delta_far_wake =
                                delta0 + middle_term + (self.ad + self.bd * (x - x_src));
                        }

                        deflection[[fi, ti, iy, iz]] = delta_near_wake + delta_far_wake;
                    }
                }
            }
        }

        Ok(deflection)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gauss_deflection_creation() {
        let gauss = GaussVelocityDeflection::new(0.01, 0.05, 0.0, 0.58, 0.077, 0.38, 0.004, 1.0);
        assert_eq!(gauss.ad, 0.05);
        assert_eq!(gauss.alpha, 0.58);
        assert_eq!(gauss.beta, 0.077);
        assert_eq!(gauss.dm, 1.0);
        // kd is stored in base.parameters
        assert_eq!(gauss.base.parameters.get("kd"), Some(&0.01));
    }

    #[test]
    fn test_gauss_deflection() {
        // Golden value from Python FLORIS v4.6.4 GaussVelocityDeflection:
        // yaw=+5, ct=0.78, ti=0.1, D=126, kd=0.05, ad=bd=0, x_i=0, evaluated
        // 1 D downstream. Positive yaw deflects toward negative y (Python
        // flips the yaw sign internally).
        let gauss = GaussVelocityDeflection::new(0.05, 0.0, 0.0, 0.58, 0.077, 0.38, 0.004, 1.0);
        let x_i = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        let y_i = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        let yaw = Array1::from_vec(vec![5.0]);
        let ti = Array::from_shape_vec((1, 1, 1, 1), vec![0.1]).unwrap();
        let ct = Array1::from_vec(vec![0.78]);
        // Evaluation point 1 D downstream of the source
        let grid_x = Array::from_shape_vec((1, 1, 1, 1), vec![126.0]).unwrap();

        let result = gauss
            .function(
                &x_i,
                &y_i,
                &yaw,
                &ti,
                &ct,
                &Array2::from_elem((1, 1), 126.0),
                &grid_x,
                &HashMap::new(),
            )
            .unwrap();
        assert_relative_eq!(result[[0, 0, 0, 0]], -1.74781831, epsilon = 1e-6);

        // Zero yaw gives exactly zero deflection
        let yaw0 = Array1::from_vec(vec![0.0]);
        let zero = gauss
            .function(
                &x_i,
                &y_i,
                &yaw0,
                &ti,
                &ct,
                &Array2::from_elem((1, 1), 126.0),
                &grid_x,
                &HashMap::new(),
            )
            .unwrap();
        assert_relative_eq!(zero[[0, 0, 0, 0]], 0.0, epsilon = 1e-12);
    }
}
