/// Gauss wake deflection model
///
/// Based on Bastankhah and Porte-Agel (2016) - wake deflection for Gauss velocity model
use crate::types::{Float, Array1, Array2};
use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{Grid, FlowField};
use std::collections::HashMap;
use ndarray::Array;

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
    pub fn new(kd: Float, ad: Float, alpha: Float, beta: Float, dm: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("kd".to_string(), kd);
        params.insert("ad".to_string(), ad);
        params.insert("alpha".to_string(), alpha);
        params.insert("beta".to_string(), beta);
        params.insert("dm".to_string(), dm);
        // Add additional parameters with default values
        params.insert("bd".to_string(), 0.0);
        params.insert("ka".to_string(), 0.38);
        params.insert("kb".to_string(), 0.004);
        params.insert("eps_gain".to_string(), 0.2);
        params.insert("use_secondary_steering".to_string(), 1.0); // Convert bool to float for BaseModel compatibility

        Self {
            base: BaseModel::new(params, "wind_vector"),
            ad,
            bd: 0.0,
            alpha,
            beta,
            ka: 0.38,
            kb: 0.004,
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
        params.insert("wind_veer".to_string(), Array1::from_elem(grid.x_sorted().shape()[0], wind_veer_value));
        
        Ok(params)
    }

    fn function(
        &self,
        x: Array2,
        y: Array2,
        yaw_angle: Float,
        turbulence_intensity: Float,
        thrust_coefficient: Float,
        rotor_diameter: Float,
        model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array2> {
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
        
        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];

        let mut deflection = Array::zeros((n_findex, n_turbines));

        // Opposite sign convention in this model
        let yaw = -yaw_angle;

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let x_i = x[[fi, ti]];
                let y_i = y[[fi, ti]];

                if x_i <= 0.0 {
                    deflection[[fi, ti]] = 0.0;
                    continue;
                }

                // Calculate initial velocity deficits
                let u_initial = 1.0; // normalized
                let ct = thrust_coefficient;
                let tilt = 0.0; // assuming no tilt
                
                let ct_yaw_tilt = ct * crate::utilities::cosd(tilt) * crate::utilities::cosd(yaw);
                let sqrt_ct_yaw_tilt = (1.0 - ct_yaw_tilt).sqrt();
                let uR = u_initial * ct_yaw_tilt 
                    / (2.0 * (1.0 - sqrt_ct_yaw_tilt));
                let u0 = u_initial * sqrt_ct_yaw_tilt;

                // Length of near wake
                let x0 = rotor_diameter * (crate::utilities::cosd(yaw) * (1.0 + sqrt_ct_yaw_tilt))
                    / (2.0_f64.sqrt() * (
                        4.0 * self.alpha * turbulence_intensity + 2.0 * self.beta * (1.0 - sqrt_ct_yaw_tilt)
                    )) + x_i;

                // Wake expansion parameters
                let ky = self.ka * turbulence_intensity + self.kb;
                let kz = self.ka * turbulence_intensity + self.kb;

                let C0 = 1.0 - u0 / u_initial;
                let M0 = C0 * (2.0 - C0);
                let E0 = C0 * C0 - 3.0 * (1.0 / 12.0_f64).exp() * C0 + 3.0 * (1.0 / 3.0_f64).exp();

                // Initial Gaussian wake expansion
                let sigma_z0 = rotor_diameter * 0.5 * (uR / (u_initial + u0)).sqrt();
                let sigma_y0 = sigma_z0 * crate::utilities::cosd(yaw) * crate::utilities::cosd(wind_veer);

                // Yaw parameters (skew angle and distance from centerline)
                let theta_c0 = self.dm * (0.3 * yaw.to_radians() / crate::utilities::cosd(yaw));
                let theta_c0 = theta_c0 * (1.0 - (1.0 - ct * crate::utilities::cosd(yaw)).sqrt());
                let delta0 = theta_c0.tan() * (x0 - x_i); // initial wake deflection

                // Deflection in the near wake
                let xR = x_i; // Reference x location
                let mut delta_near_wake = 0.0;
                if x_i >= xR && x_i <= x0 && (x0 - xR) > 1e-10 {
                    let near_wake_factor = (x_i - xR) / (x0 - xR);
                    delta_near_wake = near_wake_factor * delta0 + (self.ad + self.bd * (x_i - x_i));
                    
                    // Apply mask: only apply near wake deflection in the near wake region
                    if x_i < xR || x_i > x0 {
                        delta_near_wake = 0.0;
                    }
                }

                // Deflection in the far wake
                let mut sigma_y = sigma_y0;
                let mut sigma_z = sigma_z0;
                let downstream_dist = if x_i > x0 { x_i - x0 } else { 0.0 };
                if x_i >= x0 {
                    sigma_y = ky * downstream_dist + sigma_y0;
                    sigma_z = kz * downstream_dist + sigma_z0;
                }

                let mut delta_far_wake = 0.0;
                if x_i > x0 {
                    let M0_sqrt = M0.sqrt();
                    let middle_term_arg = ((sigma_y * sigma_z) / (sigma_y0 * sigma_z0)).sqrt();
                    let ln_delta_num = (1.6 + M0_sqrt) * (1.6 * middle_term_arg - M0_sqrt);
                    let ln_delta_den = (1.6 - M0_sqrt) * (1.6 * middle_term_arg + M0_sqrt);
                    
                    let mut middle_term = 0.0;
                    if ln_delta_den.abs() > 1e-10 && ln_delta_num > 0.0 && ln_delta_den > 0.0 {
                        let log_val = (ln_delta_num / ln_delta_den).ln();
                        middle_term = theta_c0 * E0 / 5.0 
                            * ((sigma_y0 * sigma_z0) / (ky * kz * M0)).sqrt()
                            * log_val;
                    }
                    
                    delta_far_wake = delta0 + middle_term + (self.ad + self.bd * (x_i - x_i));
                }

                deflection[[fi, ti]] = delta_near_wake + delta_far_wake;
            }
        }

        Ok(deflection)
    }
}

impl GaussVelocityDeflection {
    fn calculate_axial_induction(&self, ct: Float) -> Float {
        if ct < 0.96 {
            0.5 * (1.0 - (1.0 - ct).sqrt())
        } else {
            0.143 + (0.0203 - 0.6427 * (0.889 - ct).sqrt()).max(0.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gauss_deflection_creation() {
        let gauss = GaussVelocityDeflection::new(0.01, 0.05, 0.58, 0.077, 1.0);
        assert_eq!(gauss.ad, 0.05);
        assert_eq!(gauss.alpha, 0.58);
        assert_eq!(gauss.beta, 0.077);
        assert_eq!(gauss.dm, 1.0);
        // kd is stored in base.parameters
        assert_eq!(gauss.base.parameters.get("kd"), Some(&0.01));
    }

    #[test]
    fn test_gauss_deflection() {
        let gauss = GaussVelocityDeflection::new(0.01, 0.05, 0.58, 0.077, 1.0);
        let x = Array::from_shape_vec((1, 1), vec![100.0]).unwrap();
        let y = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        
        let result = gauss.function(
            x, y, 25.0, 0.06, 0.33, 126.0,
            &HashMap::new()
        );
        assert!(result.is_ok());
        assert!(result.unwrap()[[0, 0]] > 0.0);
    }
}
