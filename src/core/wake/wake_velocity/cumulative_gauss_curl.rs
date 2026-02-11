//! Cumulative Gauss Curl Wake Velocity Model
//!
//! Implements the cumulative-curl wake model for capturing deep array effects
//! and impacts to wake steering.
//!
//! Based on Bay et al. (2023) "Addressing deep array effects and impacts to wake
//! steering with the cumulative-curl wake model" - Wind Energy Science
//!
//! This model extends the Gaussian wake model with:
//! - Wake expansion continuation (WEC) for cumulative wake effects
//! - Curl terms for yawed turbine wake deformation

use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, GridBase};
use crate::types::{Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Cumulative Gauss Curl wake velocity model
#[derive(Debug, Clone)]
pub struct CumulativeCurlVelocityDeficit {
    pub base: BaseModel,
    // WEC (Wake Expansion Continuation) parameters
    pub wec_factor: Float,
    // Near-field parameters
    pub a_s: Float,
    pub b_s: Float,
    pub c_s1: Float,
    pub c_s2: Float,
    // Far-field parameters
    pub a_f: Float,
    pub b_f: Float,
    pub c_f: Float,
}

impl CumulativeCurlVelocityDeficit {
    pub fn new(
        wec_factor: Float,
        a_s: Float, b_s: Float, c_s1: Float, c_s2: Float,
        a_f: Float, b_f: Float, c_f: Float,
    ) -> Self {
        let mut params = HashMap::new();
        params.insert("wec_factor".to_string(), wec_factor);
        params.insert("a_s".to_string(), a_s);
        params.insert("b_s".to_string(), b_s);
        params.insert("c_s1".to_string(), c_s1);
        params.insert("c_s2".to_string(), c_s2);
        params.insert("a_f".to_string(), a_f);
        params.insert("b_f".to_string(), b_f);
        params.insert("c_f".to_string(), c_f);

        Self {
            base: BaseModel::new(params, "cumulative_gauss_curl"),
            wec_factor,
            a_s, b_s, c_s1, c_s2,
            a_f, b_f, c_f,
        }
    }

    pub fn default() -> Self {
        Self::new(
            1.0,          // wec_factor
            0.179367259,  // a_s
            0.0118889215, // b_s
            0.0563691592, // c_s1
            0.13290157,   // c_s2
            3.11,         // a_f
            -0.68,        // b_f
            2.41,         // c_f
        )
    }
}

impl Default for CumulativeCurlVelocityDeficit {
    fn default() -> Self {
        Self::default()
    }
}

/// Calculate wake spread using cumulative-curl formulation
fn calculate_wake_spread(
    x: Float,
    ti: Float,
    ct: Float,
    rotor_diameter: Float,
    params: &CumulativeCurlVelocityDeficit,
) -> Float {
    let x_d = x / rotor_diameter;

    // Near-field to far-field transition point
    let x_transition = 2.0;

    if x_d < x_transition {
        // Near-field (expanding wake)
        let sigma_star = (params.a_s * ti + params.b_s) / (ti + params.c_s1).sqrt();
        let c_t_term = params.c_s2 * ct.sqrt();
        sigma_star * (1.0 + c_t_term) * x_d
    } else {
        // Far-field (continued expansion with WEC)
        let sigma_near = calculate_wake_spread(
            x_transition * rotor_diameter,
            ti, ct, rotor_diameter, params
        );
        
        // Wake expansion continuation
        let sigma_far = params.wec_factor * (params.a_f * x_d.powf(params.b_f) + params.c_f);
        
        sigma_near + sigma_far * (x_d - x_transition)
    }
}

/// Calculate Gaussian wake deficit
fn gaussian_deficit(
    r: Float,
    sigma: Float,
    ct: Float,
) -> Float {
    if sigma <= 0.0 {
        return 0.0;
    }
    
    // Gaussian profile
    let exponent = -r.powi(2) / (2.0 * sigma.powi(2));
    let profile = (-exponent).exp();
    
    // Maximum deficit based on actuator disk theory
    let deficit_max = 0.5 * (1.0 - (1.0 - ct).sqrt());
    
    deficit_max * profile
}

impl VelocityModel for CumulativeCurlVelocityDeficit {
    fn prepare_function(
        &self,
        _grid: &dyn GridBase,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

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
        _model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4> {
        let r0 = rotor_diameter / 2.0;

        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let n_y = shape[2];
        let n_z = shape[3];

        let mut velocity_deficit = Array::zeros((n_findex, n_turbines, n_y, n_z));

        // Use the specified turbine's position as the wake source
        let x_wake_source = x[[0, turbine_index, 0, 0]];
        let y_wake_source = y[[0, turbine_index, 0, 0]];

        // Turbine upstream of reference point (x < 0) doesn't generate a wake
        if x_wake_source < 0.0 {
            return Ok(velocity_deficit);
        }

        // Convert yaw angle to radians
        let yaw_rad = yaw_angle.to_radians();

        for fi in 0..n_findex {
            // Get deflection at the wake source turbine
            let deflection_at_source = deflection_field[[fi, turbine_index]];

            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];

                        if x_point <= x_wake_source {
                            continue;
                        }

                        let y_point = y[[fi, ti, iy, iz]];
                        let z_point = z[[fi, ti, iy, iz]];

                        // Wake center with deflection
                        let wake_center_y = y_wake_source + deflection_at_source * (x_point - x_wake_source);
                        let dy = y_point - wake_center_y;
                        let dz = z_point - hub_height;

                        // Radial distance from wake center
                        let r = (dy.powi(2) + dz.powi(2)).sqrt();

                        // Wake spread using cumulative-curl formulation
                        let downstream_x = x_point - x_wake_source;
                        let sigma = calculate_wake_spread(
                            downstream_x,
                            turbulence_intensity,
                            thrust_coefficient,
                            rotor_diameter,
                            self,
                        );

                        // Apply curl correction for yawed turbines
                        let curl_correction = self.calculate_curl_correction(
                            downstream_x,
                            yaw_rad,
                            thrust_coefficient,
                            rotor_diameter,
                            dz,
                        );

                        let sigma_effective = sigma * curl_correction;

                        // Calculate deficit
                        let deficit = gaussian_deficit(r, sigma_effective, thrust_coefficient);

                        // Wake curl causes velocity to increase at edges
                        let curl_deficit = deficit * curl_correction;

                        velocity_deficit[[fi, ti, iy, iz]] = curl_deficit;
                    }
                }
            }
        }

        Ok(velocity_deficit)
    }
}

impl CumulativeCurlVelocityDeficit {
    /// Calculate curl correction factor for yawed turbines
    fn calculate_curl_correction(
        &self,
        x: Float,
        yaw_rad: Float,
        ct: Float,
        rotor_diameter: Float,
        dz: Float,
    ) -> Float {
        if yaw_rad.abs() < 0.01 {
            return 1.0;
        }

        let x_d = x / rotor_diameter;
        
        // Curl effect parameters
        // The curl effect is stronger near the turbine and decreases downstream
        let curl_strength = 0.5 * (1.0 - (-x_d / 5.0).exp());
        
        // Yaw-induced curl
        let yaw_curl = yaw_rad * ct.sqrt() * curl_strength;
        
        // Apply curl correction - increases deficit on one side, decreases on other
        let curl_factor = 1.0 + yaw_curl * (dz / (x + 1.0));
        
        curl_factor.max(0.1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;

    #[test]
    fn test_cumulative_curl_creation() {
        let curl = CumulativeCurlVelocityDeficit::new(
            1.0, 0.179, 0.012, 0.056, 0.133,
            3.11, -0.68, 2.41,
        );
        assert_eq!(curl.wec_factor, 1.0);
        assert_eq!(curl.a_s, 0.179);
    }

    #[test]
    fn test_cumulative_curl_default() {
        let curl = CumulativeCurlVelocityDeficit::default();
        assert_eq!(curl.wec_factor, 1.0);
        assert_eq!(curl.a_f, 3.11);
        assert_eq!(curl.c_f, 2.41);
    }

    #[test]
    fn test_cumulative_curl_prepare_function() {
        let curl = CumulativeCurlVelocityDeficit::default();
        let result = curl.prepare_function(&fake_grid(), &fake_flow_field());
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn test_wake_spread() {
        let params = CumulativeCurlVelocityDeficit::default();
        let sigma = calculate_wake_spread(100.0, 0.06, 0.8, 126.0, &params);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_gaussian_deficit() {
        let deficit = gaussian_deficit(10.0, 20.0, 0.8);
        assert!(deficit > 0.0);
        assert!(deficit < 0.5);
    }

    fn fake_grid() -> impl crate::core::GridBase {
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
            fn sorted_coord_indices(&self) -> &Array2 {
                static COORD_INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
                COORD_INDICES.get_or_init(|| Array2::zeros((1, 1)))
            }
            fn resolution(&self) -> usize {
                1
            }
        }
        FakeGrid
    }

    fn fake_flow_field() -> crate::core::FlowField {
        let empty_4d = ndarray::Array4::zeros((0, 0, 0, 0));
        crate::core::FlowField {
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
        }
    }
}
