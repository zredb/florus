use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, Grid};
/// Gauss wake velocity model
///
/// Gaussian wake model based on Bastankhah and Porte-Agel (2014, 2016)
/// Matches Python FLORIS implementation
use crate::types::{Array2, Array4, Float};
use crate::utilities::{cosd, sind};
use ndarray::Array;
use std::collections::HashMap;

/// Gauss wake velocity model parameters
#[derive(Debug, Clone)]
pub struct GaussVelocity {
    pub base: BaseModel,
    pub alpha: Float,
    pub beta: Float,
    pub ka: Float,
    pub kb: Float,
}

impl GaussVelocity {
    pub fn new(alpha: Float, beta: Float, ka: Float, kb: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("alpha".to_string(), alpha);
        params.insert("beta".to_string(), beta);
        params.insert("ka".to_string(), ka);
        params.insert("kb".to_string(), kb);

        Self {
            base: BaseModel::new(params, "wind_vector"),
            alpha,
            beta,
            ka,
            kb,
        }
    }
}

impl VelocityModel for GaussVelocity {
    fn prepare_function(
        &self,
        grid: &dyn Grid,
        flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        let mut params = HashMap::new();
        
        // Add wind_veer to parameters
        let wind_veer_value = flow_field.wind_veer;
        let wind_veer_array = Array4::from_elem(grid.x_sorted().dim(), wind_veer_value);
        params.insert("wind_veer".to_string(), wind_veer_array);
        
        Ok(params)
    }

    fn function(
        &self,
        x: Array4,
        y: Array4,
        z: Array4,
        _axial_induction: Float,
        deflection_field: Array2,
        yaw_angle: Float,
        turbulence_intensity: Float,
        thrust_coefficient: Float,
        hub_height: Float,
        rotor_diameter: Float,
        turbine_index: usize,
        model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4> {
        let r0 = rotor_diameter / 2.0;
        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let n_y = shape[2];
        let n_z = shape[3];

        let mut velocity_deficit = Array::zeros((n_findex, n_turbines, n_y, n_z));

        // Get wind veer from model_args
        let wind_veer = if let Some(wind_veer_arr) = model_args.get("wind_veer") {
            wind_veer_arr[[0, 0, 0, 0]] // Assuming it's constant across all dimensions
        } else {
            0.0 // Default value
        };

        // Calculate turbine center coordinates as mean of all grid points
        // This matches Python FLORIS behavior: np.mean(grid.x_sorted[:, i:i+1], axis=(2, 3))
        let mut x_i_array = Array::zeros((n_findex, n_turbines));
        let mut y_i_array = Array::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let mut x_sum = 0.0;
                let mut y_sum = 0.0;
                let mut count = 0;
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        x_sum += x[[fi, ti, iy, iz]];
                        y_sum += y[[fi, ti, iy, iz]];
                        count += 1;
                    }
                }
                x_i_array[[fi, ti]] = x_sum / count as f64;
                y_i_array[[fi, ti]] = y_sum / count as f64;
            }
        }

        let x_i = x_i_array[[0, turbine_index]];
        let y_i = y_i_array[[0, turbine_index]];

        // Turbine upstream of reference point (x < 0) doesn't generate a wake
        if x_i < 0.0 {
            return Ok(velocity_deficit);
        }

        // Opposite sign convention in this model
        let yaw = -yaw_angle;

        let ct = thrust_coefficient;
        if ct <= 0.0 {
            return Ok(velocity_deficit);
        }

        let sqrt_one_minus_ct = (1.0 - ct).sqrt();

        // Calculate initial wake parameters
        // uR and u0 are used to determine initial wake width
        let u_initial = 1.0; // Normalized
        let uR = u_initial * ct / (2.0 * (1.0 - sqrt_one_minus_ct));
        let u0 = u_initial * sqrt_one_minus_ct;

        // Initial wake width (sigma at x0)
        let sigma_z0 = rotor_diameter * 0.5 * (uR / (u_initial + u0)).sqrt();
        let sigma_y0 = sigma_z0 * cosd(yaw) * cosd(wind_veer);  // Apply wind veer correction

        // Start of far wake (x0)
        let x0 = rotor_diameter * cosd(yaw) * (1.0 + sqrt_one_minus_ct)
            / (2.0_f64.sqrt()
                * (4.0 * self.alpha * turbulence_intensity
                    + 2.0 * self.beta * (1.0 - sqrt_one_minus_ct)));
        let x0 = x0 + x_i;  // Add turbine location to x0

        // Wake expansion rate
        let ky = self.ka * turbulence_intensity + self.kb;
        let kz = self.ka * turbulence_intensity + self.kb;

        for fi in 0..n_findex {
            let deflection_at_source = deflection_field[[fi, turbine_index]];

            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];

                        // Only apply if this point is downstream of the wake source
                        if x_point <= x_i + 0.1 {
                            continue;
                        }

                        let downstream_dist = x_point - x_i;

                        // Determine if in near wake or far wake
                        let (sigma_y, sigma_z) = if x_point < x_i + x0 {
                            // Near wake region - linear interpolation
                            let near_wake_ramp_up = (x_point - x_i) / (x0 - x_i);
                            let near_wake_ramp_down = (x0 - x_point) / (x0 - x_i);

                            let sy_base = 0.501 * rotor_diameter * (ct / 2.0).sqrt();
                            let sz_base = 0.501 * rotor_diameter * (ct / 2.0).sqrt();
                            
                            let sy = near_wake_ramp_down * sy_base + near_wake_ramp_up * sigma_y0;
                            let sy = sy * ((x_point >= x_i) as i32 as f64) + sy_base * ((x_point < x_i) as i32 as f64);
                            
                            let sz = near_wake_ramp_down * sz_base + near_wake_ramp_up * sigma_z0;
                            let sz = sz * ((x_point >= x_i) as i32 as f64) + sz_base * ((x_point < x_i) as i32 as f64);
                            
                            (sy, sz)
                        } else {
                            // Far wake region
                            let sy = ky * (downstream_dist - x0) + sigma_y0;
                            let sz = kz * (downstream_dist - x0) + sigma_z0;
                            (sy, sz)
                        };

                        let y_point = y[[fi, ti, iy, iz]];
                        let z_point = z[[fi, ti, iy, iz]];

                        // Wake center at this x position (with deflection)
                        let wake_center_y = y_i + deflection_at_source;

                        // Calculate Gaussian wake profile
                        let dy = y_point - wake_center_y;
                        let dz = z_point - hub_height;

                        // Elliptic Gaussian profile with wind veer consideration
                        let wind_veer_rad = wind_veer.to_radians();
                        let cos_veer_sq = cosd(wind_veer).powi(2);
                        let sin_veer_sq = sind(wind_veer).powi(2);
                        let sin_2veer = sind(2.0 * wind_veer);
                        
                        let a = (cos_veer_sq / (2.0 * sigma_y * sigma_y)) + (sin_veer_sq / (2.0 * sigma_z * sigma_z));
                        let b = (-sin_2veer / (4.0 * sigma_y * sigma_y)) + (sin_2veer / (4.0 * sigma_z * sigma_z));
                        let c_coeff = (sin_veer_sq / (2.0 * sigma_y * sigma_y)) + (cos_veer_sq / (2.0 * sigma_z * sigma_z));
                        
                        let dy_offset = dy;
                        let dz_offset = dz;
                        
                        let r_squared = a * (dy_offset * dy_offset)
                            - 2.0 * b * dy_offset * dz_offset
                            + c_coeff * (dz_offset * dz_offset);
                        
                        // C coefficient for velocity deficit
                        let d = 1.0 - (ct * cosd(yaw) / (8.0 * sigma_y * sigma_z / (rotor_diameter * rotor_diameter)));
                        let d_clipped = d.max(0.0).min(1.0);  // Use max/min instead of clamp for compatibility
                        let c = 1.0 - d_clipped.sqrt();
                        
                        // Gaussian function
                        let deficit = c * (-r_squared).exp();

                        if deficit > 1e-10 {
                            velocity_deficit[[fi, ti, iy, iz]] = deficit;
                        }
                    }
                }
            }
        }

        Ok(velocity_deficit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{core::AveragingMethod, types::Array1};

    #[test]
    fn test_gauss_velocity_creation() {
        let gauss = GaussVelocity::new(0.58, 0.077, 0.38, 0.004);
        assert_eq!(gauss.alpha, 0.58);
        assert_eq!(gauss.beta, 0.077);
        assert_eq!(gauss.ka, 0.38);
        assert_eq!(gauss.kb, 0.004);
    }

    #[test]
    fn test_gauss_velocity_prepare_function() {
        let gauss = GaussVelocity::new(0.58, 0.077, 0.38, 0.004);
        let result = gauss.prepare_function(&fake_grid(), &fake_flow_field());
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    fn fake_grid() -> impl crate::core::Grid {
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
            fn average_method(&self) -> AveragingMethod {
                AveragingMethod::CubicMean
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
