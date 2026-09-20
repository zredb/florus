//! Empirical Gauss Wake Velocity Model
//!
//! Extended Gauss model with empirical corrections for different conditions.
//! This is a more tunable version of the Gaussian wake model.

use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, Grid};
use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::cosd;
use ndarray::Array;
use std::collections::HashMap;

/// Sigmoid integral for the smoothed piecewise-linear wake expansion
/// (port of Python sigmoid_integral)
fn sigmoid_integral(x: Float, center: Float, width: Float) -> Float {
    if (x - center) > width / 2.0 {
        x - center
    } else if (x - center) >= -width / 2.0 {
        let z = (x - center) / width + 0.5;
        width * (z.powi(6) - 3.0 * z.powi(5) + 2.5 * z.powi(4))
    } else {
        0.0
    }
}

/// Piecewise-linear wake expansion with smoothing
/// (port of Python empirical_gauss_model_wake_width)
fn empirical_gauss_wake_width(
    x: Float,
    wake_expansion_rates: &[Float],
    breakpoints: &[Float],
    sigma_0: Float,
    smoothing_length: Float,
    mixing_final: Float,
) -> Float {
    let mut sigma = (wake_expansion_rates[0] + mixing_final) * x + sigma_0;
    for (ib, &b) in breakpoints.iter().enumerate() {
        sigma += (wake_expansion_rates[ib + 1] - wake_expansion_rates[ib])
            * sigmoid_integral(x, b, smoothing_length);
    }
    sigma
}

/// r and C for the empirical gauss (port of Python rCalt).
/// dy = y - y_i - delta_y; dz = z - HH - delta_z
fn r_c_alt(
    wind_veer: Float,
    sigma_y: Float,
    sigma_z: Float,
    dy: Float,
    dz: Float,
    ct: Float,
    yaw: Float,
    tilt: Float,
    sigma_y0: Float,
    sigma_z0: Float,
) -> (Float, Float) {
    let veer_rad = wind_veer.to_radians();
    let cos_v = veer_rad.cos();
    let sin_v = veer_rad.sin();
    let sin_2v = (2.0 * veer_rad).sin();

    let a = cos_v * cos_v / (2.0 * sigma_y * sigma_y) + sin_v * sin_v / (2.0 * sigma_z * sigma_z);
    let b = -sin_2v / (4.0 * sigma_y * sigma_y) + sin_2v / (4.0 * sigma_z * sigma_z);
    let c = sin_v * sin_v / (2.0 * sigma_y * sigma_y) + cos_v * cos_v / (2.0 * sigma_z * sigma_z);

    let r = a * (dy * dy) - 2.0 * b * dy * dz + c * (dz * dz);
    let d = 1.0 - ct * (sigma_y0 * sigma_z0) / (sigma_y * sigma_z) * cosd(yaw) * cosd(tilt);
    let c_val = 1.0 - d.max(0.0).sqrt();
    (r, c_val)
}

/// Empirical Gauss wake velocity model
///
/// Gaussian profile expanding in a (smoothed) piecewise linear fashion,
/// with wake-induced mixing gain (port of Python EmpiricalGaussVelocityDeficit).
#[derive(Debug, Clone)]
pub struct EmpiricalGaussVelocityDeficit {
    pub base: BaseModel,
    pub wake_expansion_rates: Vec<Float>,
    pub breakpoints_d: Vec<Float>,
    pub sigma_0_d: Float,
    pub smoothing_length_d: Float,
    pub mixing_gain_velocity: Float,
    pub awc_wake_exp: Float,
    pub awc_wake_denominator: Float,
}

impl EmpiricalGaussVelocityDeficit {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        wake_expansion_rates: Vec<Float>,
        breakpoints_d: Vec<Float>,
        sigma_0_d: Float,
        smoothing_length_d: Float,
        mixing_gain_velocity: Float,
        awc_wake_exp: Float,
        awc_wake_denominator: Float,
    ) -> Self {
        let mut params = HashMap::new();
        params.insert("sigma_0_d".to_string(), sigma_0_d);
        params.insert("smoothing_length_d".to_string(), smoothing_length_d);
        params.insert("mixing_gain_velocity".to_string(), mixing_gain_velocity);

        Self {
            base: BaseModel::new(params, "empirical_gauss"),
            wake_expansion_rates,
            breakpoints_d,
            sigma_0_d,
            smoothing_length_d,
            mixing_gain_velocity,
            awc_wake_exp,
            awc_wake_denominator,
        }
    }

    pub fn default() -> Self {
        Self::new(vec![0.023, 0.008], vec![10.0], 0.28, 2.0, 2.0, 1.2, 400.0)
    }
}

impl Default for EmpiricalGaussVelocityDeficit {
    fn default() -> Self {
        Self::default()
    }
}

impl VelocityModel for EmpiricalGaussVelocityDeficit {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

    fn clone_box(&self) -> Box<dyn VelocityModel> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn function(
        &self,
        x: Array4,
        y: Array4,
        z: Array4,
        axial_induction: Array2,
        deflection_field: Array2,
        yaw_angle: Array2,
        turbulence_intensity: Array2,
        thrust_coefficient: Array2,
        hub_height: Float,
        rotor_diameter: Float,
        turbine_index: usize,
        _model_args: &HashMap<String, Array4>,
    ) -> anyhow::Result<Array4> {
        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let n_y = shape[2];
        let n_z = shape[3];

        let mut velocity_deficit = Array::zeros((n_findex, n_turbines, n_y, n_z));

        for fi in 0..n_findex {
            // Per-findex wake source position (rotated coordinates differ per findex)
            let x_i = x[[fi, turbine_index, 0, 0]];
            let y_i = y[[fi, turbine_index, 0, 0]];

            if x_i < 0.0 {
                continue;
            }

            let deflection_y_i = deflection_field[[fi, turbine_index]];
            let yaw_i = yaw_angle[[fi, turbine_index]];
            let ct_i = thrust_coefficient[[fi, turbine_index]];

            let sigma_y0 = self.sigma_0_d * rotor_diameter * cosd(yaw_i);
            let sigma_z0 = self.sigma_0_d * rotor_diameter;
            let breakpoints: Vec<Float> = self
                .breakpoints_d
                .iter()
                .map(|&b| b * rotor_diameter)
                .collect();
            let smoothing_length = self.smoothing_length_d * rotor_diameter;

            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];

                        let downstream_mask = x_point > x_i + 0.1;
                        if !downstream_mask {
                            continue;
                        }

                        let y_point = y[[fi, ti, iy, iz]];
                        let z_point = z[[fi, ti, iy, iz]];

                        // Wake width without mixing (mixing_i unavailable here)
                        let sigma_y = empirical_gauss_wake_width(
                            x_point - x_i,
                            &self.wake_expansion_rates,
                            &breakpoints,
                            sigma_y0,
                            smoothing_length,
                            0.0,
                        )
                        .max(sigma_y0);
                        let sigma_z = empirical_gauss_wake_width(
                            x_point - x_i,
                            &self.wake_expansion_rates,
                            &breakpoints,
                            sigma_z0,
                            smoothing_length,
                            0.0,
                        )
                        .max(sigma_z0);

                        let dy = y_point - y_i - deflection_y_i;
                        let dz = z_point - hub_height;

                        let (r, c_raw) = r_c_alt(
                            0.0, sigma_y, sigma_z, dy, dz, ct_i, yaw_i, 0.0, sigma_y0, sigma_z0,
                        );
                        let c_val = c_raw / (8.0 * self.sigma_0_d * self.sigma_0_d);

                        let deficit_real = c_val * (-r).exp();

                        // Mirror wake (ground reflection)
                        let r_mirr = ((y_point - y_i - deflection_y_i).powi(2)
                            + (z_point + hub_height).powi(2))
                        .sqrt();
                        let (r_m, c_mirr_raw) = r_c_alt(
                            0.0,
                            sigma_y,
                            sigma_z,
                            y_point - y_i - deflection_y_i,
                            z_point + hub_height,
                            ct_i,
                            yaw_i,
                            0.0,
                            sigma_y0,
                            sigma_z0,
                        );
                        let _ = r_mirr;
                        let c_mirr = c_mirr_raw / (8.0 * self.sigma_0_d * self.sigma_0_d);
                        let deficit_mirror = c_mirr * (-r_m).exp();

                        let deficit =
                            (deficit_real * deficit_real + deficit_mirror * deficit_mirror).sqrt();

                        velocity_deficit[[fi, ti, iy, iz]] = deficit;
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
    use crate::types::Array1;

    #[test]
    fn test_empirical_gauss_creation() {
        let gauss = EmpiricalGaussVelocityDeficit::new(0.1, 0.05, 0.5);
        assert_eq!(gauss.ka, 0.1);
        assert_eq!(gauss.kb, 0.05);
        assert_eq!(gauss.initial_wake_width, 0.5);
    }

    #[test]
    fn test_empirical_gauss_default() {
        let gauss = EmpiricalGaussVelocityDeficit::default();
        assert_eq!(gauss.ka, 0.1);
        assert_eq!(gauss.kb, 0.05);
        assert_eq!(gauss.initial_wake_width, 0.5);
    }

    #[test]
    fn test_empirical_gauss_prepare_function() {
        let gauss = EmpiricalGaussVelocityDeficit::default();
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
            dudz_initial_sorted: empty_4d.clone(),
            turbulence_intensity_field_sorted: empty_4d,
        }
    }
}
