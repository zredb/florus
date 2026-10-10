/// TurbOParkGauss Wake Velocity Model
///
/// TurbOPark model with Gaussian wake profile.
use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, Grid};
use crate::types::{Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct TurbOParkGaussVelocityDeficit {
    pub base: BaseModel,
    pub a: Float,
    pub include_mirror_wake: bool,
}

impl TurbOParkGaussVelocityDeficit {
    pub fn new(a: Float, include_mirror_wake: bool) -> Self {
        let mut params = HashMap::new();
        params.insert("a".to_string(), a);
        params.insert(
            "include_mirror_wake".to_string(),
            if include_mirror_wake { 1.0 } else { 0.0 },
        );

        Self {
            base: BaseModel::new(params, "turboparkgauss"),
            a,
            include_mirror_wake,
        }
    }

    pub fn default() -> Self {
        Self::new(0.04, true)
    }
}

impl Default for TurbOParkGaussVelocityDeficit {
    fn default() -> Self {
        Self::default()
    }
}

#[inline]
fn calculate_epsilon(ct: Float) -> Float {
    let sqrt_one_minus_ct = (1.0 - ct).max(0.0).sqrt();
    if sqrt_one_minus_ct < 1e-10 {
        return 0.25 * (3.0_f64.sqrt() as Float);
    }
    let ratio = 0.5 * (1.0 + sqrt_one_minus_ct) / sqrt_one_minus_ct;
    let epsilon_base = 0.25 * ratio.sqrt();
    epsilon_base.min(0.25 * (3.0_f64.sqrt() as Float))
}

fn turboparkgauss_wake_width(x_dist: Float, ti: Float, ct: Float, a: Float) -> Float {
    let c1 = 1.5;
    let c2 = 0.8;

    let alpha = ti * c1;
    let beta = (c2 * ti / ct.sqrt()).max(1e-10);

    let term1 = ((alpha + beta * x_dist).powi(2) + 1.0).sqrt();
    let term2 = (1.0 + alpha.powi(2)).sqrt();

    let dw = a * ti / beta
        * (term1
            - term2
            - ((term1 + 1.0) * alpha / ((term2 + 1.0) * (alpha + beta * x_dist))).ln());

    dw.max(0.0)
}

impl VelocityModel for TurbOParkGaussVelocityDeficit {
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
        _axial_induction: Array2,
        _deflection_field: Array4,
        _yaw_angle: Array2,
        turbulence_intensity: Array4,
        thrust_coefficient: Array2,
        _hub_height: Float,
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
            let ct_i = thrust_coefficient[[fi, turbine_index]];

            // Compute mean position of the source turbine (matches Python solver.py
            // x_i = np.mean(grid.x_sorted[:, i:i+1], axis=(2, 3), keepdims=True))
            let mut x_source_sum = 0.0;
            let mut y_source_sum = 0.0;
            let mut z_source_sum = 0.0;
            let grid_points = (n_y * n_z) as Float;
            for iy in 0..n_y {
                for iz in 0..n_z {
                    x_source_sum += x[[fi, turbine_index, iy, iz]];
                    y_source_sum += y[[fi, turbine_index, iy, iz]];
                    z_source_sum += z[[fi, turbine_index, iy, iz]];
                }
            }
            let x_source_mean = x_source_sum / grid_points;
            let y_source_mean = y_source_sum / grid_points;
            let z_source_mean = z_source_sum / grid_points;

            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];
                        let dx = x_point - x_source_mean;

                        if dx <= 0.0 {
                            continue;
                        }

                        let y_point = y[[fi, ti, iy, iz]];
                        let z_point = z[[fi, ti, iy, iz]];

                        let x_dist = dx / rotor_diameter;

                        // Per-point TI (source's slice of the evolving 4D field)
                        let ti_i = turbulence_intensity[[fi, turbine_index, iy, iz]];
                        let dw = turboparkgauss_wake_width(x_dist, ti_i, ct_i, self.a);
                        let epsilon = calculate_epsilon(ct_i);
                        let sigma = (epsilon + dw).max(1e-10);
                        let wake_width = sigma * rotor_diameter;

                        if wake_width <= 0.0 {
                            continue;
                        }

                        let ct_clipped = ct_i.clamp(0.0001, 0.9999);
                        let c_val =
                            1.0 - (1.0 - ct_clipped / (8.0 * sigma * sigma)).max(0.0).sqrt();

                        let dy = y_point - y_source_mean;
                        let dz = z_point - z_source_mean;
                        let r_real_sq = dy * dy + dz * dz;
                        let delta_real = if x_dist > 0.0 {
                            c_val * (-r_real_sq / (2.0 * wake_width * wake_width)).exp()
                        } else {
                            0.0
                        };

                        if self.include_mirror_wake {
                            let dz_mirror = z_point - 3.0 * z_source_mean;
                            let r_image_sq = dy * dy + dz_mirror * dz_mirror;
                            let delta_image = if x_dist > 0.0 {
                                c_val * (-r_image_sq / (2.0 * wake_width * wake_width)).exp()
                            } else {
                                0.0
                            };
                            velocity_deficit[[fi, ti, iy, iz]] =
                                (delta_real * delta_real + delta_image * delta_image).sqrt();
                        } else {
                            velocity_deficit[[fi, ti, iy, iz]] = delta_real;
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
    use crate::types::Array1;

    #[test]
    fn test_turboparkgauss_creation() {
        let park = TurbOParkGaussVelocityDeficit::new(0.04, true);
        assert_eq!(park.a, 0.04);
        assert!(park.include_mirror_wake);
    }

    #[test]
    fn test_turboparkgauss_default() {
        let park = TurbOParkGaussVelocityDeficit::default();
        assert_eq!(park.a, 0.04);
        assert!(park.include_mirror_wake);
    }

    #[test]
    fn test_turboparkgauss_prepare_function() {
        let park = TurbOParkGaussVelocityDeficit::default();
        let result = park.prepare_function(&fake_grid(), &fake_flow_field());
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn test_turboparkgauss_wake_width() {
        let width = turboparkgauss_wake_width(5.0, 0.06, 0.8, 0.04);
        assert!(width > 0.0);
    }

    #[test]
    fn test_epsilon() {
        let eps = calculate_epsilon(0.8);
        assert!(eps > 0.0 && eps < 1.0);
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
