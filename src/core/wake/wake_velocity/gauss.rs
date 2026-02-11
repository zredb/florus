use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, GridBase};
/// Gauss wake velocity model
///
/// Gaussian wake model based on Bastankhah and Porte-Agel (2014, 2016)
use crate::types::{Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Gauss wake velocity model parameters
#[derive(Debug, Clone)]
pub struct GaussVelocity {
    pub base: BaseModel,
    pub ka: Float,
    pub kb: Float,
    pub initial_wake_width: Float,
}

impl GaussVelocity {
    pub fn new(ka: Float, kb: Float, initial_wake_width: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("ka".to_string(), ka);
        params.insert("kb".to_string(), kb);
        params.insert("initial_wake_width".to_string(), initial_wake_width);

        Self {
            base: BaseModel::new(params, "wind_vector"),
            ka,
            kb,
            initial_wake_width,
        }
    }
}

impl VelocityModel for GaussVelocity {
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
        _yaw_angle: Float,
        turbulence_intensity: Float,
        _thrust_coefficient: Float,
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

        for fi in 0..n_findex {
            // Get deflection at the wake source turbine
            let deflection_at_source = deflection_field[[fi, turbine_index]];

            // Apply deficit to all grid points that are downstream
            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];

                        // Only apply if this point is downstream of the wake source
                        if x_point <= x_wake_source {
                            continue;
                        }

                        // Calculate distance from wake source
                        let downstream_x = x_point - x_wake_source;

                        // Calculate wake parameters at this downstream position
                        let sigma_y = self.calculate_sigma_y(downstream_x, r0, turbulence_intensity);
                        let sigma_z = self.calculate_sigma_z(downstream_x, r0, turbulence_intensity);

                        let deficit_0 =
                            self.calculate_deficit(axial_induction, downstream_x, r0, turbulence_intensity);

                        let y_point = y[[fi, ti, iy, iz]];
                        let z_point = z[[fi, ti, iy, iz]];

                        // Wake center at this x position
                        let wake_center_y =
                            y_wake_source + deflection_at_source * (x_point - x_wake_source);

                        let dy = y_point - wake_center_y;
                        let dz = z_point - hub_height;

                        // Gaussian wake profile
                        let exponent =
                            -0.5 * (dy.powi(2) / sigma_y.powi(2) + dz.powi(2) / sigma_z.powi(2));
                        let profile = (-exponent).exp();

                        if profile > 1e-10 {
                            velocity_deficit[[fi, ti, iy, iz]] = deficit_0 * profile;
                        }
                    }
                }
            }
        }

        Ok(velocity_deficit)
    }
}

impl GaussVelocity {
    fn calculate_sigma_y(&self, x: Float, r0: Float, _ti: Float) -> Float {
        let epsilon = 0.2 * r0;
        let sigma_y0 = self.initial_wake_width * r0;
        sigma_y0 + self.ka * x + epsilon
    }

    fn calculate_sigma_z(&self, x: Float, r0: Float, _ti: Float) -> Float {
        let epsilon = 0.2 * r0;
        let sigma_z0 = self.initial_wake_width * r0;
        sigma_z0 + self.kb * x + epsilon
    }

    fn calculate_deficit(&self, axial_induction: Float, x: Float, r0: Float, turbulence_intensity: Float) -> Float {
        // At x=0 (right at turbine), deficit is maximum (1 - c1)
        // As x increases, deficit decreases
        if axial_induction <= 0.0 {
            return 0.0;
        }

        let sigma = self.calculate_sigma_y(x, r0, turbulence_intensity);
        let c1 = 1.0 - axial_induction;

        // Bastankhah-Porte-Agel model
        let denominator = 1.0 + 2.0 * sigma.powi(2) / (4.0 * r0.powi(2));
        let deficit = (1.0 - c1) / denominator;

        deficit.max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{core::AveragingMethod, types::Array1};
    // use approx::assert_relative_eq;

    #[test]
    fn test_gauss_velocity_creation() {
        let gauss = GaussVelocity::new(0.1, 0.05, 0.5);
        assert_eq!(gauss.ka, 0.1);
        assert_eq!(gauss.kb, 0.05);
        assert_eq!(gauss.initial_wake_width, 0.5);
    }

    #[test]
    fn test_gauss_velocity_prepare_function() {
        let gauss = GaussVelocity::new(0.1, 0.05, 0.5);
        let result = gauss.prepare_function(&fake_grid(), &fake_flow_field());
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn test_sigma_calculation() {
        let gauss = GaussVelocity::new(0.1, 0.05, 0.5);
        let sigma_y = gauss.calculate_sigma_y(100.0, 63.0, 0.06);
        let sigma_z = gauss.calculate_sigma_z(100.0, 63.0, 0.06);
        assert!(sigma_y > sigma_z);
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
