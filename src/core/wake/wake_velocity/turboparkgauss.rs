/// TurbOParkGauss Wake Velocity Model
///
/// TurbOPark model with Gaussian wake profile.

use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, GridBase};
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
        params.insert("include_mirror_wake".to_string(), if include_mirror_wake { 1.0 } else { 0.0 });

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
        * (term1 - term2 - ((term1 + 1.0) * alpha / ((term2 + 1.0) * (alpha + beta * x_dist))).ln());

    dw.max(0.0)
}

impl VelocityModel for TurbOParkGaussVelocityDeficit {
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
        thrust_coefficient: Float,
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

        let r0 = rotor_diameter / 2.0;
        let mut velocity_deficit = Array::zeros((n_findex, n_turbines, n_y, n_z));

        // Use the specified turbine's position as the wake source
        let x_wake_source = x[[0, turbine_index, 0, 0]];
        let y_wake_source = y[[0, turbine_index, 0, 0]];

        // Turbine upstream of reference point (x < 0) doesn't generate a wake
        if x_wake_source < 0.0 {
            return Ok(velocity_deficit);
        }

        let epsilon = calculate_epsilon(thrust_coefficient);

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

                        let wake_center_y = y_wake_source + deflection_at_source * (x_point - x_wake_source);
                        let dy = (y_point - wake_center_y).abs();
                        let dz = z_point - hub_height;
                        let r = (dy.powi(2) + dz.powi(2)).sqrt();

                        let x_dist = (x_point - x_wake_source) / rotor_diameter;
                        let dw = turboparkgauss_wake_width(x_dist, turbulence_intensity, thrust_coefficient, self.a);
                        let sigma = (self.a * turbulence_intensity + epsilon + dw).max(epsilon);
                        let wake_width = sigma * rotor_diameter;

                        if wake_width > 0.0 {
                            let exponent = -r.powi(2) / (2.0 * wake_width.powi(2));
                            let profile = (-exponent).exp();

                            let deficit_max = 2.0 * axial_induction / (1.0 + axial_induction);
                            velocity_deficit[[fi, ti, iy, iz]] = deficit_max * profile;

                            if self.include_mirror_wake {
                                let z_ground = -hub_height;
                                let dz_mirror = z_point - z_ground;
                                let r_mirror = (dy.powi(2) + dz_mirror.powi(2)).sqrt();
                                let exponent_mirror = -r_mirror.powi(2) / (2.0 * wake_width.powi(2));
                                let profile_mirror = (-exponent_mirror).exp();
                                velocity_deficit[[fi, ti, iy, iz]] += deficit_max * 0.5 * profile_mirror;
                            }
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

    fn fake_grid() -> impl crate::core::GridBase {
        struct FakeGrid;
        impl crate::core::GridBase for FakeGrid {
            fn n_turbines(&self) -> usize { 1 }
            fn n_findex(&self) -> usize { 1 }
            fn x_sorted(&self) -> &Array4 { panic!() }
            fn y_sorted(&self) -> &Array4 { panic!() }
            fn z_sorted(&self) -> &Array4 { panic!() }
            fn x_sorted_inertial_frame(&self) -> &Array4 { panic!() }
            fn y_sorted_inertial_frame(&self) -> &Array4 { panic!() }
            fn z_sorted_inertial_frame(&self) -> &Array4 { panic!() }
            fn cubature_weights(&self) -> Option<&Array2> { None }
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
            fn resolution(&self) -> usize { 1 }
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
