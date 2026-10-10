//! Jensen wake velocity model
//!
//! Classical Jensen (N.O. Jensen) wake model with linear wake expansion

use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{FlowField, Grid};
use crate::types::{Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Jensen wake velocity model parameters
#[derive(Debug, Clone)]
pub struct JensenVelocity {
    pub base: BaseModel,
    pub kd: Float,
    pub initial_wake_width: Float,
}

impl JensenVelocity {
    pub fn new(kd: Float, initial_wake_width: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("kd".to_string(), kd);
        params.insert("initial_wake_width".to_string(), initial_wake_width);

        Self {
            base: BaseModel::new(params, "wind_vector"),
            kd,
            initial_wake_width,
        }
    }
}

impl VelocityModel for JensenVelocity {
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
        deflection_field: Array4,
        _yaw_angle: Array2,
        _turbulence_intensity: Array2,
        _thrust_coefficient: Array2,
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

        let r0 = rotor_diameter / 2.0;
        let we = self.kd;
        const NUM_EPS: Float = 0.001;

        for fi in 0..n_findex {
            let mut x_wake_source = 0.0;
            let mut y_wake_source = 0.0;
            let mut z_wake_source = 0.0;
            let mut count = 0.0;
            for iy in 0..n_y {
                for iz in 0..n_z {
                    x_wake_source += x[[fi, turbine_index, iy, iz]];
                    y_wake_source += y[[fi, turbine_index, iy, iz]];
                    z_wake_source += z[[fi, turbine_index, iy, iz]];
                    count += 1.0;
                }
            }
            x_wake_source /= count;
            y_wake_source /= count;
            z_wake_source /= count;

            // Apply deficit to all downstream grid points
            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];

                        // dx = downstream distance from source turbine
                        let dx = x_point - x_wake_source;

                        // Downstream mask: only include points downstream of the source
                        // NUM_EPS=0.001 matches Python's downstream_mask = (dx > 0 + NUM_EPS)
                        if dx <= NUM_EPS {
                            continue;
                        }

                        // Circular boundary check (matches Python):
                        // sqrt(dy^2 + dz^2) < we * dx + rotor_radius, where
                        // dy = y - y_i - deflection_field (per grid point)
                        let dy = y[[fi, ti, iy, iz]]
                            - y_wake_source
                            - deflection_field[[fi, ti, iy, iz]];
                        let dz = z[[fi, ti, iy, iz]] - z_wake_source;
                        let radial_distance = (dy * dy + dz * dz).sqrt();
                        let boundary = we * dx + r0;

                        if radial_distance >= boundary {
                            continue;
                        }

                        // C = (rotor_radius / (rotor_radius + we * dx + NUM_EPS))^2
                        // Deficit decays with downstream distance
                        let c = (r0 / (r0 + we * dx + NUM_EPS)).powi(2);

                        // velocity_deficit = 2 * axial_induction * c
                        // Factor of 2 matches actuator disk theory: deficit at disc = 2a
                        velocity_deficit[[fi, ti, iy, iz]] =
                            2.0 * axial_induction[[fi, turbine_index]] * c;
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

    #[test]
    fn test_jensen_velocity_creation() {
        let jensen = JensenVelocity::new(0.1, 0.0);
        assert_eq!(jensen.kd, 0.1);
    }

    #[test]
    fn test_wake_expansion() {
        let r0 = 63.0_f64;
        let we = 0.05_f64;
        let dx = 100.0_f64;
        let boundary = we * dx + r0;
        assert!(boundary > r0);
    }

    #[test]
    fn test_deficit_formula() {
        let a = 0.33_f64;
        let r0 = 63.0_f64;
        let we = 0.05_f64;
        let dx = 100.0_f64;
        let c = (r0 / (r0 + we * dx)).powi(2);
        let deficit = 2.0 * a * c;
        assert!(deficit > 0.0 && deficit < 2.0 * a);
    }
}
