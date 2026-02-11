//! Jensen wake velocity model
//!
//! Classical Jensen (N.O. Jensen) wake model with linear wake expansion

use crate::types::{Float, Array2, Array4};
use crate::core::wake::{BaseModel, VelocityModel};
use crate::core::{GridBase, FlowField};
use std::collections::HashMap;
use ndarray::Array;

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
        _grid: &dyn GridBase,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array4>> {
        Ok(HashMap::new())
    }

    fn function(
        &self,
        x: Array4,
        y: Array4,
        _z: Array4,
        axial_induction: Float,
        deflection_field: Array2,
        _yaw_angle: Float,
        turbulence_intensity: Float,
        _thrust_coefficient: Float,
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
            
            // Calculate wake radius at the source
            let r0 = rotor_diameter / 2.0;

            // Calculate base deficit at source
            let deficit_0 = self.calculate_deficit(axial_induction, x_wake_source, r0, turbulence_intensity);

            // Apply deficit to all downstream grid points
            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let x_point = x[[fi, ti, iy, iz]];
                        
                        // Only apply if this point is downstream
                        if x_point <= x_wake_source {
                            continue;
                        }

                        // Calculate wake radius at this downstream position
                        let dx = x_point - x_wake_source;
                        let wake_radius = r0 + self.kd * dx + self.initial_wake_width * r0;
                        
                        // Calculate distance from wake center
                        let y_point = y[[fi, ti, iy, iz]];
                        let wake_center_y = y_wake_source + deflection_at_source * dx;
                        let lateral_distance = (y_point - wake_center_y).abs();

                        // Jensen model: uniform deficit within wake radius
                        if lateral_distance < wake_radius {
                            velocity_deficit[[fi, ti, iy, iz]] = deficit_0;
                        }
                    }
                }
            }
        }

        Ok(velocity_deficit)
    }
}

impl JensenVelocity {
    fn calculate_wake_radius(&self, x: Float, r0: Float, _ti: Float) -> Float {
        r0 + self.kd * x + self.initial_wake_width * r0
    }

    fn calculate_deficit(
        &self,
        axial_induction: Float,
        x: Float,
        r0: Float,
        ti: Float,
    ) -> Float {
        // At x=0 (right at turbine), maximum deficit
        // Deficit decreases as wake expands
        if axial_induction <= 0.0 {
            return 0.0;
        }

        let wake_radius = self.calculate_wake_radius(x, r0, ti);
        let area_ratio = (r0 / wake_radius).powi(2);
        
        let deficit = axial_induction * area_ratio;
        
        deficit.max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    //use approx::assert_relative_eq;

    #[test]
    fn test_jensen_velocity_creation() {
        let jensen = JensenVelocity::new(0.1, 0.0);
        assert_eq!(jensen.kd, 0.1);
    }

    #[test]
    fn test_wake_radius() {
        let jensen = JensenVelocity::new(0.1, 0.0);
        let radius = jensen.calculate_wake_radius(100.0, 63.0, 0.06);
        assert!(radius > 63.0);
    }

    #[test]
    fn test_jensen_deficit() {
        let jensen = JensenVelocity::new(0.1, 0.0);
        let deficit = jensen.calculate_deficit(0.33, 100.0, 63.0, 0.06);
        assert!(deficit > 0.0 && deficit < 0.33);
    }
}
