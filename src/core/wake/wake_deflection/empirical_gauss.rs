use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{FlowField, Grid};
/// Empirical Gauss wake deflection model
///
/// Extended Gauss model with empirical corrections for different conditions
use crate::types::{Array1, Array2, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Empirical Gauss wake deflection model parameters
///
/// Simplified formulation independent from the velocity deficit model
/// (port of Python EmpiricalGaussVelocityDeflection).
#[derive(Debug, Clone)]
pub struct EmpiricalGaussVelocityDeflection {
    pub base: BaseModel,
    pub horizontal_deflection_gain_d: Float,
    pub vertical_deflection_gain_d: Float,
    pub deflection_rate: Float,
    pub mixing_gain_deflection: Float,
    pub yaw_added_mixing_gain: Float,
}

impl EmpiricalGaussVelocityDeflection {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        horizontal_deflection_gain_d: Float,
        vertical_deflection_gain_d: Float,
        deflection_rate: Float,
        mixing_gain_deflection: Float,
        yaw_added_mixing_gain: Float,
    ) -> Self {
        let mut params = HashMap::new();
        params.insert(
            "horizontal_deflection_gain_d".to_string(),
            horizontal_deflection_gain_d,
        );
        params.insert("deflection_rate".to_string(), deflection_rate);

        Self {
            base: BaseModel::new(params, "wind_vector"),
            horizontal_deflection_gain_d,
            vertical_deflection_gain_d,
            deflection_rate,
            mixing_gain_deflection,
            yaw_added_mixing_gain,
        }
    }

    pub fn default() -> Self {
        Self::new(3.0, -1.0, 22.0, 0.0, 0.0)
    }
}

impl Default for EmpiricalGaussVelocityDeflection {
    fn default() -> Self {
        Self::default()
    }
}

impl DeflectionModel for EmpiricalGaussVelocityDeflection {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        Ok(HashMap::new())
    }

    fn clone_box(&self) -> Box<dyn DeflectionModel> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn function(
        &self,
        x: Array2,
        _y: Array2,
        yaw_angle: Float,
        _turbulence_intensity: Float,
        thrust_coefficient: Float,
        rotor_diameter: Float,
        _model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array2> {
        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];

        let mut deflection = Array::zeros((n_findex, n_turbines));

        let deflection_gain_y = self.horizontal_deflection_gain_d * rotor_diameter;
        // CW yaw for consistency with other models
        let yaw_r = -yaw_angle.to_radians();

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let x_downstream = x[[fi, ti]];

                if x_downstream <= 0.0 {
                    deflection[[fi, ti]] = 0.0;
                    continue;
                }

                let x_normalized = x_downstream / rotor_diameter;
                // mixing_i unavailable in this interface: mixing = 0
                let a_y = (deflection_gain_y * thrust_coefficient * yaw_r)
                    / (1.0 + self.mixing_gain_deflection * 0.0);
                let log_term = ((x_normalized - self.deflection_rate)
                    / (x_normalized + self.deflection_rate)
                    + 2.0)
                    .ln();
                deflection[[fi, ti]] = a_y * log_term;
            }
        }

        Ok(deflection)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empirical_gauss_deflection_creation() {
        let empirical = EmpiricalGaussVelocityDeflection::new(0.01, 0.05, 0.5);
        assert_eq!(empirical.kd, 0.01);
        assert_eq!(empirical.ad, 0.05);
        assert_eq!(empirical.bd, 0.5);
    }
}
