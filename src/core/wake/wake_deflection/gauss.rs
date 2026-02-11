/// Gauss wake deflection model
///
/// Based on Bastankhah and Porte-Agel (2016) - wake deflection for Gauss velocity model
use crate::types::{Float, Array1, Array2};
use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{GridBase, FlowField};
use std::collections::HashMap;
use ndarray::Array;

/// Gauss wake deflection model parameters
#[derive(Debug, Clone)]
pub struct GaussVelocityDeflection {
    pub base: BaseModel,
    pub kd: Float,
    pub ad: Float,
    pub alpha: Float,
    pub beta: Float,
    pub dm: Float,
}

impl GaussVelocityDeflection {
    pub fn new(kd: Float, ad: Float, alpha: Float, beta: Float, dm: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("kd".to_string(), kd);
        params.insert("ad".to_string(), ad);
        params.insert("alpha".to_string(), alpha);
        params.insert("beta".to_string(), beta);
        params.insert("dm".to_string(), dm);

        Self {
            base: BaseModel::new(params, "wind_vector"),
            kd,
            ad,
            alpha,
            beta,
            dm,
        }
    }
}

impl DeflectionModel for GaussVelocityDeflection {
    fn prepare_function(
        &self,
        _grid: &dyn GridBase,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        Ok(HashMap::new())
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
        let yaw_rad = yaw_angle.to_radians();
        
        let shape = x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];

        let mut deflection = Array::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let x_downstream = x[[fi, ti]];
                
                if x_downstream <= 0.0 {
                    deflection[[fi, ti]] = 0.0;
                    continue;
                }

                let axial_induction = self.calculate_axial_induction(thrust_coefficient);
                
                let c = self.ad / (self.kd * 2.0);
                let exp_term = (-self.kd * x_downstream / rotor_diameter).exp();
                
                deflection[[fi, ti]] = c * axial_induction * (1.0 - exp_term) * yaw_rad;
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
        assert_eq!(gauss.kd, 0.01);
        assert_eq!(gauss.ad, 0.05);
        assert_eq!(gauss.alpha, 0.58);
        assert_eq!(gauss.beta, 0.077);
        assert_eq!(gauss.dm, 1.0);
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
