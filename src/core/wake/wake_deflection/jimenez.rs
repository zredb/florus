/// Jimenez wake deflection model
///
/// Based on Jimenez et al. (2010) - wake deflection due to yaw misalignment
use crate::types::{Float, Array1, Array2};
use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{GridBase, FlowField};
use std::collections::HashMap;
use ndarray::Array;

/// Jimenez wake deflection model parameters
#[derive(Debug, Clone)]
pub struct JimenezVelocityDeflection {
    pub base: BaseModel,
    pub kd: Float,
    pub ad: Float,
}

impl JimenezVelocityDeflection {
    pub fn new(kd: Float, ad: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("kd".to_string(), kd);
        params.insert("ad".to_string(), ad);
        
        Self {
            base: BaseModel::new(params, "wind_vector"),
            kd,
            ad,
        }
    }
}

impl DeflectionModel for JimenezVelocityDeflection {
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
        _thrust_coefficient: Float,
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
                
                if x_downstream > 0.0 {
                    let exp_term = (-self.kd * x_downstream / rotor_diameter).exp();
                    deflection[[fi, ti]] = (self.ad / self.kd) * (1.0 - exp_term) * yaw_rad;
                } else {
                    deflection[[fi, ti]] = 0.0;
                }
            }
        }

        Ok(deflection)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_jimenez_deflection_creation() {
        let jimenez = JimenezVelocityDeflection::new(0.01, 0.05);
        assert_eq!(jimenez.kd, 0.01);
        assert_eq!(jimenez.ad, 0.05);
    }

    #[test]
    fn test_jimenez_deflection_at_zero_x() {
        let jimenez = JimenezVelocityDeflection::new(0.01, 0.05);
        let x = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        let y = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        
        let result = jimenez.function(
            x, y, 25.0, 0.06, 0.33, 126.0, 
            &HashMap::new()
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_jimenez_deflection_proportional_to_yaw() {
        let jimenez = JimenezVelocityDeflection::new(0.01, 0.05);
        let x = Array::from_shape_vec((1, 1), vec![100.0]).unwrap();
        let y = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        
        let result_10 = jimenez.function(x.clone(), y.clone(), 10.0, 0.06, 0.33, 126.0, &HashMap::new()).unwrap();
        let result_20 = jimenez.function(x.clone(), y.clone(), 20.0, 0.06, 0.33, 126.0, &HashMap::new()).unwrap();
        
        let ratio = result_20[[0, 0]] / result_10[[0, 0]];
        assert_relative_eq!(ratio, 2.0, epsilon = 0.01);
    }
}
