/// Empirical Gauss wake deflection model
///
/// Extended Gauss model with empirical corrections for different conditions
use crate::types::{Float, Array1, Array2};
use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{GridBase, FlowField};
use std::collections::HashMap;
use ndarray::Array;

/// Empirical Gauss wake deflection model parameters
#[derive(Debug, Clone)]
pub struct EmpiricalGaussVelocityDeflection {
    pub base: BaseModel,
    pub kd: Float,
    pub ad: Float,
    pub bd: Float,
}

impl EmpiricalGaussVelocityDeflection {
    pub fn new(kd: Float, ad: Float, bd: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("kd".to_string(), kd);
        params.insert("ad".to_string(), ad);
        params.insert("bd".to_string(), bd);
        
        Self {
            base: BaseModel::new(params, "wind_vector"),
            kd,
            ad,
            bd,
        }
    }
}

impl DeflectionModel for EmpiricalGaussVelocityDeflection {
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
        turbulence_intensity: Float,
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
                
                let c1 = self.ad / (self.kd * 2.0);
                let exp_term = (-self.kd * x_downstream / rotor_diameter).exp();
                let base_deflection = c1 * axial_induction * (1.0 - exp_term) * yaw_rad;
                
                let ti_correction = 1.0 + self.bd * (turbulence_intensity - 0.1).max(0.0);
                
                deflection[[fi, ti]] = base_deflection * ti_correction;
            }
        }

        Ok(deflection)
    }
}

impl EmpiricalGaussVelocityDeflection {
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
    fn test_empirical_gauss_deflection_creation() {
        let empirical = EmpiricalGaussVelocityDeflection::new(0.01, 0.05, 0.5);
        assert_eq!(empirical.kd, 0.01);
        assert_eq!(empirical.ad, 0.05);
        assert_eq!(empirical.bd, 0.5);
    }
}
