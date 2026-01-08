/// Crespo-Hernandez wake turbulence model
///
/// Based on Crespo et al. with Hernandez - wake-added turbulence intensity

use crate::types::{Float, Array1};
use crate::core::wake::{BaseModel, TurbulenceModel};
use crate::core::{GridBase, FlowField};
use std::collections::HashMap;

/// Crespo-Hernandez wake turbulence model parameters
#[derive(Debug, Clone)]
pub struct CrespoHernandez {
    pub base: BaseModel,
    pub c: Float,
    pub p: Float,
}

impl CrespoHernandez {
    pub fn new(c: Float, p: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("c".to_string(), c);
        params.insert("p".to_string(), p);
        
        Self {
            base: BaseModel::new(params, "scalar"),
            c,
            p,
        }
    }
}

impl TurbulenceModel for CrespoHernandez {
    fn prepare_function(
        &self,
        _grid: &dyn GridBase,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        Ok(HashMap::new())
    }

    fn function(
        &self,
        ambient_turbulence_intensity: Float,
        rotor_diameter_eff: Float,
        downstream_distance: Array1,
        turbine_type_parameters: &HashMap<String, Float>,
        _model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array1> {
        let c = turbine_type_parameters.get("c")
            .copied()
            .unwrap_or(self.c);
        let _p = turbine_type_parameters.get("p")
            .copied()
            .unwrap_or(self.p);

        let shape = downstream_distance.shape();
        let n_distances = shape[0];

        let mut wake_added_ti = Array1::zeros(n_distances);

        for i in 0..n_distances {
            let x = downstream_distance[i];
            
            if x <= 0.0 {
                wake_added_ti[i] = 0.0;
                continue;
            }

            let x_normalized = x / rotor_diameter_eff;
            
            if x_normalized > 0.0 {
                let added_ti = c * (x_normalized).sqrt();
                wake_added_ti[i] = added_ti.min(0.5);
            }
        }

        let total_ti: Array1 = wake_added_ti.mapv(|ti_added| {
            (ambient_turbulence_intensity.powi(2) + ti_added.powi(2)).sqrt()
        });

        Ok(total_ti)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_crespo_hernandez_creation() {
        let crespo = CrespoHernandez::new(0.9, 0.9);
        assert_eq!(crespo.c, 0.9);
        assert_eq!(crespo.p, 0.9);
    }

    #[test]
    fn test_crespo_hernandez_turbulence() {
        let crespo = CrespoHernandez::new(0.9, 0.9);
        let ambient_ti = 0.06;
        let rotor_d = 126.0;
        let downstream = Array1::from_vec(vec![100.0, 200.0, 500.0]);
        
        let result = crespo.function(
            ambient_ti,
            rotor_d,
            downstream,
            &HashMap::new(),
            &HashMap::new()
        );
        assert!(result.is_ok());
        
        let ti = result.unwrap();
        for &val in ti.iter() {
            assert!(val >= ambient_ti);
        }
    }

    #[test]
    fn test_crespo_hernandez_at_zero_distance() {
        let crespo = CrespoHernandez::new(0.9, 0.9);
        let ambient_ti = 0.06;
        let rotor_d = 126.0;
        let downstream = Array1::from_vec(vec![0.0]);
        
        let result = crespo.function(
            ambient_ti,
            rotor_d,
            downstream,
            &HashMap::new(),
            &HashMap::new()
        );
        assert!(result.is_ok());
        
        let ti = result.unwrap();
        assert_relative_eq!(ti[0], ambient_ti);
    }
}
