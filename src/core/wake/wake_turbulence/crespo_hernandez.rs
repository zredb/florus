//! Crespo-Hernandez wake turbulence model
//!
//! Wake-added turbulence intensity following Crespo et al., ported from
//! FLORIS v4.6.4 `core/wake_turbulence/crespo_hernandez.py`.

use crate::types::{Array1, Array4, Float};
use crate::core::wake::{BaseModel, TurbulenceModel};
use crate::core::{FlowField, Grid};
use std::collections::HashMap;

/// Crespo-Hernandez wake turbulence model parameters
#[derive(Debug, Clone)]
pub struct CrespoHernandez {
    pub base: BaseModel,
    pub initial: Float,
    pub constant: Float,
    pub ai: Float,
    pub downstream: Float,
}

impl CrespoHernandez {
    pub fn new(initial: Float, constant: Float, ai: Float, downstream: Float) -> Self {
        let mut params = HashMap::new();
        params.insert("initial".to_string(), initial);
        params.insert("constant".to_string(), constant);
        params.insert("ai".to_string(), ai);
        params.insert("downstream".to_string(), downstream);

        Self {
            base: BaseModel::new(params, "scalar"),
            initial,
            constant,
            ai,
            downstream,
        }
    }
}

impl TurbulenceModel for CrespoHernandez {
    fn prepare_function(
        &self,
        _grid: &dyn Grid,
        _flow_field: &FlowField,
    ) -> anyhow::Result<HashMap<String, Array1>> {
        Ok(HashMap::new())
    }

    fn clone_box(&self) -> Box<dyn TurbulenceModel> {
        Box::new(self.clone())
    }

    fn function(
        &self,
        ambient_turbulence_intensity: &Array4,
        x: &Array4,
        x_i: &Array4,
        rotor_diameter: Float,
        axial_induction: &Array4,
    ) -> anyhow::Result<Array4> {
        // delta_x = x - x_i  (x_i broadcasts over the [findex, 1, 1, 1] axis)
        let delta_x = x - x_i;

        // Upstream points (delta_x <= 0.1) are pinned to 1.0 so the power law
        // stays finite; downstream points keep their true distance.
        let upstream_mask = delta_x.mapv(|v| v <= 0.1);
        let downstream_mask = delta_x.mapv(|v| v > -0.1);
        let delta_x_eff = &delta_x * downstream_mask.mapv(|b| if b { 1.0 } else { 0.0 })
            + upstream_mask.mapv(|b| if b { 1.0 } else { 0.0 });

        let axial = axial_induction.mapv(|v| v.powf(self.ai));
        let ambient = ambient_turbulence_intensity.mapv(|v| v.powf(self.initial));
        let ratio = (&delta_x_eff / rotor_diameter).mapv(|v| v.powf(self.downstream));

        let ti = (axial * ambient * ratio).mapv(|v| self.constant * v);

        // Zero out the upstream region (delta_x <= -0.1); keep downstream added TI.
        let result = ti * downstream_mask.mapv(|b| if b { 1.0 } else { 0.0 });
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_crespo_hernandez_creation() {
        let crespo = CrespoHernandez::new(0.1, 0.9, 0.8, -0.32);
        assert_eq!(crespo.initial, 0.1);
        assert_eq!(crespo.constant, 0.9);
        assert_eq!(crespo.ai, 0.8);
        assert_eq!(crespo.downstream, -0.32);
    }

    #[test]
    fn test_crespo_hernandez_downstream_positive() {
        // Single findex, single turbine grid point downstream of source.
        let crespo = CrespoHernandez::new(0.1, 0.9, 0.8, -0.32);
        let ambient = Array4::from_elem((1, 1, 1, 1), 0.06);
        let x = Array4::from_elem((1, 1, 1, 1), 500.0);
        let x_i = Array4::from_elem((1, 1, 1, 1), 0.0);
        let ai = Array4::from_elem((1, 1, 1, 1), 0.25);
        let result = crespo.function(&ambient, &x, &x_i, 126.0, &ai).unwrap();
        // delta_x = 500 > 0 -> downstream, ti = 0.9 * 0.25^0.8 * 0.06^0.1 * (500/126)^(-0.32)
        let expected = 0.9 * 0.25_f64.powf(0.8) * 0.06_f64.powf(0.1) * (500.0 / 126.0_f64).powf(-0.32);
        assert_relative_eq!(result[[0, 0, 0, 0]], expected, epsilon = 1e-12);
    }

    #[test]
    fn test_crespo_hernandez_upstream_zero() {
        // Point upstream of the source gets masked to 0.
        let crespo = CrespoHernandez::new(0.1, 0.9, 0.8, -0.32);
        let ambient = Array4::from_elem((1, 1, 1, 1), 0.06);
        let x = Array4::from_elem((1, 1, 1, 1), -100.0);
        let x_i = Array4::from_elem((1, 1, 1, 1), 0.0);
        let ai = Array4::from_elem((1, 1, 1, 1), 0.25);
        let result = crespo.function(&ambient, &x, &x_i, 126.0, &ai).unwrap();
        assert_relative_eq!(result[[0, 0, 0, 0]], 0.0, epsilon = 1e-12);
    }
}
