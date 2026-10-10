use crate::core::wake::{BaseModel, DeflectionModel};
use crate::core::{FlowField, Grid};
/// Jimenez wake deflection model
///
/// Based on Jimenez et al. (2010) - wake deflection due to yaw misalignment
use crate::types::{Array1, Array2, Array4, Float};
use ndarray::Array;
use std::collections::HashMap;

/// Jimenez wake deflection model parameters
#[derive(Debug, Clone)]
pub struct JimenezVelocityDeflection {
    pub base: BaseModel,
    pub kd: Float,
    pub ad: Float,
    pub bd: Float,
}

impl JimenezVelocityDeflection {
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

impl DeflectionModel for JimenezVelocityDeflection {
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
        x_i: &Array2,
        _y_i: &Array2,
        yaw_i: &Array1,
        _turbulence_intensity_i: &Array4,
        thrust_coefficient_i: &Array1,
        rotor_diameter: Float,
        grid_x: &Array4,
        _model_args: &HashMap<String, Array1>,
    ) -> anyhow::Result<Array4> {
        let shape = grid_x.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];
        let grid_y_dim = shape[2];
        let grid_z_dim = shape[3];

        let mut deflection = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));

        for fi in 0..n_findex {
            // Port of Python jimenez.py: xi_init = cosd(yaw) * sind(yaw) * ct / 2
            let yaw = yaw_i[fi];
            let xi_init = crate::utilities::cosd(yaw)
                * crate::utilities::sind(yaw)
                * thrust_coefficient_i[fi]
                / 2.0;
            let x_src = x_i[[fi, 0]];

            for ti in 0..n_turbines {
                for iy in 0..grid_y_dim {
                    for iz in 0..grid_z_dim {
                        let delta_x = grid_x[[fi, ti, iy, iz]] - x_src;

                        // Numexpr port of the Python reference implementation
                        let a = 15.0 * (2.0 * self.kd * delta_x / rotor_diameter + 1.0).powi(4)
                            + xi_init.powi(2);
                        let b = (30.0 * self.kd / rotor_diameter)
                            * (2.0 * self.kd * delta_x / rotor_diameter + 1.0).powi(5);
                        let c = xi_init * rotor_diameter * (15.0 + xi_init.powi(2));
                        let d = 30.0 * self.kd;

                        let y_yaw_init = (xi_init * a / b) - (c / d);
                        deflection[[fi, ti, iy, iz]] = y_yaw_init + self.ad + self.bd * delta_x;
                    }
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

    fn make_model() -> JimenezVelocityDeflection {
        JimenezVelocityDeflection::new(0.05, 0.0, 0.0)
    }

    fn deflection_at(
        model: &JimenezVelocityDeflection,
        yaw: f64,
        x_downstream: f64,
        ct: f64,
    ) -> f64 {
        let x_i = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        let y_i = Array::from_shape_vec((1, 1), vec![0.0]).unwrap();
        let yaw_arr = Array1::from_vec(vec![yaw]);
        let ti = Array::from_shape_vec((1, 1, 1, 1), vec![0.1]).unwrap();
        let ct = Array1::from_vec(vec![ct]);
        let grid_x = Array::from_shape_vec((1, 1, 1, 1), vec![x_downstream]).unwrap();

        let result = model
            .function(
                &x_i,
                &y_i,
                &yaw_arr,
                &ti,
                &ct,
                126.0,
                &grid_x,
                &HashMap::new(),
            )
            .unwrap();
        result[[0, 0, 0, 0]]
    }

    #[test]
    fn test_jimenez_deflection_creation() {
        let jimenez = make_model();
        assert_eq!(jimenez.kd, 0.05);
        assert_eq!(jimenez.ad, 0.0);
        assert_eq!(jimenez.bd, 0.0);
    }

    #[test]
    fn test_jimenez_deflection_at_zero_x() {
        // At the rotor plane the yaw-induced term is canceled by the -C/D term
        let jimenez = make_model();
        assert_relative_eq!(
            deflection_at(&jimenez, 10.0, 0.0, 0.33),
            0.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_jimenez_deflection_grows_downstream() {
        // Golden values from Python FLORIS v4.6.4 JimenezVelocityDeflection:
        // yaw=+5, ct=0.78, D=126, kd=0.05, x_i=0. Positive yaw deflects the
        // wake toward negative y, and the deflection magnitude grows
        // downstream.
        let jimenez = make_model();
        let d_1d = deflection_at(&jimenez, 5.0, 126.0, 0.78);
        let d_5d = deflection_at(&jimenez, 5.0, 5.0 * 126.0, 0.78);
        assert_relative_eq!(d_1d, -3.87990514, epsilon = 1e-6);
        assert_relative_eq!(d_5d, -14.2246176, epsilon = 1e-6);
        assert!(d_5d < d_1d); // more negative further downstream
    }

    #[test]
    fn test_jimenez_deflection_odd_in_yaw() {
        // xi_init is odd in yaw, so deflection reverses sign with the yaw
        let jimenez = make_model();
        let pos = deflection_at(&jimenez, 10.0, 5.0 * 126.0, 0.33);
        let neg = deflection_at(&jimenez, -10.0, 5.0 * 126.0, 0.33);
        assert_relative_eq!(pos, -neg, epsilon = 1e-12);
    }
}
