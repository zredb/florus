/// Turbine model and operations
///
/// Corresponds to core/turbine/ module in Python implementation
use crate::{
    core::{
        rotor_effective_velocity,
        turbines::{cp_ct_table::TableConditions, turbine_type::TurbineType},
        AveragingMethod,
    },
    types::{Array2, Array4, Float},
};
use ndarray::Array;

/// Wind turbine representation with embedded turbine type
#[derive(Debug, Clone)]
pub struct Turbine {
    /// Turbine type containing all turbine parameters and curves
    pub turbine_type: TurbineType,
}

impl Turbine {
    /// Calculate power output for given velocities
    ///
    /// # Arguments
    /// * `velocities` - Wind velocities at turbine rotor (findex, n_turbines, n_points, n_grid)
    /// * `yaw_angles` - Yaw angles [degrees]
    /// * `tilt_angles` - Tilt angles [degrees]
    /// * `average_method` - Method for averaging velocities
    pub fn calculate_power(
        &self,
        velocities: &Array4,
        _air_density: Float,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        average_method: AveragingMethod,
    ) -> crate::Result<Array2> {
        // Calculate rotor effective velocity
        // Create ref_tilt array for single turbine
        let n_findex = velocities.shape()[0];
        let n_turbines = velocities.shape()[1];
        let ref_tilt = Array::from_elem(
            (n_findex, n_turbines),
            self.turbine_type.power_thrust_table.ref_tilt.unwrap_or(5.0),
        );
        let correct_cp_ct_for_tilt = Array::from_elem(
            (n_findex, n_turbines),
            self.turbine_type.correct_cp_ct_for_tilt,
        );

        let rotor_velocities = rotor_effective_velocity(
            velocities,
            average_method,
            None,
            yaw_angles,
            tilt_angles,
            Some(&ref_tilt),
            None,
            None,
            Some(&correct_cp_ct_for_tilt),
        )?;

        let shape = rotor_velocities.shape();
        let mut power_output = Array::zeros((shape[0], shape[1]));

        // Calculate power for each turbine at each condition
        for fi in 0..shape[0] {
            for ti in 0..shape[1] {
                let v = rotor_velocities[[fi, ti, 0]];

                let wind_speed = self
                    .turbine_type
                    .power_thrust_table
                    .cp_ct_table
                    .wind_speeds();
                if fi == 0 && ti == 0 {
                    log::warn!("DEBUG: v={}, ws[0]={}, ws[-1]={}", v, wind_speed[0], wind_speed[wind_speed.len() - 1]);
                }

                if v < wind_speed[0] || v > wind_speed[wind_speed.len() - 1] {
                    power_output[[fi, ti]] = 0.0;
                } else {
                    let mut conditions = TableConditions::default();
                    conditions.wind_speed = v;
                    let power_kw = self.turbine_type.power_thrust_table.cp_ct_table.get_cp(&conditions).unwrap();
                    let power_watts = power_kw * 1000.0;
                    if fi == 0 && ti == 0 {
                        eprintln!("DEBUG: v={}, power_kw={}, power_watts={}", v, power_kw, power_watts);
                    }

                    // Apply yaw correction from operation model
                    if let Some(yaw) = yaw_angles {
                        let yaw_rad = yaw[[fi, ti]].to_radians();
                        let operation_model = self.turbine_type.operation_model.as_ref();
                        let loss_factor = operation_model.power_loss_factor(yaw_rad, 2.0);
                        power_output[[fi, ti]] = power_watts * loss_factor;
                    } else {
                        power_output[[fi, ti]] = power_watts;
                    }
                }
            }
        }

        Ok(power_output)
    }

    /// Calculate thrust coefficient for given velocities
    pub fn calculate_thrust_coefficient(
        &self,
        velocities: &Array4,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        average_method: AveragingMethod,
    ) -> crate::Result<Array2> {
        // Create ref_tilt array for single turbine
        let n_findex = velocities.shape()[0];
        let n_turbines = velocities.shape()[1];
        let ref_tilt = Array::from_elem(
            (n_findex, n_turbines),
            self.turbine_type.power_thrust_table.ref_tilt.unwrap_or(5.0),
        );
        let correct_cp_ct_for_tilt = Array::from_elem(
            (n_findex, n_turbines),
            self.turbine_type.correct_cp_ct_for_tilt,
        );

        let rotor_velocities = rotor_effective_velocity(
            velocities,
            average_method,
            None,
            yaw_angles,
            tilt_angles,
            Some(&ref_tilt),
            None,
            None,
            Some(&correct_cp_ct_for_tilt),
        )?;

        let shape = rotor_velocities.shape();
        let mut ct_output = Array::zeros((shape[0], shape[1]));

        for fi in 0..shape[0] {
            for ti in 0..shape[1] {
                let v = rotor_velocities[[fi, ti, 0]];
                let mut conditions = TableConditions::default();
                conditions.wind_speed = v;
                ct_output[[fi, ti]] = self
                    .turbine_type
                    .power_thrust_table
                    .cp_ct_table
                    .get_ct(&conditions)
                    .unwrap_or(0.0);
            }
        }

        Ok(ct_output)
    }

    /// Calculate axial induction factor from thrust coefficient
    pub fn calculate_axial_induction(&self, ct: Float) -> Float {
        // Using momentum theory relationship
        if ct < 0.96 {
            0.5 * (1.0 - (1.0 - ct).sqrt())
        } else {
            // High thrust region - empirical relationship
            0.143 + (0.0203 - 0.6427 * (0.889 - ct)).sqrt()
        }
    }
}
