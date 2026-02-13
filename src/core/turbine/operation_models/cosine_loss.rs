//! Cosine loss turbine operation model
//!
//! Actuator disk model with yaw/tilt misalignment handling
//! Power loss to yawing is defined by cos(yaw)^p where p is cosine_loss_exponent_yaw

use crate::types::Float;
use crate::types::{Array2, Array4};
use crate::core::turbine::operation_models::base::*;
use crate::core::turbine::operation_models::helpers::*;

#[derive(Debug, Clone, Default)]
pub struct CosineLossTurbine;

impl OperationModel for CosineLossTurbine {
    fn model_name(&self) -> &'static str {
        "cosine-loss"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let yaw_angles = ctx.yaw_angles.ok_or_else(|| anyhow::anyhow!("Yaw angles required for CosineLossTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            ctx.velocities,
            ctx.average_method,
            ctx.cubature_weights,
        )?;

        let mut power = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j, 0]];
                let air_density_correction = (ctx.air_density[i] / params.ref_air_density).powf(1.0 / 3.0);
                let mut effective_vel = vel * air_density_correction;

                // Apply yaw cosine correction
                let yaw = yaw_angles[[i, j]];
                let pw = params.cosine_loss_exponent_yaw / 3.0;
                effective_vel *= crate::utilities::cosd(yaw).powf(pw);

                power[[i, j]] = params.power_table.interpolate(effective_vel) * 1000.0;

                // Apply tilt correction if available
                if let Some(tilt_angles) = ctx.tilt_angles {
                    let ref_tilt_rad = params.ref_tilt.to_radians();
                    let tilt = tilt_angles[[i, j]].to_radians();
                    effective_vel *= crate::utilities::cosd(tilt).powf(pw) / crate::utilities::cosd(ref_tilt_rad);
                }

                // Apply power setpoints if available
                if let Some(power_setpoints) = ctx.power_setpoints {
                    if power_setpoints[[i, j]] < POWER_SETPOINT_DEFAULT {
                        power[[i, j]] = power[[i, j]].min(power_setpoints[[i, j]]);
                    }
                }
            }
        }

        Ok(power)
    }

    fn thrust_coefficient(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let yaw_angles = ctx.yaw_angles.ok_or_else(|| anyhow::anyhow!("Yaw angles required for CosineLossTurbine"))?;
        let tilt_angles = ctx.tilt_angles.ok_or_else(|| anyhow::anyhow!("Tilt angles required for CosineLossTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            ctx.velocities,
            ctx.average_method,
            ctx.cubature_weights,
        )?;

        let mut thrust_coeff = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j, 0]];
                let ct = params.thrust_table.interpolate(vel);
                thrust_coeff[[i, j]] = ct.clamp(0.0001, 0.9999);

                // Apply yaw correction
                let yaw = yaw_angles[[i, j]];
                thrust_coeff[[i, j]] *= crate::utilities::cosd(yaw);

                // Apply tilt correction if available
                let tilt = tilt_angles[[i, j]].to_radians();
                let ref_tilt_rad = params.ref_tilt.to_radians();
                thrust_coeff[[i, j]] *= crate::utilities::cosd(tilt) / crate::utilities::cosd(ref_tilt_rad);
            }
        }

        Ok(thrust_coeff)
    }

    fn axial_induction(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let ct = self.thrust_coefficient(params, ctx)?;
        let yaw_angles = ctx.yaw_angles.ok_or_else(|| anyhow::anyhow!("Yaw angles required for CosineLossTurbine"))?;
        let tilt_angles = ctx.tilt_angles.ok_or_else(|| anyhow::anyhow!("Tilt angles required for CosineLossTurbine"))?;

        let n_findex = ct.nrows();
        let n_turbines = ct.ncols();

        // Compute misalignment loss factor
        let mut misalignment_loss = ndarray::Array::ones((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let yaw = yaw_angles[[i, j]].to_radians();
                let tilt = tilt_angles[[i, j]].to_radians();
                let ref_tilt_rad = params.ref_tilt.to_radians();
                misalignment_loss[[i, j]] = crate::utilities::cosd(yaw.to_degrees()) * crate::utilities::cosd(tilt.to_degrees()) / crate::utilities::cosd(ref_tilt_rad.to_degrees());
            }
        }

        // Unified axial induction formula for yawed actuator disks
        let mut ai = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let ct = ct[[i, j]];
                ai[[i, j]] = 0.5 / misalignment_loss[[i, j]] * (1.0 - (ct * misalignment_loss[[i, j]]).sqrt());
            }
        }

        Ok(ai)
    }
}
