//! Simple turbine operation model
//!
//! Basic actuator disk turbine model with no yaw or tilt handling

use crate::types::Float;
use crate::types::Array2;
use crate::core::turbine::operation_models::base::*;
use crate::core::turbine::operation_models::helpers::*;

#[derive(Debug, Clone, Default)]
pub struct SimpleTurbine;

impl OperationModel for SimpleTurbine {
    fn model_name(&self) -> &'static str {
        "simple"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
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
                let effective_vel = vel * air_density_correction;
                power[[i, j]] = params.power_table.interpolate(effective_vel) * 1000.0;
            }
        }

        Ok(power)
    }

    fn thrust_coefficient(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
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
                
                // Apply air density correction for thrust (uses 1/2 power because thrust ~ v²)
                let air_density_correction = (ctx.air_density[i] / params.ref_air_density).powf(1.0 / 2.0);
                let effective_vel = vel * air_density_correction;
                
                let ct = params.thrust_table.interpolate(effective_vel);
                thrust_coeff[[i, j]] = ct.clamp(0.0001, 0.9999);
            }
        }

        Ok(thrust_coeff)
    }

    fn axial_induction(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let ct = self.thrust_coefficient(params, ctx)?;
        Ok(axial_induction_from_ct(&ct))
    }
}
