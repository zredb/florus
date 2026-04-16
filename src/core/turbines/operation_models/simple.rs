//! Simple turbine operation model
//!
//! Basic actuator disk turbine model with no yaw or tilt handling

use crate::core::turbines::operation_models::helpers::*;
use crate::core::turbines::operation_models::OperationModel;
use crate::core::TurbineContext;
use crate::core::TurbineParameters;
use crate::types::Array2;
use crate::types::Float;

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
                let air_density_correction =
                    (ctx.air_density[i] / params.ref_air_density).powf(1.0 / 3.0);
                let effective_vel = vel * air_density_correction;
                power[[i, j]] = params.power_table.interpolate(effective_vel) * 1000.0;
            }
        }

        Ok(power)
    }

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
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
                let air_density_correction =
                    (ctx.air_density[i] / params.ref_air_density).powf(1.0 / 2.0);
                let effective_vel = vel * air_density_correction;

                let ct = params.thrust_table.interpolate(effective_vel);
                thrust_coeff[[i, j]] = ct.clamp(0.0001, 0.9999);
            }
        }

        Ok(thrust_coeff)
    }
    ///衡量风机让风减速程度的系数
    /// 最佳值：约 0.33（对应最大发电效率）
    /// 重要性：它是连接风机性能（发了多少电）和尾流效应（对后面风机影响多大）的桥梁。
    fn axial_induction(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let ct = self.thrust_coefficient(params, ctx)?;
        Ok(axial_induction_from_ct(&ct))
    }

    fn clone_box(&self) -> Box<dyn OperationModel> {
        Box::new(self.clone())
    }
}
