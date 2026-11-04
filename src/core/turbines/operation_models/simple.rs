//! Simple turbine operation model
//!
//! Basic actuator disk turbine model with no yaw or tilt handling

use crate::core::turbines::operation_models::helpers::*;
use crate::core::turbines::operation_models::OperationModel;
use crate::core::TurbineContext;
use crate::core::TurbineParameters;
use crate::types::Array2;
///简单模型
/// 最基础的基准模型。它完全依据风机自带的功率曲线和推力系数曲线来计算输出，不响应偏航、降额等任何控制指令。通常用于获取未受控制干预时的“自然”发电量，作为对比基准。
#[derive(Debug, Clone, Default)]
pub struct SimpleTurbine;

impl OperationModel for SimpleTurbine {
    fn model_name(&self) -> &'static str {
        "simple"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.n_turbines();

        let velocities_dyn = ctx.velocities.view().into_dyn().to_owned();
        let cubature_weights = ctx.cubature_weights_1d();
        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            &velocities_dyn,
            None,
            &ctx.average_method,
            cubature_weights.as_ref(),
        )?;

        let mut power = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j]];
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
        let n_turbines = ctx.n_turbines();

        let velocities_dyn = ctx.velocities.view().into_dyn().to_owned();
        let cubature_weights = ctx.cubature_weights_1d();
        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            &velocities_dyn,
            None,
            &ctx.average_method,
            cubature_weights.as_ref(),
        )?;

        let mut thrust_coeff = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                // Python SimpleTurbine.thrust_coefficient has NO air-density
                // correction on the Ct lookup (operation_models.py:131 TODO);
                // only power is density-corrected (exponent 1/3).
                let vel = rotor_avg_velocities[[i, j]];

                let ct = params.thrust_table.interpolate(vel);
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
