//! Simple derating turbine operation model
//!
//! Power setpoint control that limits turbine power output
//! Thrust coefficient scales linearly with the power fraction,
//! capped at the base Ct (Python operation_models.py:358-360)

use super::SimpleTurbine;
use crate::{
    core::{
        turbines::{operation_models::OperationModel, POWER_SETPOINT_DEFAULT},
        TurbineContext, TurbineParameters,
    },
    types::Array2,
};

/// 简单降额模型
/// 允许你为每台风机指定一个功率设定点（power_setpoints），风机的实际输出会被限制在该值以内。它不涉及复杂的控制器逻辑，只是简单地“压低”功率和相应的推力。
#[derive(Debug, Clone, Default)]
pub struct SimpleDeratingTurbine;

impl OperationModel for SimpleDeratingTurbine {
    fn model_name(&self) -> &'static str {
        "simple-derating"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        // Get base power without derating
        let simple = SimpleTurbine;
        let base_power = simple.power(params, ctx)?;

        // Apply power setpoint limit if specified
        if let Some(power_setpoints) = ctx.power_setpoints {
            let mut result = base_power.clone();
            for i in 0..result.nrows() {
                for j in 0..result.ncols() {
                    if power_setpoints[[i, j]] < POWER_SETPOINT_DEFAULT {
                        result[[i, j]] = result[[i, j]].min(power_setpoints[[i, j]]);
                    }
                }
            }
            Ok(result)
        } else {
            Ok(base_power)
        }
    }

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let simple = SimpleTurbine;
        let base_ct = simple.thrust_coefficient(params, ctx)?;
        let base_power = simple.power(params, ctx)?;

        // Apply power setpoint if specified
        if let Some(power_setpoints) = ctx.power_setpoints {
            let mut result = base_ct.clone();
            for i in 0..result.nrows() {
                for j in 0..result.ncols() {
                    if power_setpoints[[i, j]] < POWER_SETPOINT_DEFAULT {
                        // Python scales thrust linearly with the power fraction
                        // and caps it at the base Ct
                        // (operation_models.py:358-360); the previous
                        // power_ratio^(2/3) heuristic overestimated derated Ct
                        // (Bug #29).
                        let power_fraction = power_setpoints[[i, j]] / base_power[[i, j]];
                        let scaled = power_fraction * base_ct[[i, j]];
                        result[[i, j]] = base_ct[[i, j]].min(scaled);
                    }
                }
            }
            Ok(result)
        } else {
            Ok(base_ct)
        }
    }

    fn axial_induction(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let ct = self.thrust_coefficient(params, ctx)?;
        Ok(crate::core::turbines::operation_models::helpers::axial_induction_from_ct(&ct))
    }

    fn clone_box(&self) -> Box<dyn OperationModel> {
        Box::new(self.clone())
    }
}
