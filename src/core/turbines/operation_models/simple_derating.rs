//! Simple derating turbine operation model
//!
//! Power setpoint control that limits turbine power output
//! When derated, thrust coefficient scales with power^(2/3)

use super::SimpleTurbine;
use crate::{
    core::{TurbineContext, TurbineParameters, turbines::{POWER_SETPOINT_DEFAULT, operation_models::OperationModel}},
    types::Array2,
};

/// Simple derating: limit power to specified setpoint
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
                        // Scale Ct by power^(2/3): P ∝ v³, T ∝ v²
                        // So if power reduces to X% of baseline, Ct reduces to (X%)^(2/3)
                        let power_ratio = power_setpoints[[i, j]] / base_power[[i, j]];
                        result[[i, j]] = base_ct[[i, j]] * power_ratio.powf(2.0 / 3.0);
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
