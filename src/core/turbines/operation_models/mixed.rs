//! Mixed operation turbine model
//!
//! Combines yaw derating and power setpoints for mixed control
//! If yaw angle is non-zero, use yaw control (CosineLossTurbine)
//! If power setpoint is set, use derating control (SimpleDeratingTurbine)
//! Otherwise use simple control (SimpleTurbine)

use crate::core::turbines::POWER_SETPOINT_DEFAULT;
use crate::core::turbines::operation_models::simple::SimpleTurbine;
use crate::core::turbines::operation_models::simple_derating::SimpleDeratingTurbine;
use crate::core::turbines::operation_models::OperationModel;
use crate::core::{TurbineContext, TurbineParameters};
use crate::types::Array2;

#[derive(Debug, Clone, Default)]
pub struct MixedOperationTurbine;

impl OperationModel for MixedOperationTurbine {
    fn model_name(&self) -> &'static str {
        "mixed"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for MixedOperationTurbine"))?;
        let power_setpoints = ctx
            .power_setpoints
            .ok_or_else(|| anyhow::anyhow!("Power setpoints required for MixedOperationTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        // Create masks
        let mut yaw_mask = ndarray::Array::from_elem((n_findex, n_turbines), false);
        let mut derating_mask = ndarray::Array::from_elem((n_findex, n_turbines), false);
        let mut simple_mask = ndarray::Array::from_elem((n_findex, n_turbines), true);

        for i in 0..n_findex {
            for j in 0..n_turbines {
                // Check if yaw angle is non-zero (use yaw control)
                if yaw_angles[[i, j]] != 0.0 {
                    yaw_mask[[i, j]] = true;
                    derating_mask[[i, j]] = false;
                    simple_mask[[i, j]] = false;
                } else {
                    simple_mask[[i, j]] = true;
                    derating_mask[[i, j]] = false;
                    yaw_mask[[i, j]] = false;
                }

                // Check if power setpoint is set (use derating control)
                if power_setpoints[[i, j]] < POWER_SETPOINT_DEFAULT {
                    yaw_mask[[i, j]] = false;
                    derating_mask[[i, j]] = true;
                    simple_mask[[i, j]] = false;
                }
            }
        }

        // Initialize result array
        let mut result = ndarray::Array::zeros((n_findex, n_turbines));

        // Get yaw-controlled power
        if yaw_mask.iter().any(|x| *x) {
            let mut yaw_ctx = ctx.clone();
            yaw_ctx.yaw_angles = Some(yaw_angles);
            yaw_ctx.power_setpoints = None;
            let cosine = crate::core::turbines::operation_models::cosine_loss::CosineLossTurbine;
            let yaw_power = cosine.power(&params, &yaw_ctx)?;

            for i in 0..n_findex {
                for j in 0..n_turbines {
                    if yaw_mask[[i, j]] {
                        result[[i, j]] = yaw_power[[i, j]];
                    }
                }
            }
        }

        // Get derating-controlled power
        if derating_mask.iter().any(|x| *x) {
            let mut derating_ctx = ctx.clone();
            derating_ctx.power_setpoints = Some(power_setpoints);
            let zero_yaw = ndarray::Array::zeros((n_findex, n_turbines));
            derating_ctx.yaw_angles = Some(&zero_yaw);
            let simple_derating =
                crate::core::turbines::operation_models::simple_derating::SimpleDeratingTurbine;
            let derating_power = simple_derating.power(&params, &derating_ctx)?;

            for i in 0..n_findex {
                for j in 0..n_turbines {
                    if derating_mask[[i, j]] {
                        result[[i, j]] = derating_power[[i, j]];
                    }
                }
            }
        }

        // Get simple power (no yaw, no derating)
        if simple_mask.iter().any(|x| *x) {
            let simple = crate::core::turbines::operation_models::simple::SimpleTurbine;
            let simple_power = simple.power(&params, ctx)?;

            for i in 0..n_findex {
                for j in 0..n_turbines {
                    if simple_mask[[i, j]] {
                        result[[i, j]] = simple_power[[i, j]];
                    }
                }
            }
        }

        Ok(result)
    }

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for MixedOperationTurbine"))?;
        let power_setpoints = ctx
            .power_setpoints
            .ok_or_else(|| anyhow::anyhow!("Power setpoints required for MixedOperationTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        // Create masks
        let mut yaw_mask = ndarray::Array::from_elem((n_findex, n_turbines), false);
        let mut derating_mask = ndarray::Array::from_elem((n_findex, n_turbines), false);

        for i in 0..n_findex {
            for j in 0..n_turbines {
                // Check if yaw angle is non-zero
                if yaw_angles[[i, j]] != 0.0 {
                    yaw_mask[[i, j]] = true;
                    derating_mask[[i, j]] = false;
                }
            }

            // Check if power setpoint is set
            for j in 0..n_turbines {
                if power_setpoints[[i, j]] < POWER_SETPOINT_DEFAULT {
                    derating_mask[[i, j]] = true;
                }
            }
        }

        // Get base thrust coefficient from SimpleTurbine
        let simple = crate::core::turbines::operation_models::simple::SimpleTurbine;
        let base_ct = simple.thrust_coefficient(params, ctx)?;

        // Apply derating effect if needed
        let mut result = base_ct.clone();
        if derating_mask.iter().any(|x| *x) {
            let base_powers = simple.power(params, ctx)?;

            for i in 0..n_findex {
                for j in 0..n_turbines {
                    if derating_mask[[i, j]] {
                        // Scale Ct by power^(2/3) as in simple derating model
                        let power_ratio = power_setpoints[[i, j]] / base_powers[[i, j]];
                        result[[i, j]] = base_ct[[i, j]] * power_ratio.powf(2.0 / 3.0);
                    }
                }
            }
        }

        Ok(result)
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
