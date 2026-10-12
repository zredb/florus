//! Mixed operation turbine model
//!
//! Port of Python MixedOperationTurbine (operation_models.py:383-507):
//! dispatches per (findex, turbine) between yaw control (CosineLossTurbine),
//! power-setpoint control (SimpleDeratingTurbine) and plain SimpleTurbine.
//! Disabled turbines (setpoint <= POWER_SETPOINT_DISABLED) get their yaw
//! zeroed; a turbine with BOTH a nonzero yaw and an active setpoint raises.

use crate::core::turbines::operation_models::cosine_loss::CosineLossTurbine;
use crate::core::turbines::operation_models::simple::SimpleTurbine;
use crate::core::turbines::operation_models::simple_derating::SimpleDeratingTurbine;
use crate::core::turbines::operation_models::OperationModel;
use crate::core::turbines::{POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED};
use crate::core::{TurbineContext, TurbineParameters};
use crate::types::Array2;

#[derive(Debug, Clone, Default)]
pub struct MixedOperationTurbine;

/// Masks returned by `handle_mixed_operation_setpoints` (Python
/// `_handle_mixed_operation_setpoints`, operation_models.py:482-507).
struct MixedMasks {
    yaw_mask: ndarray::Array2<bool>,
    setpoint_mask: ndarray::Array2<bool>,
    neither_mask: ndarray::Array2<bool>,
}

/// Zero the yaw of disabled turbines, build the dispatch masks, and raise on
/// incompatible yaw + setpoint combinations (Python raises ValueError).
fn handle_mixed_operation_setpoints(
    yaw_angles: &Array2,
    power_setpoints: &Array2,
) -> crate::Result<(Array2, MixedMasks)> {
    let mut yaw = yaw_angles.clone();
    let (n_findex, n_turbines) = yaw.dim();

    // Disabled turbines operate with zero yaw
    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            if power_setpoints[[fi, ti]] <= POWER_SETPOINT_DISABLED {
                yaw[[fi, ti]] = 0.0;
            }
        }
    }

    let mut yaw_mask = ndarray::Array2::<bool>::from_elem((n_findex, n_turbines), false);
    let mut setpoint_mask = ndarray::Array2::<bool>::from_elem((n_findex, n_turbines), false);
    let mut neither_mask = ndarray::Array2::<bool>::from_elem((n_findex, n_turbines), false);

    for fi in 0..n_findex {
        for ti in 0..n_turbines {
            let is_yaw = yaw[[fi, ti]] != 0.0;
            let is_setpoint = power_setpoints[[fi, ti]] < POWER_SETPOINT_DEFAULT;
            if is_yaw && is_setpoint {
                anyhow::bail!(
                    "Power setpoints and yaw angles are incompatible. If yaw_angles entry \
                     is nonzero, power_setpoints must be greater than or equal to {}.",
                    POWER_SETPOINT_DEFAULT
                );
            }
            yaw_mask[[fi, ti]] = is_yaw;
            setpoint_mask[[fi, ti]] = is_setpoint;
            neither_mask[[fi, ti]] = !is_yaw && !is_setpoint;
        }
    }

    Ok((
        yaw,
        MixedMasks {
            yaw_mask,
            setpoint_mask,
            neither_mask,
        },
    ))
}

impl MixedOperationTurbine {
    /// Dispatch a sub-model over the given mask: run it on the full context
    /// and scatter its values into the masked positions (Python builds the
    /// sub-model's full array and indexes it with the mask).
    fn dispatch(
        masks: &MixedMasks,
        call: impl Fn() -> crate::Result<Array2>,
        mask_of: fn(&MixedMasks) -> &ndarray::Array2<bool>,
    ) -> crate::Result<Array2> {
        let mask = mask_of(masks);
        let (n_findex, n_turbines) = mask.dim();
        let mut result = ndarray::Array::zeros((n_findex, n_turbines));
        if mask.iter().any(|x| *x) {
            let sub = call()?;
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    if mask[[fi, ti]] {
                        result[[fi, ti]] = sub[[fi, ti]];
                    }
                }
            }
        }
        Ok(result)
    }
}

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
        let (yaw, masks) = handle_mixed_operation_setpoints(yaw_angles, power_setpoints)?;

        // Yaw-controlled: CosineLossTurbine.power with the (zeroed) yaw
        let yaw_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = None;
                CosineLossTurbine.power(params, &sub_ctx)
            },
            |m| &m.yaw_mask,
        )?;

        // Setpoint-controlled: SimpleDeratingTurbine.power
        let derating_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = Some(power_setpoints);
                SimpleDeratingTurbine.power(params, &sub_ctx)
            },
            |m| &m.setpoint_mask,
        )?;

        // Neither: SimpleTurbine.power
        let simple_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = None;
                SimpleTurbine.power(params, &sub_ctx)
            },
            |m| &m.neither_mask,
        )?;

        Ok(yaw_result + derating_result + simple_result)
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
        let (yaw, masks) = handle_mixed_operation_setpoints(yaw_angles, power_setpoints)?;

        // Yaw-controlled: CosineLossTurbine.thrust_coefficient (applies the
        // cosd(yaw)*cosd(tilt)/cosd(ref_tilt) correction — previously the
        // mixed model gave yaw-controlled turbines the plain Simple Ct, Bug #36)
        let yaw_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = None;
                CosineLossTurbine.thrust_coefficient(params, &sub_ctx)
            },
            |m| &m.yaw_mask,
        )?;

        // Setpoint-controlled: SimpleDeratingTurbine.thrust_coefficient
        // (linear power-fraction scaling with a base-Ct cap, Bug #29)
        let derating_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = Some(power_setpoints);
                SimpleDeratingTurbine.thrust_coefficient(params, &sub_ctx)
            },
            |m| &m.setpoint_mask,
        )?;

        // Neither: SimpleTurbine.thrust_coefficient
        let simple_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = None;
                SimpleTurbine.thrust_coefficient(params, &sub_ctx)
            },
            |m| &m.neither_mask,
        )?;

        Ok(yaw_result + derating_result + simple_result)
    }

    fn axial_induction(
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
        let (yaw, masks) = handle_mixed_operation_setpoints(yaw_angles, power_setpoints)?;

        let yaw_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = None;
                CosineLossTurbine.axial_induction(params, &sub_ctx)
            },
            |m| &m.yaw_mask,
        )?;

        let derating_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = Some(power_setpoints);
                SimpleDeratingTurbine.axial_induction(params, &sub_ctx)
            },
            |m| &m.setpoint_mask,
        )?;

        let simple_result = Self::dispatch(
            &masks,
            || {
                let mut sub_ctx = ctx.clone();
                sub_ctx.yaw_angles = Some(&yaw);
                sub_ctx.power_setpoints = None;
                SimpleTurbine.axial_induction(params, &sub_ctx)
            },
            |m| &m.neither_mask,
        )?;

        Ok(yaw_result + derating_result + simple_result)
    }

    fn clone_box(&self) -> Box<dyn OperationModel> {
        Box::new(self.clone())
    }
}
