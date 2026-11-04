//! Cosine loss turbine operation model
//!
//! Actuator disk model with yaw/tilt misalignment handling
//! Power loss to yawing is defined by cos(yaw)^p where p is cosine_loss_exponent_yaw

use crate::{
    core::{
        turbines::{operation_models::OperationModel, POWER_SETPOINT_DEFAULT},
        TurbineContext, TurbineParameters,
    },
    types::Array2,
};
///余弦损失模型
/// FLORIS 的默认模型。它用余弦函数来近似偏航或倾斜带来的气动损失：推力按 cos γ 衰减，功率按 (cos γ)^p 衰减（p 为可配置的损失指数）。这是进行常规尾流转向优化时最常用的模型。
#[derive(Debug, Clone, Default)]
pub struct CosineLossTurbine;

impl OperationModel for CosineLossTurbine {
    fn model_name(&self) -> &'static str {
        "cosine-loss"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for CosineLossTurbine"))?;

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
                let mut effective_vel = vel * air_density_correction;

                // Apply yaw cosine correction
                let yaw = yaw_angles[[i, j]];
                let pw_yaw = params.cosine_loss_exponent_yaw / 3.0;
                effective_vel *= crate::utilities::cosd(yaw).powf(pw_yaw);

                // Apply tilt cosine correction (Python
                // rotor_velocity_tilt_cosine_correction, applied BEFORE the
                // power table lookup). For floating turbines the tilt is
                // interpolated from the current rotor-effective velocity
                // (post density/yaw) when correct_cp_ct_for_tilt is set
                // (operation_models.py:208-215).
                if let Some(tilt_angles) = ctx.tilt_angles {
                    let static_tilt = tilt_angles[[i, j]];
                    let tilt = ctx.effective_tilt(params, i, j, effective_vel, static_tilt);
                    let pw_tilt = params.cosine_loss_exponent_tilt / 3.0;
                    effective_vel *= (crate::utilities::cosd(tilt)
                        / crate::utilities::cosd(params.ref_tilt))
                    .powf(pw_tilt);
                }

                power[[i, j]] = params.power_table.interpolate(effective_vel) * 1000.0;

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

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for CosineLossTurbine"))?;
        let tilt_angles = ctx
            .tilt_angles
            .ok_or_else(|| anyhow::anyhow!("Tilt angles required for CosineLossTurbine"))?;

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
                // Python CosineLossTurbine.thrust_coefficient builds on
                // SimpleTurbine.thrust_coefficient, which has NO air-density
                // correction on the Ct lookup (Bug #35).
                let vel = rotor_avg_velocities[[i, j]];

                let ct = params.thrust_table.interpolate(vel);
                thrust_coeff[[i, j]] = ct.clamp(0.0001, 0.9999);

                // Apply yaw/tilt corrections (Python CosineLossTurbine.
                // thrust_coefficient): ct * cosd(yaw) * cosd(tilt) / cosd(ref_tilt).
                // The tilt is interpolated from the UNCORRECTED rotor
                // average velocity for floating turbines requesting the
                // correction (operation_models.py:252-261).
                let yaw = yaw_angles[[i, j]];
                thrust_coeff[[i, j]] *= crate::utilities::cosd(yaw);

                let tilt = ctx.effective_tilt(params, i, j, vel, tilt_angles[[i, j]]);
                thrust_coeff[[i, j]] *=
                    crate::utilities::cosd(tilt) / crate::utilities::cosd(params.ref_tilt);
            }
        }

        Ok(thrust_coeff)
    }

    fn axial_induction(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        // Python CosineLossTurbine.axial_induction (operation_models.py:
        // 272-298) delegates the thrust coefficient to
        // `CosineLossTurbine.thrust_coefficient` — so for floating turbines
        // the Ct carries the velocity-interpolated tilt — while the
        // `misalignment_loss` below uses the ORIGINAL (static) tilt angles
        // from the caller (a Python quirk: the local tilt reassignment
        // inside thrust_coefficient does not propagate back):
        //   misalignment_loss = cosd(yaw) * cosd(tilt) / cosd(ref_tilt)
        //   ai = 0.5 / misalignment_loss * (1 - sqrt(1 - ct * misalignment_loss))

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.n_turbines();

        // Full Ct through the thrust path (yaw factor + effective tilt).
        let ct_full = self.thrust_coefficient(params, ctx)?;

        // Python defaults: yaw = 0, tilt = ref_tilt (farm.set_tilt_to_ref_tilt)
        let zero_yaw = Array2::zeros((n_findex, n_turbines));
        let ref_tilt_arr = Array2::from_elem((n_findex, n_turbines), params.ref_tilt);
        let yaw_angles = ctx.yaw_angles.unwrap_or(&zero_yaw);
        let tilt_angles = ctx.tilt_angles.unwrap_or(&ref_tilt_arr);

        let mut axial_induction = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let misalignment_loss = crate::utilities::cosd(yaw_angles[[i, j]])
                    * crate::utilities::cosd(tilt_angles[[i, j]])
                    / crate::utilities::cosd(params.ref_tilt);
                let ct = ct_full[[i, j]];
                axial_induction[[i, j]] =
                    0.5 / misalignment_loss * (1.0 - (1.0 - ct * misalignment_loss).sqrt());
            }
        }

        Ok(axial_induction)
    }

    fn clone_box(&self) -> Box<dyn OperationModel> {
        Box::new(self.clone())
    }
}
