//! Unified Momentum Model turbine operation model
//!
//! Port of Python UnifiedMomentumModelTurbine (unified_momentum_model.py),
//! the yawed actuator-disk model of Heck et al. (2023). Assumes the modified
//! thrust coefficient C_T' is invariant to yaw:
//!
//! - `umm_rotor_axial_induction`: C_T' from Eq. (19) of Calaf et al. (2010),
//!   then the Heck fixed-point iteration for the yawed induction `an`.
//! - `umm_rotor_velocity_yaw_correction`: ratio of yaw-adjusted to yaw-aligned
//!   rotor wind speed (Eq. 2.3 / Eq. D1 of Heck et al.).
//!
//! The Heck solver is the fixed-point iteration defined at
//! unified_momentum_model.py:451 — max_iter=500, tolerance=1e-5,
//! relaxation=0.1 — with the LimitedHeck (Eq. 2.19/2.20) initial guess. No
//! defensive clamps: unconverged solves propagate NaN exactly as Python does.

use crate::core::turbines::operation_models::OperationModel;
use crate::core::{TurbineContext, TurbineParameters};
use crate::types::{Array2, Array4, Float};

/// 统一动量模型
/// 由麻省理工学院（MIT）开发。基于偏航致动盘动量理论（Heck et al., 2023），通过迭代求解来修正偏航状态下的轴向诱导因子和转子有效风速。它从动量守恒的角度出发，为偏航工况提供了一种不同于经验公式的物理建模路径。
#[derive(Debug, Clone, Default)]
pub struct UnifiedMomentumTurbine;

/// Limited case v_4 << u_4 (Eq. 2.19, 2.20 of Heck et al.). Returns
/// (a, u4, v4).
#[inline]
fn limited_heck(ct_prime: Float, yaw: Float) -> (Float, Float, Float) {
    let cos_yaw = yaw.cos();
    let sin_yaw = yaw.sin();
    let c2 = ct_prime * cos_yaw * cos_yaw;
    let a = c2 / (4.0 + c2);
    let u4 = (4.0 - c2) / (4.0 + c2);
    let v4 = -(4.0 * ct_prime * sin_yaw * cos_yaw * cos_yaw) / (4.0 + c2).powi(2);
    (a, u4, v4)
}

/// Heck fixed-point iteration on the residual of Eq. 2.15 (Python
/// `@fixedpointiteration(max_iter=500, tolerance=0.00001, relaxation=0.1)`,
/// unified_momentum_model.py:451). Returns the yawed axial induction `an`,
/// or NaN if the iteration did not converge (Python post_process).
fn heck_an(ct_prime: Float, yaw: Float) -> Float {
    let eps = 0.00001;
    let maxiter = 500;
    let relax = 0.1;

    let (mut a, mut u4, mut v4) = limited_heck(ct_prime, yaw);
    let mut converged = false;

    for _ in 0..maxiter {
        // residual(x, Ct', yaw) — Eq. 2.15
        let e_a = 1.0 - (1.0 - u4 * u4 - v4 * v4).sqrt() / (ct_prime.sqrt() * yaw.cos()) - a;
        let e_u4 = (1.0 - 0.5 * ct_prime * (1.0 - a) * yaw.cos() * yaw.cos()) - u4;
        let e_v4 =
            -0.25 * ct_prime * (1.0 - a) * (1.0 - a) * yaw.sin() * yaw.cos() * yaw.cos() - v4;

        // x0 = x0 + (1 - relax) * residual
        a += (1.0 - relax) * e_a;
        u4 += (1.0 - relax) * e_u4;
        v4 += (1.0 - relax) * e_v4;

        if e_a.abs() < eps && e_u4.abs() < eps && e_v4.abs() < eps {
            converged = true;
            break;
        }
    }

    if converged {
        a
    } else {
        Float::NAN
    }
}

/// UMM_rotor_axial_induction (unified_momentum_model.py:28-47).
#[inline]
fn umm_rotor_axial_induction(ct: Float, yaw_deg: Float) -> Float {
    let ai_yawaligned = 0.5 * (1.0 - (1.0 - ct).sqrt()); // Actuator disc theory
    let ct_prime = ct / (1.0 - ai_yawaligned).powi(2); // Eq. (19) of Calaf et al., 2010
    heck_an(ct_prime, yaw_deg.to_radians())
}

/// UMM_rotor_velocity_yaw_correction (unified_momentum_model.py:49-79).
#[allow(clippy::too_many_arguments)]
fn umm_rotor_velocity_yaw_correction(
    ct: Float,
    yaw_deg: Float,
    axial_induction: Float,
    rotor_effective_velocity: Float,
) -> Float {
    let ai_yawaligned = 0.5 * (1.0 - (1.0 - ct).sqrt());
    let ct_prime = ct / (1.0 - ai_yawaligned).powi(2);

    let yaw = yaw_deg.to_radians();
    let u_d_yawed = (1.0 - axial_induction) * yaw.cos(); // Eq. 2.3 of Heck et al.
    let u_d_aligned = 1.0 - ct_prime / (ct_prime + 4.0); // Eq. D1 of Heck et al.

    (u_d_yawed / u_d_aligned) * rotor_effective_velocity
}

impl UnifiedMomentumTurbine {
    /// Shared computation: (ct at v_eff, UMM axial induction, corrected
    /// rotor-effective velocity). Ct lookup has fill_value=0.0001 outside the
    /// table and NO clamp (Python interp1d directly).
    fn solve(
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<(Array2, Array2, Array2, Array2)> {
        let yaw_angles = ctx
            .yaw_angles
            .ok_or_else(|| anyhow::anyhow!("Yaw angles required for UnifiedMomentumTurbine"))?;

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

        // rotor_velocity_air_density_correction: v * (rho/rho_ref)^(1/3)
        let mut rotor_effective_velocities = ndarray::Array::zeros((n_findex, n_turbines));
        let mut thrust_coefficients = ndarray::Array::zeros((n_findex, n_turbines));
        let mut axial_inductions = ndarray::Array::zeros((n_findex, n_turbines));
        let mut corrected_velocities = ndarray::Array::zeros((n_findex, n_turbines));

        for i in 0..n_findex {
            let density_factor = (ctx.air_density[i] / params.ref_air_density).powf(1.0 / 3.0);
            for j in 0..n_turbines {
                let v_eff = rotor_avg_velocities[[i, j]] * density_factor;
                rotor_effective_velocities[[i, j]] = v_eff;

                let ct = params.thrust_table.interpolate(v_eff);
                thrust_coefficients[[i, j]] = ct;

                let yaw = yaw_angles[[i, j]];
                let ai = umm_rotor_axial_induction(ct, yaw);
                axial_inductions[[i, j]] = ai;

                corrected_velocities[[i, j]] =
                    umm_rotor_velocity_yaw_correction(ct, yaw, ai, v_eff);
            }
        }

        Ok((
            rotor_effective_velocities,
            thrust_coefficients,
            axial_inductions,
            corrected_velocities,
        ))
    }
}

impl OperationModel for UnifiedMomentumTurbine {
    fn model_name(&self) -> &'static str {
        "unified-momentum"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        // Python: power_interpolator(corrected_rotor_effective_velocities)*1e3
        // (no yaw/tilt cosine corrections — TODO in the Python source).
        let (.., corrected) = UnifiedMomentumTurbine::solve(params, ctx)?;
        let mut power = ndarray::Array::zeros(corrected.dim());
        for ((p, v), table) in power
            .iter_mut()
            .zip(corrected.iter())
            .zip(std::iter::repeat(&params.power_table))
        {
            *p = table.interpolate(*v) * 1000.0;
        }
        Ok(power)
    }

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        // Python: yawed_thrust_coefficients = interpolator(corrected velocity)
        let (.., corrected) = UnifiedMomentumTurbine::solve(params, ctx)?;
        let mut yawed_ct = ndarray::Array::zeros(corrected.dim());
        for ((c, v), table) in yawed_ct
            .iter_mut()
            .zip(corrected.iter())
            .zip(std::iter::repeat(&params.thrust_table))
        {
            *c = table.interpolate(*v);
        }
        Ok(yawed_ct)
    }

    fn axial_induction(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let (_, _, axial, _) = UnifiedMomentumTurbine::solve(params, ctx)?;
        Ok(axial)
    }

    fn clone_box(&self) -> Box<dyn OperationModel> {
        Box::new(self.clone())
    }
}
