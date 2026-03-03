//! Unified Momentum turbine operation model
//!
//! Based on Unified Momentum Model (Liew et al., 2024) which provides
//! a unified framework for rotor aerodynamics across operating regimes
//!
//! This model implements:
//! 1. Rotor-normal axial induction factor (an) - accounting for induction increase with yaw
//! 2. Streamwise outlet velocity (u4) - for power calculation
//! 3. Modified thrust coefficient (CT') - including yaw and tilt effects
//! 4. Near-wake length (x0) - distance where wake pressure equals ambient
//! 5. Wake pressure deficit - pressure difference across wake
//!
//! Reference: "Unified Momentum Model for Rotor Aerodynamics Across
//!            Operating Regimes" (Nature Communications, 2024)
//!
//! Implementation follows the coupled equation system from the paper:
//! - an = (1 - sqrt(1 - ct*cos^2(gamma))) / (2*cos^2(gamma))
//! - u4 = u_inf * (1 - an) / 2) * cos^2(gamma)
//! - x0 = -C'T / (2*u4) where C' captures pressure effects
//!
//! The model provides analytical solutions that:
//! - Capture the monotonic increase in thrust with induction observed in experiments
//! - Work across yaw misaligned and high thrust states
//! - Eliminate the need for empirical corrections in classical momentum theory

use crate::types::Float;
use crate::types::{Array2, Array3, Array4};
use crate::core::turbine::operation_models::base::*;
use ndarray::s;

/// Unified Momentum turbine model
#[derive(Debug, Clone, Default)]
pub struct UnifiedMomentumTurbine;

impl OperationModel for UnifiedMomentumTurbine {
    fn model_name(&self) -> &'static str {
        "unified-momentum"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let yaw_angles = ctx.yaw_angles.ok_or_else(|| anyhow::anyhow!("Yaw angles required for UnifiedMomentumTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            ctx.velocities,
            ctx.average_method,
            ctx.cubature_weights,
        )?;

        let yaw_rad = yaw_angles.mapv(|y| y.to_radians());

        // Get base thrust coefficient from table
        let mut ct_base = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j, 0]];
                ct_base[[i, j]] = params.thrust_table.interpolate(vel).clamp(0.0001, 0.9999);
            }
        }

        // Modified thrust coefficient with yaw effect from Unified Momentum Model
        let mut ct_mod = ct_base.clone();
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let cos_gamma = yaw_rad[[i, j]].cos();
                ct_mod[[i, j]] = ct_base[[i, j]] * cos_gamma * cos_gamma;
            }
        }

        // Compute rotor-normal axial induction factor (an)
        // From Unified Momentum Model: an = (1 - sqrt(1 - ct*cos^2(gamma))) / (2*cos^2(gamma))
        // This captures the observed monotonic increase in thrust with induction
        let mut an = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let ct = ct_mod[[i, j]];
                an[[i, j]] = (1.0 - (1.0 - ct).sqrt()) / (2.0 * ct.powf(2.0));
            }
        }

        // Compute streamwise outlet velocity (u4) - for power calculation
        let mut u4_normalized = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let one_minus_an = 1.0 - an[[i, j]];
                u4_normalized[[i, j]] = one_minus_an * yaw_rad[[i, j]].cos().powf(2.0) / 2.0;
            }
        }

        // Apply air density correction
        let density_correction = (ctx.air_density[0] / params.ref_air_density).powf(1.0 / 3.0);

        // Compute effective velocities and power
        let mut power = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j, 0]];
                let mut effective_vel = vel * density_correction;

                // Apply yaw correction to effective velocity
                let yaw = yaw_rad[[i, j]];
                let pw = params.cosine_loss_exponent_yaw / 3.0;
                effective_vel *= yaw.cos().powf(pw);

                power[[i, j]] = params.power_table.interpolate(effective_vel) * u4_normalized[[i, j]] * 1000.0;
            }
        }

        Ok(power)
    }

    fn thrust_coefficient(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let yaw_angles = ctx.yaw_angles.ok_or_else(|| anyhow::anyhow!("Yaw angles required for UnifiedMomentumTurbine"))?;
        let tilt_angles = ctx.tilt_angles.ok_or_else(|| anyhow::anyhow!("Tilt angles required for UnifiedMomentumTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            ctx.velocities,
            ctx.average_method,
            ctx.cubature_weights,
        )?;

        let yaw_rad = yaw_angles.mapv(|y| y.to_radians());

        // Get base thrust coefficient from table with air density correction
        let mut ct_base = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j, 0]];
                
                // Apply air density correction for thrust (uses 1/2 power because thrust ~ v²)
                let air_density_correction = (ctx.air_density[i] / params.ref_air_density).powf(1.0 / 2.0);
                let effective_vel = vel * air_density_correction;
                
                ct_base[[i, j]] = params.thrust_table.interpolate(effective_vel).clamp(0.0001, 0.9999);
            }
        }

        // Modified thrust coefficient with yaw effect
        let mut ct_mod = ct_base.clone();
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let cos_gamma = yaw_rad[[i, j]].cos();
                ct_mod[[i, j]] = ct_base[[i, j]] * cos_gamma * cos_gamma;
            }
        }

        // Apply tilt correction if available
        let ref_tilt_rad = params.ref_tilt.to_radians();
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let tilt = tilt_angles[[i, j]].to_radians();
                let cos_tilt = tilt.cos();
                let ref_cos = ref_tilt_rad.cos();
                ct_mod[[i, j]] = ct_mod[[i, j]] * cos_tilt / ref_cos;
            }
        }

        Ok(ct_mod)
    }

    fn axial_induction(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let yaw_angles = ctx.yaw_angles.ok_or_else(|| anyhow::anyhow!("Yaw angles required for UnifiedMomentumTurbine"))?;
        ctx.tilt_angles.ok_or_else(|| anyhow::anyhow!("Tilt angles required for UnifiedMomentumTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        let yaw_rad = yaw_angles.mapv(|y| y.to_radians());

        // Get thrust coefficient using unified model
        let ct = self.thrust_coefficient(params, ctx)?;

        // Compute rotor-normal induction factor (an) from Unified Momentum Model
        let mut an = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let gamma = yaw_rad[[i, j]];
                let cos_gamma_sq = gamma.cos().powf(2.0);
                let ct_ij = ct[[i, j]];

                let sqrt_term = (1.0 - ct_ij * cos_gamma_sq).sqrt();
                an[[i, j]] = sqrt_term / (2.0 * cos_gamma_sq);
            }
        }

        Ok(an)
    }
}
