//! Peak shaving turbine model
//!
//! Reduces thrust and power at or near rated wind speeds
//! Based on turbulence intensity and peak_shaving_fraction parameters

use crate::types::Float;
use crate::types::Array2;
use crate::core::turbine::operation_models::base::*;
use crate::core::turbine::operation_models::helpers::*;
use super::SimpleTurbine;

#[derive(Debug, Clone, Default)]
pub struct PeakShavingTurbine;

impl OperationModel for PeakShavingTurbine {
    fn model_name(&self) -> &'static str {
        "peak-shaving"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let turbulence_intensities = ctx.turbulence_intensities.ok_or_else(|| anyhow::anyhow!("Turbulence intensities required for PeakShavingTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        // Get base values
        let simple = SimpleTurbine;
        let base_powers = simple.power(params, ctx)?;
        let base_thrust_coefficients = simple.thrust_coefficient(params, ctx)?;
        let base_ais = crate::core::turbine::operation_models::helpers::axial_induction_from_ct(&base_thrust_coefficients);

        // Calculate peak shaving thrust limit
        let mut max_allowable_thrust = ndarray::Array::zeros((n_findex, n_turbines));
        let wind_speeds = &params.power_table.wind_speeds;
        let ct_values = &params.thrust_table.values;
        let mut peak_normal_thrust_prime = 0.0;
        for (ws, ct) in wind_speeds.iter().zip(ct_values.iter()) {
            peak_normal_thrust_prime = peak_normal_thrust_prime.max(ws * ws * ct);
        }

        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            ctx.velocities,
            ctx.average_method,
            ctx.cubature_weights,
        )?;

        let mut result = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j, 0]];

                // Replace zeros with small value to avoid division
                let safe_vel = vel.max(0.01);

                // Calculate max allowable thrust
                max_allowable_thrust[[i, j]] = (1.0 - params.peak_shaving_fraction)
                    * peak_normal_thrust_prime
                    / safe_vel.powf(2.0);

                // Check if turbulence threshold is met
                let ti = if turbulence_intensities.shape()[0] > i && turbulence_intensities.shape()[1] > j {
                    let ti_slice = turbulence_intensities.slice(s![i, j, .., ..]);
                    if ti_slice.len() > 0 {
                        ti_slice.iter().sum::<f64>() / ti_slice.len() as f64
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                let should_limit = ti >= params.peak_shaving_ti_threshold;

                if should_limit {
                    max_allowable_thrust[[i, j]] = max_allowable_thrust[[i, j]];
                } else {
                    max_allowable_thrust[[i, j]] = base_thrust_coefficients[[i, j]];
                }

                // Compute power fraction
                let base_ai = base_ais[[i, j]];
                let peak_shaving_ai = base_ais[[i, j]];

                let power_fraction = if should_limit {
                    (base_thrust_coefficients[[i, j]] * (1.0 - peak_shaving_ai))
                        / (base_thrust_coefficients[[i, j]] * (1.0 - base_ai))
                } else {
                    1.0
                };

                result[[i, j]] = base_powers[[i, j]] * power_fraction;
            }
        }

        Ok(result)
    }

    fn thrust_coefficient(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let turbulence_intensities = ctx.turbulence_intensities.ok_or_else(|| anyhow::anyhow!("Turbulence intensities required for PeakShavingTurbine"))?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.velocities.shape()[1];

        let rotor_avg_velocities = crate::core::rotor_velocity::average_velocity(
            ctx.velocities,
            ctx.average_method,
            ctx.cubature_weights,
        )?;

        let wind_speeds = &params.power_table.wind_speeds;
        let ct_values = &params.thrust_table.values;
        let mut peak_normal_thrust_prime = 0.0;
        for (ws, ct) in wind_speeds.iter().zip(ct_values.iter()) {
            peak_normal_thrust_prime = peak_normal_thrust_prime.max(ws * ws * ct);
        }

        let mut thrust_coeff = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                let vel = rotor_avg_velocities[[i, j, 0]];

                // Replace zeros with small value to avoid division
                let safe_vel = vel.max(0.01);

                // Calculate max allowable thrust
                let mut max_allowable = (1.0 - params.peak_shaving_fraction)
                    * peak_normal_thrust_prime
                    / safe_vel.powf(2.0);

                // Get base thrust coefficient
                let base_ct = params.thrust_table.interpolate(vel);

                // Check if turbulence threshold is met
                let ti = if turbulence_intensities.shape()[0] > i && turbulence_intensities.shape()[1] > j {
                    let ti_slice = turbulence_intensities.slice(s![i, j, .., ..]);
                    if ti_slice.len() > 0 {
                        ti_slice.iter().sum::<f64>() / ti_slice.len() as f64
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                let should_limit = ti >= params.peak_shaving_ti_threshold;

                if should_limit {
                    thrust_coeff[[i, j]] = max_allowable.min(base_ct);
                } else {
                    thrust_coeff[[i, j]] = base_ct.clamp(0.0001, 0.9999);
                }
            }
        }

        Ok(thrust_coeff)
    }

    fn axial_induction(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let ct = self.thrust_coefficient(params, ctx)?;
        Ok(crate::core::turbine::operation_models::helpers::axial_induction_from_ct(&ct))
    }
}
