//! Peak shaving turbine model
//!
//! Reduces thrust and power at or near rated wind speeds
//! Based on turbulence intensity and peak_shaving_fraction parameters

use crate::core::turbines::operation_models::simple::SimpleTurbine;
use crate::core::turbines::operation_models::OperationModel;
use crate::core::{TurbineContext, TurbineParameters};
use crate::types::{Array2, Float};

#[derive(Debug, Clone, Default)]
pub struct PeakShavingTurbine;

impl OperationModel for PeakShavingTurbine {
    fn model_name(&self) -> &'static str {
        "peak-shaving"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        // Fail fast if TI is missing (Python's power path receives it too; the
        // rotor-averaged value is consumed inside thrust_coefficient).
        ctx.turbulence_intensities.ok_or_else(|| {
            anyhow::anyhow!("Turbulence intensities required for PeakShavingTurbine")
        })?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.n_turbines();

        // Python PeakShavingTurbine.power (operation_models.py:628-673):
        //   base_powers = SimpleTurbine.power(...)
        //   base_ct = SimpleTurbine.thrust_coefficient(...)
        //   ps_ct = PeakShavingTurbine.thrust_coefficient(...)
        //   base_ais, ps_ais = (1 - sqrt(1 - ct)) / 2
        //   power_fraction = (ps_ct * (1 - ps_ai)) / (base_ct * (1 - base_ai))
        //   power = power_fraction * base_powers
        // Previously `peak_shaving_ai = base_ais` made the fraction always 1
        // and the peak-shaving Ct never entered the power path (Bug #35).
        let simple = SimpleTurbine;
        let base_powers = simple.power(params, ctx)?;
        let base_thrust_coefficients = simple.thrust_coefficient(params, ctx)?;
        let peak_shaving_thrust_coefficients = self.thrust_coefficient(params, ctx)?;

        let base_ais = crate::core::turbines::operation_models::helpers::axial_induction_from_ct(
            &base_thrust_coefficients,
        );
        let peak_shaving_ais =
            crate::core::turbines::operation_models::helpers::axial_induction_from_ct(
                &peak_shaving_thrust_coefficients,
            );

        let mut result = ndarray::Array::zeros((n_findex, n_turbines));
        for i in 0..n_findex {
            for j in 0..n_turbines {
                // When no shaving applies ps_ct == base_ct, so the fraction is
                // exactly 1 (same numerator and denominator), matching Python.
                let bc = base_thrust_coefficients[[i, j]];
                let pc = peak_shaving_thrust_coefficients[[i, j]];
                let power_fraction =
                    (pc * (1.0 - peak_shaving_ais[[i, j]])) / (bc * (1.0 - base_ais[[i, j]]));
                result[[i, j]] = base_powers[[i, j]] * power_fraction;
            }
        }

        Ok(result)
    }

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let turbulence_intensities = ctx.turbulence_intensities.ok_or_else(|| {
            anyhow::anyhow!("Turbulence intensities required for PeakShavingTurbine")
        })?;

        let n_findex = ctx.air_density.len();
        let n_turbines = ctx.n_turbines();

        // Base Ct via the simple model (Python: SimpleTurbine.thrust_coefficient,
        // no air-density correction on the Ct lookup).
        let simple = SimpleTurbine;
        let base_thrust_coefficients = simple.thrust_coefficient(params, ctx)?;

        // peak_normal_thrust_prime = max(ws^2 * ct) over the power-thrust table
        // (Python operation_models.py:690-693)
        let wind_speeds = &params.power_table.keys;
        let ct_values = &params.thrust_table.values;
        let mut peak_normal_thrust_prime: Float = 0.0;
        for (ws, ct) in wind_speeds.iter().zip(ct_values.iter()) {
            peak_normal_thrust_prime = peak_normal_thrust_prime.max(ws * ws * ct);
        }

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
                // np.maximum(rotor_average_velocities, 0.01) — no density correction
                let safe_vel = rotor_avg_velocities[[i, j]].max(0.01);

                // max_allowable_thrust_coefficient =
                //     (1 - peak_shaving_fraction) * peak_normal_thrust_prime / v^2
                let max_allowable = (1.0 - params.peak_shaving_fraction) * peak_normal_thrust_prime
                    / safe_vel.powf(2.0);

                // Rotor-averaged TI of this turbine (Python:
                // turbulence_intensities.mean(axis=(2, 3)))
                let mut ti_sum = 0.0;
                let mut ti_count = 0.0;
                if turbulence_intensities.shape()[0] > i && turbulence_intensities.shape()[1] > j {
                    for iy in 0..turbulence_intensities.shape()[2] {
                        for iz in 0..turbulence_intensities.shape()[3] {
                            ti_sum += turbulence_intensities[[i, j, iy, iz]];
                            ti_count += 1.0;
                        }
                    }
                }
                let ti = if ti_count > 0.0 {
                    ti_sum / ti_count
                } else {
                    0.0
                };

                // np.where(mean TI >= peak_shaving_TI_threshold, max_allowable, base_ct)
                let capped = if ti >= params.peak_shaving_ti_threshold {
                    max_allowable
                } else {
                    base_thrust_coefficients[[i, j]]
                };

                // np.minimum(base_ct, capped)
                thrust_coeff[[i, j]] = base_thrust_coefficients[[i, j]].min(capped);
            }
        }

        Ok(thrust_coeff)
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
