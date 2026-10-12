//! Active Wake Control (AWC) turbine operation model
//!
//! Port of Python AWCTurbine (operation_models.py:510-623): for turbines in
//! `helix` mode, base power/Ct are modified as
//! `base * (1 - (b + c*base) * amplitude^a)` with the helix tuning parameters
//! from the turbine YAML; `baseline` turbines pass through unchanged.
//! Python raises UserWarning (an exception) for unknown modes and for
//! at/near-rated operation — mirrored here as errors.

use crate::core::turbines::operation_models::simple::SimpleTurbine;
use crate::core::turbines::operation_models::OperationModel;
use crate::core::{TurbineContext, TurbineParameters};
use crate::types::{Array2, Float};

#[derive(Debug, Clone, Default)]
pub struct AWCTurbine;

/// Python AWCTurbine.AWC_model (operation_models.py:522-523)
#[inline]
fn awc_model(a: Float, b: Float, c: Float, base_value: Float, awc_amplitude: Float) -> Float {
    base_value * (1.0 - (b + c * base_value) * awc_amplitude.powf(a))
}

impl AWCTurbine {
    fn check_modes(awc_modes: &ndarray::Array2<String>) -> crate::Result<()> {
        for mode in awc_modes.iter() {
            if mode != "helix" && mode != "baseline" {
                anyhow::bail!(
                    "Active wake mixing strategies other than the `helix` strategy have not \
                     yet been implemented in FLORIS. Returning baseline power."
                );
            }
        }
        Ok(())
    }
}

impl OperationModel for AWCTurbine {
    fn model_name(&self) -> &'static str {
        "awc"
    }

    fn power(&self, params: &TurbineParameters, ctx: &TurbineContext) -> crate::Result<Array2> {
        let awc_modes = ctx
            .awc_modes
            .ok_or_else(|| anyhow::anyhow!("AWC modes required for AWCTurbine"))?;
        let awc_amplitudes = ctx
            .awc_amplitudes
            .ok_or_else(|| anyhow::anyhow!("AWC amplitudes required for AWCTurbine"))?;
        AWCTurbine::check_modes(awc_modes)?;

        let base_powers = SimpleTurbine.power(params, ctx)?;

        let (n_findex, n_turbines) = base_powers.dim();
        let mut powers = base_powers.clone();

        for i in 0..n_findex {
            for j in 0..n_turbines {
                if awc_modes[[i, j]] != "helix" {
                    continue;
                }
                // Python raises if the helix turbine operates at/near rated:
                // np.isclose(base_powers[helix]/1000, max(power)) with default
                // rtol=1e-5, atol=1e-8.
                let base_kw = base_powers[[i, j]] / 1000.0;
                let max_kw = params
                    .power_table
                    .values
                    .iter()
                    .cloned()
                    .fold(Float::NEG_INFINITY, Float::max);
                if (base_kw - max_kw).abs() <= 1e-8 + 1e-5 * max_kw.abs() {
                    anyhow::bail!(
                        "The selected wind speed is above or near rated wind speed. AWCTurbine \
                         operation model is not designed or verified for above-rated conditions."
                    );
                }
                powers[[i, j]] = awc_model(
                    params.helix_a,
                    params.helix_power_b,
                    params.helix_power_c,
                    base_powers[[i, j]],
                    awc_amplitudes[[i, j]],
                );
            }
        }

        Ok(powers)
    }

    fn thrust_coefficient(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        let awc_modes = ctx
            .awc_modes
            .ok_or_else(|| anyhow::anyhow!("AWC modes required for AWCTurbine"))?;
        let awc_amplitudes = ctx
            .awc_amplitudes
            .ok_or_else(|| anyhow::anyhow!("AWC amplitudes required for AWCTurbine"))?;
        AWCTurbine::check_modes(awc_modes)?;

        let base_thrust_coefficients = SimpleTurbine.thrust_coefficient(params, ctx)?;

        let (n_findex, n_turbines) = base_thrust_coefficients.dim();
        let mut thrust_coefficients = base_thrust_coefficients.clone();

        for i in 0..n_findex {
            for j in 0..n_turbines {
                if awc_modes[[i, j]] != "helix" {
                    continue;
                }
                thrust_coefficients[[i, j]] = awc_model(
                    params.helix_a,
                    params.helix_thrust_b,
                    params.helix_thrust_c,
                    base_thrust_coefficients[[i, j]],
                    awc_amplitudes[[i, j]],
                );
            }
        }

        Ok(thrust_coefficients)
    }

    fn axial_induction(
        &self,
        params: &TurbineParameters,
        ctx: &TurbineContext,
    ) -> crate::Result<Array2> {
        // Python: (1 - sqrt(1 - ct)) / 2 on the AWC-adjusted Ct
        let ct = self.thrust_coefficient(params, ctx)?;
        Ok(crate::core::turbines::operation_models::helpers::axial_induction_from_ct(&ct))
    }

    fn clone_box(&self) -> Box<dyn OperationModel> {
        Box::new(self.clone())
    }
}
