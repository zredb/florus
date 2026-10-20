/// Turbine model and operations
///
/// Corresponds to core/turbine/ module in Python implementation
use crate::{
    core::{
        rotor_effective_velocity,
        turbines::{
            cp_ct_table::{MultiDimConditionTables, TableConditions},
            operation_models::{TurbineContext, TurbineParameters},
            turbine_type::TurbineType,
        },
        AveragingMethod,
    },
    types::{Array1, Array2, Array4, Float},
};
use ndarray::Array;
use std::collections::HashMap;
use std::sync::Arc;

/// The multi-dimensional conditions of a model: condition name -> one value
/// per wind condition (Python `flow_field.multidim_conditions`).
pub type MultidimConditionMap = HashMap<String, Array1>;

/// Wind turbine representation with embedded turbine type
#[derive(Debug, Clone)]
pub struct Turbine {
    /// Turbine type containing all turbine parameters and curves
    /// Using Arc for shared ownership - TurbineType is immutable configuration data
    pub turbine_type: Arc<TurbineType>,
}

impl Turbine {
    /// Build the operation-model parameter set (lookup tables + scalars)
    /// from the turbine definition.
    ///
    /// Bails for multi-dimensional Cp/Ct tables: the operation-model dispatch
    /// only supports the 1-D `power`/`thrust_coefficient` curves (Python
    /// turbine.py routes multi-dim tables through a separate branch that
    /// FLORUS has not ported yet).
    /// The floating tilt interpolation table of this turbine, if any.
    ///
    /// Python `Turbine.tilt_interp` (turbine.py:642-648): built only when the
    /// turbine type requests the tilt correction. interp1d semantics: linear
    /// in range; below the first wind speed -> 0.0; above the last -> the
    /// last tilt angle.
    pub fn tilt_interp(&self) -> Option<crate::core::turbines::operation_models::LookupTable> {
        if !self.turbine_type.correct_cp_ct_for_tilt {
            return None;
        }
        self.turbine_type.floating_tilt_table.as_ref().map(|t| {
            let last_tilt = t.tilt_angle[t.tilt_angle.len() - 1];
            crate::core::turbines::operation_models::LookupTable {
                keys: t.wind_speed.clone(),
                values: t.tilt_angle.clone(),
                fill_value: 0.0,
                fill_value_low: Some(0.0),
                fill_value_high: Some(last_tilt),
            }
        })
    }

    fn operation_model_parameters(&self) -> crate::Result<TurbineParameters> {
        self.operation_model_parameters_with_tables(None)
    }

    /// Build the operation-model parameters, optionally overriding the
    /// power/thrust lookup tables.
    ///
    /// The override is used by the multi-dimensional power/thrust tables,
    /// where the table depends on the wind condition
    /// (Python `floris/core/turbine/turbine.py:power`).
    fn operation_model_parameters_with_tables(
        &self,
        tables: Option<(&Array1, &Array1, &Array1)>,
    ) -> crate::Result<TurbineParameters> {
        // Extract the wind speeds, CP and CT values of the table
        let (wind_speeds, cp_values, ct_values) = match tables {
            Some((wind_speeds, cp_values, ct_values)) => (wind_speeds, cp_values, ct_values),
            None => match &self.turbine_type.power_thrust_table.cp_ct_table {
                crate::core::turbines::cp_ct_table::CpCtTable::OneDimensional(table) => {
                    (&table.wind_speeds, &table.cp_values, &table.ct_values)
                }
                _ => {
                    return Err(anyhow::anyhow!(
                        "Multi-dimensional Cp/Ct tables require the per-condition tables"
                    ));
                }
            },
        };

        let power_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: cp_values.clone(),
            fill_value: 0.0,
            ..Default::default()
        };

        let thrust_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: ct_values.clone(),
            fill_value: 0.0001,
            ..Default::default()
        };

        // Get ref_tilt value
        let ref_tilt = self.turbine_type.power_thrust_table.ref_tilt.unwrap_or(5.0);

        // Floating-tilt interpolation table (Python `Turbine.tilt_interp`,
        // turbine.py:642-648): built only when the turbine type requests the
        // tilt correction.
        let tilt_interp = self.tilt_interp();

        // Create turbine parameters
        Ok(TurbineParameters {
            power_table,
            thrust_table,
            ref_air_density: self
                .turbine_type
                .power_thrust_table
                .ref_air_density
                .unwrap_or(1.225),
            cosine_loss_exponent_yaw: self
                .turbine_type
                .power_thrust_table
                .cosine_loss_exponent_yaw
                .unwrap_or(2.0),
            cosine_loss_exponent_tilt: self
                .turbine_type
                .power_thrust_table
                .cosine_loss_exponent_tilt
                .unwrap_or(2.0),
            ref_tilt,
            peak_shaving_fraction: self
                .turbine_type
                .power_thrust_table
                .peak_shaving_fraction
                .unwrap_or(0.0),
            peak_shaving_ti_threshold: self
                .turbine_type
                .power_thrust_table
                .peak_shaving_ti_threshold
                .unwrap_or(0.0),
            helix_a: self
                .turbine_type
                .power_thrust_table
                .helix
                .as_ref()
                .map(|h| h.a)
                .unwrap_or(0.0),
            helix_power_b: self
                .turbine_type
                .power_thrust_table
                .helix
                .as_ref()
                .map(|h| h.helix_power_b)
                .unwrap_or(0.0),
            helix_power_c: self
                .turbine_type
                .power_thrust_table
                .helix
                .as_ref()
                .map(|h| h.helix_power_c)
                .unwrap_or(0.0),
            helix_thrust_b: self
                .turbine_type
                .power_thrust_table
                .helix
                .as_ref()
                .map(|h| h.helix_thrust_b)
                .unwrap_or(0.0),
            helix_thrust_c: self
                .turbine_type
                .power_thrust_table
                .helix
                .as_ref()
                .map(|h| h.helix_thrust_c)
                .unwrap_or(0.0),
            controller_dependent: self
                .turbine_type
                .power_thrust_table
                .controller_dependent_turbine_parameters
                .as_ref()
                .and_then(|p| {
                    // All scalar parameters and the loaded Cp/Ct surface must
                    // be present; otherwise the model stays unavailable.
                    let surface = p.cp_ct_surface.as_ref()?;
                    Some(
                        crate::core::turbines::operation_models::ControllerDependentParams {
                            sigma: p.rotor_solidity?,
                            beta: p.beta?,
                            cd: p.cd?,
                            cl_alfa: p.cl_alfa?,
                            rotor_radius: p.rotor_diameter? / 2.0,
                            rated_rpm: p.rated_rpm?,
                            generator_efficiency: p.generator_efficiency?,
                            rated_power_kw: p.rated_power?,
                            ref_air_density: self
                                .turbine_type
                                .power_thrust_table
                                .ref_air_density
                                .unwrap_or(1.225),
                            cp_lut: surface.cp_lut.clone(),
                            ct_lut: surface.ct_lut.clone(),
                            pitch_lut: surface.pitch_lut.clone(),
                            tsr_lut: surface.tsr_lut.clone(),
                        },
                    )
                }),
            tilt_interp,
        })
    }

    /// Calculate power output for given velocities
    ///
    /// # Arguments
    /// * `velocities` - Wind velocities at turbine rotor (findex, n_turbines, n_points, n_grid)
    /// * `air_density` - Air density [kg/m^3]
    /// * `yaw_angles` - Yaw angles [degrees]
    /// * `tilt_angles` - Tilt angles [degrees]
    /// * `average_method` - Method for averaging velocities
    /// * `cubature_weights` - Grid cubature weights (TurbineCubatureGrid only;
    ///   Python turbine.py passes grid.cubature_weights through to the
    ///   operation models' average_velocity calls)
    /// * `correct_cp_ct_for_tilt` - Per-(findex, turbine) request for the
    ///   floating-tilt Cp/Ct correction (Python farm.correct_cp_ct_for_tilt)
    pub fn power(
        &self,
        velocities: &Array4,
        air_density: Float,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        turbulence_intensities: Option<&Array4>,
        power_setpoints: Option<&Array2>,
        awc_modes: Option<&ndarray::Array2<String>>,
        awc_amplitudes: Option<&Array2>,
        average_method: AveragingMethod,
        cubature_weights: Option<&Array2>,
        correct_cp_ct_for_tilt: Option<&ndarray::Array2<bool>>,
        multidim_conditions: Option<&MultidimConditionMap>,
    ) -> crate::Result<Array2> {
        let n_findex = velocities.shape()[0];

        // Multi-dimensional power/thrust tables: the table is selected per
        // wind condition (Python floris/core/turbine/turbine.py:power loops
        // over the specified conditions and masks the result).
        if let Some(tables) = self.multidim_tables() {
            let selection = self.multidim_selection(tables, multidim_conditions, n_findex)?;
            let air_density_array = Array::from_elem(n_findex, air_density);
            let ctx = self.turbine_context(
                velocities,
                &air_density_array,
                yaw_angles,
                tilt_angles,
                power_setpoints,
                turbulence_intensities,
                awc_modes,
                awc_amplitudes,
                cubature_weights,
                correct_cp_ct_for_tilt,
                average_method,
            );
            let mut power = Array2::zeros((n_findex, ctx.n_turbines()));
            for group in unique_groups(&selection) {
                let params = self.operation_model_parameters_with_tables(Some((
                    &tables.entries[group].wind_speeds,
                    &tables.entries[group].power,
                    &tables.entries[group].thrust_coefficient,
                )))?;
                let group_power = self.turbine_type.operation_model.power(&params, &ctx)?;
                for fi in 0..n_findex {
                    if selection[fi] != group {
                        continue;
                    }
                    for ti in 0..ctx.n_turbines() {
                        power[[fi, ti]] = group_power[[fi, ti]];
                    }
                }
            }
            return Ok(power);
        }

        let params = self.operation_model_parameters()?;

        // Create air density array
        let air_density_array = Array::from_elem(n_findex, air_density);

        // Create turbine context
        let ctx = TurbineContext {
            velocities,
            air_density: &air_density_array,
            yaw_angles,
            tilt_angles,
            power_setpoints,
            turbulence_intensities,
            awc_modes,
            awc_amplitudes,
            cubature_weights,
            correct_cp_ct_for_tilt,
            average_method,
        };

        // Call the operation model's power method
        self.turbine_type.operation_model.power(&params, &ctx)
    }

    /// Calculate thrust coefficient for given velocities
    pub fn thrust_coefficient(
        &self,
        velocities: &Array4,
        air_density: Float,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        turbulence_intensities: Option<&Array4>,
        power_setpoints: Option<&Array2>,
        awc_modes: Option<&ndarray::Array2<String>>,
        awc_amplitudes: Option<&Array2>,
        average_method: AveragingMethod,
        cubature_weights: Option<&Array2>,
        correct_cp_ct_for_tilt: Option<&ndarray::Array2<bool>>,
        multidim_conditions: Option<&MultidimConditionMap>,
    ) -> crate::Result<Array2> {
        let n_findex = velocities.shape()[0];

        // Multi-dimensional power/thrust tables (Python turbine.py:thrust_coefficient)
        if let Some(tables) = self.multidim_tables() {
            let selection = self.multidim_selection(tables, multidim_conditions, n_findex)?;
            let air_density_array = Array::from_elem(n_findex, air_density);
            let ctx = self.turbine_context(
                velocities,
                &air_density_array,
                yaw_angles,
                tilt_angles,
                power_setpoints,
                turbulence_intensities,
                awc_modes,
                awc_amplitudes,
                cubature_weights,
                correct_cp_ct_for_tilt,
                average_method,
            );
            let mut thrust = Array2::zeros((n_findex, ctx.n_turbines()));
            for group in unique_groups(&selection) {
                let params = self.operation_model_parameters_with_tables(Some((
                    &tables.entries[group].wind_speeds,
                    &tables.entries[group].power,
                    &tables.entries[group].thrust_coefficient,
                )))?;
                let group_thrust = self
                    .turbine_type
                    .operation_model
                    .thrust_coefficient(&params, &ctx)?;
                for fi in 0..n_findex {
                    if selection[fi] != group {
                        continue;
                    }
                    for ti in 0..ctx.n_turbines() {
                        thrust[[fi, ti]] = group_thrust[[fi, ti]];
                    }
                }
            }
            return Ok(thrust);
        }

        let params = self.operation_model_parameters()?;

        // Create air density array
        let air_density_array = Array::from_elem(n_findex, air_density);

        // Create turbine context
        let ctx = TurbineContext {
            velocities,
            air_density: &air_density_array,
            yaw_angles,
            tilt_angles,
            power_setpoints,
            turbulence_intensities,
            awc_modes,
            awc_amplitudes,
            cubature_weights,
            correct_cp_ct_for_tilt,
            average_method,
        };

        // Call the operation model's thrust_coefficient method
        self.turbine_type
            .operation_model
            .thrust_coefficient(&params, &ctx)
    }

    /// Calculate axial induction factor for given velocities
    pub fn axial_induction(
        &self,
        velocities: &Array4,
        air_density: Float,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        turbulence_intensities: Option<&Array4>,
        power_setpoints: Option<&Array2>,
        awc_modes: Option<&ndarray::Array2<String>>,
        awc_amplitudes: Option<&Array2>,
        average_method: AveragingMethod,
        cubature_weights: Option<&Array2>,
        correct_cp_ct_for_tilt: Option<&ndarray::Array2<bool>>,
        multidim_conditions: Option<&MultidimConditionMap>,
    ) -> crate::Result<Array2> {
        let n_findex = velocities.shape()[0];

        // Multi-dimensional power/thrust tables (Python turbine.py:axial_induction)
        if let Some(tables) = self.multidim_tables() {
            let selection = self.multidim_selection(tables, multidim_conditions, n_findex)?;
            let air_density_array = Array::from_elem(n_findex, air_density);
            let ctx = self.turbine_context(
                velocities,
                &air_density_array,
                yaw_angles,
                tilt_angles,
                power_setpoints,
                turbulence_intensities,
                awc_modes,
                awc_amplitudes,
                cubature_weights,
                correct_cp_ct_for_tilt,
                average_method,
            );
            let mut induction = Array2::zeros((n_findex, ctx.n_turbines()));
            for group in unique_groups(&selection) {
                let params = self.operation_model_parameters_with_tables(Some((
                    &tables.entries[group].wind_speeds,
                    &tables.entries[group].power,
                    &tables.entries[group].thrust_coefficient,
                )))?;
                let group_induction = self
                    .turbine_type
                    .operation_model
                    .axial_induction(&params, &ctx)?;
                for fi in 0..n_findex {
                    if selection[fi] != group {
                        continue;
                    }
                    for ti in 0..ctx.n_turbines() {
                        induction[[fi, ti]] = group_induction[[fi, ti]];
                    }
                }
            }
            return Ok(induction);
        }

        let params = self.operation_model_parameters()?;

        // Create air density array
        let air_density_array = Array::from_elem(n_findex, air_density);

        // Create turbine context
        let ctx = TurbineContext {
            velocities,
            air_density: &air_density_array,
            yaw_angles,
            tilt_angles,
            power_setpoints,
            turbulence_intensities,
            awc_modes,
            awc_amplitudes,
            cubature_weights,
            correct_cp_ct_for_tilt,
            average_method,
        };

        // Call the operation model's axial_induction method
        self.turbine_type
            .operation_model
            .axial_induction(&params, &ctx)
    }

    /// The multi-dimensional power/thrust tables of this turbine, if any.
    fn multidim_tables(&self) -> Option<&MultiDimConditionTables> {
        match &self.turbine_type.power_thrust_table.cp_ct_table {
            crate::core::turbines::cp_ct_table::CpCtTable::MultiDimensionalConditions(tables) => {
                Some(tables)
            }
            _ => None,
        }
    }

    /// Select the power/thrust table of every wind condition.
    fn multidim_selection(
        &self,
        tables: &MultiDimConditionTables,
        conditions: Option<&MultidimConditionMap>,
        n_findex: usize,
    ) -> crate::Result<Vec<usize>> {
        let conditions = conditions.ok_or_else(|| {
            anyhow::anyhow!(
                "multidim_condition must be provided if using multidimensional \
                 power/thrust_coefficient"
            )
        })?;
        tables
            .select_entries(conditions, n_findex)
            .map_err(|e| anyhow::anyhow!("{e}"))
    }

    /// Build the turbine context shared by the operation models.
    #[allow(clippy::too_many_arguments)]
    fn turbine_context<'a>(
        &'a self,
        velocities: &'a Array4,
        air_density: &'a Array1,
        yaw_angles: Option<&'a Array2>,
        tilt_angles: Option<&'a Array2>,
        power_setpoints: Option<&'a Array2>,
        turbulence_intensities: Option<&'a Array4>,
        awc_modes: Option<&'a ndarray::Array2<String>>,
        awc_amplitudes: Option<&'a Array2>,
        cubature_weights: Option<&'a Array2>,
        correct_cp_ct_for_tilt: Option<&'a ndarray::Array2<bool>>,
        average_method: AveragingMethod,
    ) -> TurbineContext<'a> {
        TurbineContext {
            velocities,
            air_density,
            yaw_angles,
            tilt_angles,
            power_setpoints,
            turbulence_intensities,
            awc_modes,
            awc_amplitudes,
            cubature_weights,
            correct_cp_ct_for_tilt,
            average_method,
        }
    }
}

/// The unique group indices of a selection, in ascending order.
fn unique_groups(selection: &[usize]) -> Vec<usize> {
    let mut groups: Vec<usize> = selection.to_vec();
    groups.sort_unstable();
    groups.dedup();
    groups
}
