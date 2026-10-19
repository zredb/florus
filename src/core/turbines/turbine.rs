/// Turbine model and operations
///
/// Corresponds to core/turbine/ module in Python implementation
use crate::{
    core::{
        rotor_effective_velocity,
        turbines::{
            cp_ct_table::TableConditions,
            operation_models::{TurbineContext, TurbineParameters},
            turbine_type::TurbineType,
        },
        AveragingMethod,
    },
    types::{Array2, Array4, Float},
};
use ndarray::Array;
use std::sync::Arc;

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
    fn operation_model_parameters(&self) -> crate::Result<TurbineParameters> {
        // Create lookup tables from cp_ct_table
        let wind_speeds = self
            .turbine_type
            .power_thrust_table
            .cp_ct_table
            .wind_speeds();

        // Extract CP and CT values from the table
        let (cp_values, ct_values) = match &self.turbine_type.power_thrust_table.cp_ct_table {
            crate::core::turbines::cp_ct_table::CpCtTable::OneDimensional(table) => {
                (&table.cp_values, &table.ct_values)
            }
            _ => {
                return Err(anyhow::anyhow!(
                    "Multi-dimensional Cp/Ct tables are not yet supported by the operation-model dispatch (Turbine::power/thrust_coefficient/axial_induction)"
                ));
            }
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
        // tilt correction. interp1d semantics: linear in range; below the
        // first wind speed -> 0.0; above the last -> the last tilt angle.
        let tilt_interp = if self.turbine_type.correct_cp_ct_for_tilt {
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
        } else {
            None
        };

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
    ) -> crate::Result<Array2> {
        let params = self.operation_model_parameters()?;

        // Create air density array
        let n_findex = velocities.shape()[0];
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
    ) -> crate::Result<Array2> {
        let params = self.operation_model_parameters()?;

        // Create air density array
        let n_findex = velocities.shape()[0];
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
    ) -> crate::Result<Array2> {
        let params = self.operation_model_parameters()?;

        // Create air density array
        let n_findex = velocities.shape()[0];
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
}
