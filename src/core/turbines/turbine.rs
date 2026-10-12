/// Turbine model and operations
///
/// Corresponds to core/turbine/ module in Python implementation
use crate::{
    core::{
        rotor_effective_velocity,
        turbines::{cp_ct_table::TableConditions, turbine_type::TurbineType},
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
    /// Calculate power output for given velocities
    ///
    /// # Arguments
    /// * `velocities` - Wind velocities at turbine rotor (findex, n_turbines, n_points, n_grid)
    /// * `air_density` - Air density [kg/m^3]
    /// * `yaw_angles` - Yaw angles [degrees]
    /// * `tilt_angles` - Tilt angles [degrees]
    /// * `average_method` - Method for averaging velocities
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
    ) -> crate::Result<Array2> {
        use crate::core::turbines::operation_models::{
            OperationModel, TurbineContext, TurbineParameters,
        };

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
                // For multi-dimensional tables, we need to handle differently
                // For now, return empty arrays (this should be improved)
                return Err(anyhow::anyhow!(
                    "Multi-dimensional Cp/Ct tables not yet supported in Turbine::power"
                ));
            }
        };

        let power_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: cp_values.clone(),
            fill_value: 0.0,
        };

        let thrust_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: ct_values.clone(),
            fill_value: 0.0001,
        };

        // Get ref_tilt value
        let ref_tilt = self.turbine_type.power_thrust_table.ref_tilt.unwrap_or(5.0);

        // Create turbine parameters
        let params = TurbineParameters {
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
        };

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
            cubature_weights: None,
            correct_cp_ct_for_tilt: None,
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
    ) -> crate::Result<Array2> {
        use crate::core::turbines::operation_models::{
            OperationModel, TurbineContext, TurbineParameters,
        };

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
                // For multi-dimensional tables, we need to handle differently
                return Err(anyhow::anyhow!("Multi-dimensional Cp/Ct tables not yet supported in Turbine::thrust_coefficient"));
            }
        };

        let power_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: cp_values.clone(),
            fill_value: 0.0,
        };

        let thrust_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: ct_values.clone(),
            fill_value: 0.0001,
        };

        // Get ref_tilt value
        let ref_tilt = self.turbine_type.power_thrust_table.ref_tilt.unwrap_or(5.0);

        // Create turbine parameters
        let params = TurbineParameters {
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
        };

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
            cubature_weights: None,
            correct_cp_ct_for_tilt: None,
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
    ) -> crate::Result<Array2> {
        use crate::core::turbines::operation_models::{
            OperationModel, TurbineContext, TurbineParameters,
        };

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
                // For multi-dimensional tables, we need to handle differently
                return Err(anyhow::anyhow!(
                    "Multi-dimensional Cp/Ct tables not yet supported in Turbine::axial_induction"
                ));
            }
        };

        let power_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: cp_values.clone(),
            fill_value: 0.0,
        };

        let thrust_table = crate::core::turbines::operation_models::LookupTable {
            keys: ndarray::Array1::from_vec(wind_speeds.to_vec()),
            values: ct_values.clone(),
            fill_value: 0.0001,
        };

        // Get ref_tilt value
        let ref_tilt = self.turbine_type.power_thrust_table.ref_tilt.unwrap_or(5.0);

        // Create turbine parameters
        let params = TurbineParameters {
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
        };

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
            cubature_weights: None,
            correct_cp_ct_for_tilt: None,
            average_method,
        };

        // Call the operation model's axial_induction method
        self.turbine_type
            .operation_model
            .axial_induction(&params, &ctx)
    }
}
