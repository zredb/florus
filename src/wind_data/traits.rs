//checked by zhb on 2026-1-16

//! WindData trait and common types for wind data structures.
//!
//! This module defines the base trait that all wind data sources must implement.
//!
use crate::heterogeneous_map::{HeterogeneousInflowConfig, MultidimConditions};
use crate::types::{Array1, Array2, Float};

/// Base trait for wind data sources
///
/// This trait defines the common interface that all wind data structures
/// must implement to be used with FLORIS for wind farm simulations.
pub trait WindData {
    /// Get wind speeds
    fn wind_speeds(&self) -> Array1;

    /// Get wind directions
    fn wind_directions(&self) -> Array1;

    /// Get turbulence intensities
    fn turbulence_intensities(&self) -> Array1;

    /// Get number of conditions
    fn n_conditions(&self) -> usize {
        self.wind_directions().len()
            * self.wind_speeds().len()
            * self.turbulence_intensities().len()
    }

    /// Get frequencies for AEP calculations
    ///
    /// Returns frequency table [n_conditions, n_turbines]
    fn frequencies(&self) -> Array2;

    /// Set heterogeneous inflow configuration
    fn heterogeneous_inflow_config(&self) -> HeterogeneousInflowConfig;

    /// Set turbine layout
    fn set_layout(&mut self, layout_x: &Option<Array1>, layout_y: &Option<Array1>);

    /// Unpack wind conditions for simulation
    ///
    /// Returns:
    /// - wind_directions: Array1
    /// - wind_speeds: Array1
    /// - turbulence_intensities: Array1
    /// - frequency table [n_conditions, n_turbines]
    /// - value table [n_conditions, n_turbines]
    /// - heterogeneous_inflow_config
    fn unpack(
        &self,
    ) -> (
        Array1,
        Array1,
        Array1,
        Array2,
        Array2,
        HeterogeneousInflowConfig,
    );

    /// Unpack for reinitialization
    ///
    /// Provides wind conditions without frequency/value tables.
    fn unpack_for_reinitialize(&self) -> (Array1, Array1, Array1, HeterogeneousInflowConfig) {
        let (
            wind_directions_unpack,
            wind_speeds_unpack,
            ti_table_unpack,
            _,
            _,
            heterogeneous_inflow_config,
        ) = self.unpack().clone();
        (
            wind_directions_unpack,
            wind_speeds_unpack,
            ti_table_unpack,
            heterogeneous_inflow_config,
        )
    }

    fn unpack_freq(&self) -> Array2 {
        let (_, _, _, frequency_table, _, _) = self.unpack().clone();
        frequency_table
    }
    fn unpack_value(&self) -> Array2 {
        let (_, _, _, _, value_table, _) = self.unpack().clone();
        value_table
    }
    fn unpack_multidim_conditions(&self) -> MultidimConditions {
        //   NOTE: This is a temporary method for backwards compatibility and will be removed in a
        // future release, when multidim_conditions are included in the unpack() method of child
        // classes.
        unimplemented!()
    }
    // rust 的类型系统可以保证heterogeneous_inflow_config 正确, 所以不需要check方法
    // fn check_heterogeneous_inflow_config(&self, heterogeneous_inflow_config: HeterogeneousInflowConfig) -> anyhow::Result<()> {
    //     Ok(())
    // }
}

/// Turbulence intensity parameters for IEC method
///
/// Used to calculate TI based on wind speed using the IEC 61400-1
/// normal turbulence model.
#[derive(Debug, Clone, Copy)]
pub struct TIParams {
    /// Reference turbulence level at 15 m/s
    pub iref: Float,
    /// Offset value from IEC standard
    pub offset: Float,
}

impl Default for TIParams {
    fn default() -> Self {
        Self {
            iref: 0.07,  // Default Iref (lower than IEC classes for realistic TI values)
            offset: 3.8, // IEC standard offset
        }
    }
}

impl TIParams {
    /// Calculate TI using IEC method
    pub fn calculate_ti(&self, wind_speed: Float) -> Float {
        // IEC 61400-1 normal turbulence model
        // TI = iref * (15 + offset) / (wind_speed + offset)
        let ti = self.iref * (15.0 + self.offset) / (wind_speed + self.offset);
        ti.clamp(0.0, 1.0)
    }
}
