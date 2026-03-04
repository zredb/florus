//! Wind rose with turbulence intensity as an additional dimension.
//!
//! WindTIRose extends WindRose by adding turbulence intensity as an explicit
//! dimension, allowing for more detailed characterization of wind resources.

use crate::heterogeneous_map::HeterogeneousInflowConfig;
use crate::types::{Array1, Array2, Array3, Float};
use crate::wind_data::traits::{TIParams, WindData};
use crate::Result;
use serde::{Deserialize, Serialize};

/// WindTIRose - Wind rose with TI as an additional dimension
///
/// WindTIRose extends the standard WindRose by adding turbulence intensity
/// as an explicit dimension, allowing for joint probability distributions
/// of wind direction, speed, and turbulence intensity.
///
/// Corresponds to WindTIRose in Python FLORIS v4.6
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindTIRose {
    /// Wind directions [n_directions]
    pub wind_directions: Array1,
    /// Wind speeds [n_speeds]
    pub wind_speeds: Array1,
    /// Turbulence intensities [n_tis]
    pub turbulence_intensities: Array1,
    /// TI table [n_directions, n_speeds, n_tis]
    pub ti_table: Array3,
    /// Frequency table [n_directions, n_speeds]
    pub freq_table: Option<Array2>,
    /// Value table [n_directions, n_speeds]
    pub value_table: Option<Array2>,
}

impl Default for WindTIRose {
    fn default() -> Self {
        Self {
            wind_directions: Array1::from_vec(vec![]),
            wind_speeds: Array1::from_vec(vec![]),
            turbulence_intensities: Array1::from_vec(vec![]),
            ti_table: Array3::from_shape_vec((0, 0, 0), vec![]).unwrap(),
            freq_table: None,
            value_table: None,
        }
    }
}

impl WindTIRose {
    /// Create a new WindTIRose
    pub fn new(
        wind_directions: Array1,
        wind_speeds: Array1,
        turbulence_intensities: Array1,
        ti_table: Array3,
        freq_table: Option<Array2>,
        value_table: Option<Array2>,
    ) -> Result<Self> {
        let n_dir = wind_directions.len();
        let n_ws = wind_speeds.len();
        let n_ti = turbulence_intensities.len();

        if ti_table.shape() != &[n_dir, n_ws, n_ti] {
            anyhow::bail!("ti_table must have shape ({}, {}, {})", n_dir, n_ws, n_ti);
        }

        Ok(Self {
            wind_directions,
            wind_speeds,
            turbulence_intensities,
            ti_table,
            freq_table,
            value_table,
        })
    }

    /// Assign value using a function of wind direction, wind speed, and TI
    pub fn assign_value_using_wd_ws_ti_function<F>(&mut self, func: F, normalize: bool)
    where
        F: Fn(Float, Float, Float) -> Float,
    {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let mut value_table = Array2::from_shape_fn((n_dir, n_ws), |(i, j)| {
            func(
                self.wind_directions[i],
                self.wind_speeds[j],
                self.turbulence_intensities[0],
            )
        });

        if normalize {
            let mean: Float = value_table.iter().sum::<Float>() / (n_dir * n_ws) as Float;
            if mean > 0.0 {
                for val in &mut value_table {
                    *val /= mean;
                }
            }
        }

        self.value_table = Some(value_table);
    }

    /// Assign value using piecewise linear function
    pub fn assign_value_piecewise_linear(
        &mut self,
        value_zero_ws: Float,
        ws_knee: Float,
        slope_1: Float,
        slope_2: Float,
        limit_to_zero: bool,
        normalize: bool,
    ) {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let mut value_table = Array2::from_shape_fn((n_dir, n_ws), |(_, j)| {
            let ws = self.wind_speeds[j];
            if ws <= ws_knee {
                value_zero_ws + slope_1 * (ws - 3.0)
            } else {
                value_zero_ws + slope_1 * (ws_knee - 3.0) + slope_2 * (ws - ws_knee)
            }
        });

        if limit_to_zero {
            for val in &mut value_table {
                *val = val.max(0.0);
            }
        }

        if normalize {
            let mean: Float = value_table.iter().sum::<Float>() / (n_dir * n_ws) as Float;
            if mean > 0.0 {
                for val in &mut value_table {
                    *val /= mean;
                }
            }
        }

        self.value_table = Some(value_table);
    }

    /// Assign TI using IEC method
    pub fn assign_ti_using_iec_method(&mut self, params: Option<TIParams>) {
        let params = params.unwrap_or_default();
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let n_ti = self.turbulence_intensities.len();

        for i in 0..n_dir {
            for j in 0..n_ws {
                for k in 0..n_ti {
                    self.ti_table[[i, j, k]] = params.calculate_ti(self.wind_speeds[j]);
                }
            }
        }
    }
}

impl WindData for WindTIRose {
    fn wind_speeds(&self) -> Array1 {
        self.wind_speeds.clone()
    }

    fn wind_directions(&self) -> Array1 {
        self.wind_directions.clone()
    }

    fn turbulence_intensities(&self) -> Array1 {
        self.turbulence_intensities.clone()
    }

    fn frequencies(&self) -> Array2 {
        // Use unpack_freq to get frequency table
        self.unpack_freq()
    }

    fn heterogeneous_inflow_config(&self) -> HeterogeneousInflowConfig {
        let n_conditions = self.n_conditions();
        HeterogeneousInflowConfig {
            x: Array1::from_vec(vec![]),
            y: Array1::from_vec(vec![]),
            z: None,
            wind_speeds: Some(self.wind_speeds.clone()),
            wind_directions: Some(self.wind_directions.clone()),
            speed_multipliers: Array2::from_shape_vec((n_conditions, 0), vec![]).unwrap(),
        }
    }

    fn set_layout(&mut self, _layout_x: &Option<Array1>, _layout_y: &Option<Array1>) {
        // WindTIRose doesn't support layout changes
    }

    fn unpack(
        &self,
    ) -> (
        Array1,
        Array1,
        Array1,
        Array2,
        Array2,
        HeterogeneousInflowConfig,
    ) {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let n_ti = self.turbulence_intensities.len();

        let default_freq = Array2::from_elem((n_dir, n_ws), 1.0);
        let default_value = Array2::from_elem((n_dir, n_ws), 1.0);
        let freq = self.freq_table.as_ref().unwrap_or(&default_freq);
        let value = self.value_table.as_ref().unwrap_or(&default_value);

        let mut wind_directions = Vec::new();
        let mut wind_speeds = Vec::new();
        let mut turbulence_intensities = Vec::new();
        let mut frequencies = Vec::new();
        let mut values = Vec::new();

        for i in 0..n_dir {
            for j in 0..n_ws {
                for k in 0..n_ti {
                    if freq[[i, j]] > 0.0 {
                        wind_directions.push(self.wind_directions[i]);
                        wind_speeds.push(self.wind_speeds[j]);
                        turbulence_intensities.push(self.ti_table[[i, j, k]]);
                        frequencies.push(freq[[i, j]]);
                        values.push(value[[i, j]]);
                    }
                }
            }
        }

        // Convert to 2D arrays with single column for single turbine
        let n = frequencies.len();
        let freq_2d = Array2::from_shape_vec((n, 1), frequencies).unwrap();
        let value_2d = Array2::from_shape_vec((n, 1), values).unwrap();

        (
            Array1::from_vec(wind_directions),
            Array1::from_vec(wind_speeds),
            Array1::from_vec(turbulence_intensities),
            freq_2d,
            value_2d,
            self.heterogeneous_inflow_config(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Array1, Array2, Array3};
    use crate::wind_data::traits::WindData;

    // ============================================================================
    // WindTIRose Creation Tests
    // ============================================================================

    #[test]
    fn test_wind_ti_rose_creation_basic() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0, 270.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0, 15.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08, 0.10]);
        let ti_table = Array3::from_elem((4, 3, 3), 0.08);

        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        assert_eq!(wtr.wind_directions.len(), 4);
        assert_eq!(wtr.wind_speeds.len(), 3);
        assert_eq!(wtr.turbulence_intensities.len(), 3);
    }

    #[test]
    fn test_wind_ti_rose_creation_with_freq_and_value() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 12.0]);
        let ti = Array1::from_vec(vec![0.06, 0.10]);
        let ti_table = Array3::from_elem((2, 2, 2), 0.08);
        let freq = Array2::from_shape_vec((2, 2), vec![0.2, 0.3, 0.2, 0.3]).unwrap();
        let value = Array2::from_shape_vec((2, 2), vec![1.0, 2.0, 3.0, 4.0]).unwrap();

        let wtr = WindTIRose::new(wd, ws, ti, ti_table, Some(freq), Some(value)).unwrap();

        assert!(wtr.freq_table.is_some());
        assert!(wtr.value_table.is_some());

        let freq = wtr.freq_table.unwrap();
        let value = wtr.value_table.unwrap();
        assert_eq!(freq.shape(), &[2, 2]);
        assert_eq!(value.shape(), &[2, 2]);
    }

    #[test]
    fn test_wind_ti_rose_default() {
        let wtr = WindTIRose::default();

        assert!(wtr.wind_directions.is_empty());
        assert!(wtr.wind_speeds.is_empty());
        assert!(wtr.turbulence_intensities.is_empty());
        assert_eq!(wtr.ti_table.shape(), &[0, 0, 0]);
        assert!(wtr.freq_table.is_none());
        assert!(wtr.value_table.is_none());
    }

    #[test]
    fn test_wind_ti_rose_creation_invalid_ti_shape() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08]);
        let ti_table = Array3::from_elem((2, 2, 2), 0.08); // Wrong shape

        let result = WindTIRose::new(wd, ws, ti, ti_table, None, None);

        assert!(result.is_err());
    }

    #[test]
    fn test_wind_ti_rose_creation_single_ti() {
        let wd = Array1::from_vec(vec![0.0, 180.0]);
        let ws = Array1::from_vec(vec![10.0]);
        let ti = Array1::from_vec(vec![0.08]);
        let ti_table = Array3::from_shape_vec((2, 1, 1), vec![0.08, 0.08]).unwrap();

        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        assert_eq!(wtr.turbulence_intensities.len(), 1);
    }

    // ============================================================================
    // Value Assignment Tests
    // ============================================================================

    #[test]
    fn test_assign_value_using_wd_ws_ti_function() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.10]);
        let ti_table = Array3::from_elem((3, 2, 2), 0.08);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        // Assign value based on wind speed only (TI not used in value calculation)
        wtr.assign_value_using_wd_ws_ti_function(|_wd, ws, _ti| ws / 10.0, false);

        assert!(wtr.value_table.is_some());
        let vt = wtr.value_table.unwrap();
        // At ws=5.0, value = 0.5
        assert!((vt[[0, 0]] - 0.5).abs() < 1e-10);
        // At ws=10.0, value = 1.0
        assert!((vt[[0, 1]] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_assign_value_using_wd_ws_ti_function_with_normalize() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 2, 1), 0.08);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        // All values are 2.0, after normalization they should be 1.0
        wtr.assign_value_using_wd_ws_ti_function(|_wd, _ws, _ti| 2.0, true);

        let vt = wtr.value_table.unwrap();
        for val in vt.iter() {
            assert!((val - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_assign_value_using_wd_ws_ti_function_complex() {
        let wd = Array1::from_vec(vec![0.0, 180.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0, 15.0]);
        let ti = Array1::from_vec(vec![0.05, 0.10, 0.15]);
        let ti_table = Array3::from_elem((2, 3, 3), 0.08);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        // Value = wd/100 + ws/10 + ti
        wtr.assign_value_using_wd_ws_ti_function(|wd, ws, ti| wd / 100.0 + ws / 10.0 + ti, false);

        let vt = wtr.value_table.unwrap();
        // At wd=0, ws=5, ti (first) = 0.05: value = 0 + 0.5 + 0.05 = 0.55
        assert!((vt[[0, 0]] - 0.55).abs() < 1e-10);
        // At wd=180, ws=15, ti (first) = 0.05: value = 1.8 + 1.5 + 0.05 = 3.35
        assert!((vt[[1, 2]] - 3.35).abs() < 1e-10);
    }

    #[test]
    fn test_assign_value_piecewise_linear_basic() {
        let wd = Array1::from_vec(vec![0.0]);
        let ws = Array1::from_vec(vec![3.0, 6.0, 9.0, 12.0]);
        let ti = Array1::from_vec(vec![0.08]);
        let ti_table = Array3::from_elem((1, 4, 1), 0.08);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        // value_zero_ws = 1.0, ws_knee = 8.0, slope_1 = 0.1, slope_2 = 0.05
        wtr.assign_value_piecewise_linear(1.0, 8.0, 0.1, 0.05, false, false);

        let vt = wtr.value_table.unwrap();
        // At ws=3.0: value = 1.0 + 0.1 * (3 - 3) = 1.0
        assert!((vt[[0, 0]] - 1.0).abs() < 1e-10);
        // At ws=6.0: value = 1.0 + 0.1 * (6 - 3) = 1.3
        assert!((vt[[0, 1]] - 1.3).abs() < 1e-10);
        // At ws=9.0: value = 1.0 + 0.1 * (8 - 3) + 0.05 * (9 - 8) = 1.55
        assert!((vt[[0, 2]] - 1.55).abs() < 1e-10);
        // At ws=12.0: value = 1.0 + 0.1 * (8 - 3) + 0.05 * (12 - 8) = 1.7
        assert!((vt[[0, 3]] - 1.7).abs() < 1e-10);
    }

    #[test]
    fn test_assign_value_piecewise_linear_limit_to_zero() {
        let wd = Array1::from_vec(vec![0.0]);
        let ws = Array1::from_vec(vec![3.0, 4.0, 5.0]);
        let ti = Array1::from_vec(vec![0.08]);
        let ti_table = Array3::from_elem((1, 3, 1), 0.08);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        // Negative slope that would produce negative values
        wtr.assign_value_piecewise_linear(0.0, 10.0, -0.5, -0.1, true, false);

        let vt = wtr.value_table.unwrap();
        // All values should be >= 0
        for val in vt.iter() {
            assert!(*val >= 0.0);
        }
    }

    #[test]
    fn test_assign_value_piecewise_linear_with_normalize() {
        let wd = Array1::from_vec(vec![0.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0]);
        let ti = Array1::from_vec(vec![0.08]);
        let ti_table = Array3::from_elem((1, 2, 1), 0.08);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        wtr.assign_value_piecewise_linear(1.0, 8.0, 0.1, 0.05, false, true);

        let vt = wtr.value_table.unwrap();
        let mean: Float = vt.iter().sum::<Float>() / 2.0;
        // After normalization, mean should be approximately 1.0
        assert!((mean - 1.0).abs() < 1e-10);
    }

    // ============================================================================
    // TI Assignment Tests
    // ============================================================================

    #[test]
    fn test_assign_ti_using_iec_method_default_params() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 12.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08]);
        let ti_table = Array3::from_elem((2, 2, 2), 0.0);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        wtr.assign_ti_using_iec_method(None);

        // IEC method should produce valid TI values
        for val in wtr.ti_table.iter() {
            assert!(*val > 0.0 && *val <= 1.0);
        }
    }

    #[test]
    fn test_assign_ti_using_iec_method_consistent_per_speed() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0]);
        let ws = Array1::from_vec(vec![8.0, 12.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08, 0.10]);
        let ti_table = Array3::from_elem((3, 2, 3), 0.0);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        wtr.assign_ti_using_iec_method(None);

        // TI should be the same for all directions and TI bins at the same wind speed
        // Check that all values for a given wind speed are equal
        let ti_at_ws_8 = wtr.ti_table[[0, 0, 0]];
        let ti_at_ws_12 = wtr.ti_table[[0, 1, 0]];

        for i in 0..3 {
            for k in 0..3 {
                assert!((wtr.ti_table[[i, 0, k]] - ti_at_ws_8).abs() < 1e-10);
                assert!((wtr.ti_table[[i, 1, k]] - ti_at_ws_12).abs() < 1e-10);
            }
        }
    }

    // ============================================================================
    // WindData Trait Tests
    // ============================================================================

    #[test]
    fn test_wind_data_trait_wind_speeds() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0, 15.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 3, 1), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        let speeds = wtr.wind_speeds();
        assert_eq!(speeds.len(), 3);
        assert!((speeds[0] - 5.0).abs() < 1e-10);
        assert!((speeds[2] - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_wind_data_trait_wind_directions() {
        let wd = Array1::from_vec(vec![0.0, 45.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((3, 2, 1), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        let directions = wtr.wind_directions();
        assert_eq!(directions.len(), 3);
        assert!((directions[0] - 0.0).abs() < 1e-10);
        assert!((directions[2] - 90.0).abs() < 1e-10);
    }

    #[test]
    fn test_wind_data_trait_turbulence_intensities() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.05, 0.08, 0.12]);
        let ti_table = Array3::from_elem((2, 2, 3), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        let tis = wtr.turbulence_intensities();
        assert_eq!(tis.len(), 3);
        assert!((tis[0] - 0.05).abs() < 1e-10);
        assert!((tis[2] - 0.12).abs() < 1e-10);
    }

    #[test]
    fn test_wind_data_trait_n_conditions() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08]);
        let ti_table = Array3::from_elem((3, 2, 2), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        // n_conditions = 3 * 2 * 2 = 12
        assert_eq!(wtr.n_conditions(), 12);
    }

    #[test]
    fn test_wind_data_trait_frequencies() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 2, 1), 0.08);
        let freq = Array2::from_shape_vec((2, 2), vec![0.1, 0.2, 0.3, 0.4]).unwrap();
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, Some(freq), None).unwrap();

        let frequencies = wtr.frequencies();
        assert_eq!(frequencies.shape(), &[wtr.n_conditions(), 1]);
    }

    #[test]
    fn test_wind_data_trait_heterogeneous_inflow_config() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 2, 1), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        let config = wtr.heterogeneous_inflow_config();

        // Should have empty points since no heterogeneous map
        assert!(config.x.is_empty());
        assert!(config.y.is_empty());
    }

    #[test]
    fn test_wind_data_trait_set_layout() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 2, 1), 0.08);
        let mut wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        // set_layout should be a no-op for WindTIRose
        let layout_x = Some(Array1::from_vec(vec![0.0, 100.0]));
        let layout_y = Some(Array1::from_vec(vec![0.0, 100.0]));
        wtr.set_layout(&layout_x, &layout_y);

        // Data should remain unchanged
        assert_eq!(wtr.wind_directions.len(), 2);
    }

    #[test]
    fn test_wind_data_trait_unpack() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.10]);
        let ti_table = Array3::from_shape_vec(
            (2, 2, 2),
            vec![0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12],
        )
        .unwrap();
        let freq = Array2::from_shape_vec((2, 2), vec![0.2, 0.3, 0.2, 0.3]).unwrap();
        let value = Array2::from_shape_vec((2, 2), vec![1.0, 2.0, 3.0, 4.0]).unwrap();
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, Some(freq), Some(value)).unwrap();

        let (wd_unpack, ws_unpack, ti_unpack, freq_2d, value_2d, _het_config) = wtr.unpack();

        // Should have 2 * 2 * 2 = 8 conditions (all with freq > 0)
        assert_eq!(wd_unpack.len(), 8);
        assert_eq!(ws_unpack.len(), 8);
        assert_eq!(ti_unpack.len(), 8);
        assert_eq!(freq_2d.shape(), &[8, 1]);
        assert_eq!(value_2d.shape(), &[8, 1]);
    }

    #[test]
    fn test_wind_data_trait_unpack_with_zero_freq() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 2, 1), 0.08);
        // One bin has zero frequency
        let freq = Array2::from_shape_vec((2, 2), vec![0.3, 0.0, 0.3, 0.4]).unwrap();
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, Some(freq), None).unwrap();

        let (wd_unpack, ws_unpack, ti_unpack, freq_2d, _value_2d, _het_config) = wtr.unpack();

        // Should have 3 conditions (one zero freq bin excluded)
        // Actually for WindTIRose, it's 3 * 1 = 3 conditions
        assert_eq!(wd_unpack.len(), 3);
    }

    // ============================================================================
    // Edge Cases Tests
    // ============================================================================

    #[test]
    fn test_wind_ti_rose_single_direction() {
        let wd = Array1::from_vec(vec![180.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0, 15.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08]);
        let ti_table = Array3::from_elem((1, 3, 2), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        assert_eq!(wtr.wind_directions.len(), 1);
        assert_eq!(wtr.n_conditions(), 6);
    }

    #[test]
    fn test_wind_ti_rose_single_speed() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0, 270.0]);
        let ws = Array1::from_vec(vec![10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08, 0.10]);
        let ti_table = Array3::from_elem((4, 1, 3), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        assert_eq!(wtr.wind_speeds.len(), 1);
        assert_eq!(wtr.n_conditions(), 12);
    }

    #[test]
    fn test_wind_ti_rose_single_ti() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0]);
        let ws = Array1::from_vec(vec![5.0, 10.0]);
        let ti = Array1::from_vec(vec![0.08]);
        let ti_table = Array3::from_elem((3, 2, 1), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        assert_eq!(wtr.turbulence_intensities.len(), 1);
        assert_eq!(wtr.n_conditions(), 6);
    }

    #[test]
    fn test_wind_ti_rose_large_data() {
        // Create a large wind TI rose with many bins
        let n_dir = 36; // 10-degree resolution
        let n_ws = 25; // 1 m/s from 0-25
        let n_ti = 5; // 5 TI bins

        let wd: Vec<Float> = (0..n_dir).map(|i| i as Float * 10.0).collect();
        let ws: Vec<Float> = (0..n_ws).map(|i| i as Float).collect();
        let ti: Vec<Float> = (0..n_ti).map(|i| 0.05 + i as Float * 0.02).collect();

        let ti_table = Array3::from_elem((n_dir, n_ws, n_ti), 0.08);

        let wtr = WindTIRose::new(
            Array1::from_vec(wd),
            Array1::from_vec(ws),
            Array1::from_vec(ti),
            ti_table,
            None,
            None,
        )
        .unwrap();

        assert_eq!(wtr.n_conditions(), n_dir * n_ws * n_ti);
    }

    // ============================================================================
    // Serialization Tests
    // ============================================================================

    #[test]
    fn test_wind_ti_rose_serialization() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.10]);
        let ti_table = Array3::from_shape_vec(
            (2, 2, 2),
            vec![0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12],
        )
        .unwrap();
        let freq = Array2::from_shape_vec((2, 2), vec![0.25; 4]).unwrap();
        let value = Array2::from_shape_vec((2, 2), vec![1.0, 2.0, 3.0, 4.0]).unwrap();

        let original = WindTIRose::new(wd, ws, ti, ti_table, Some(freq), Some(value)).unwrap();

        // Serialize
        let json = serde_json::to_string(&original).unwrap();

        // Deserialize
        let deserialized: WindTIRose = serde_json::from_str(&json).unwrap();

        assert_eq!(
            deserialized.wind_directions.len(),
            original.wind_directions.len()
        );
        assert_eq!(deserialized.wind_speeds.len(), original.wind_speeds.len());
        assert_eq!(
            deserialized.turbulence_intensities.len(),
            original.turbulence_intensities.len()
        );
        assert_eq!(deserialized.ti_table.shape(), original.ti_table.shape());
        assert!(deserialized.freq_table.is_some());
        assert!(deserialized.value_table.is_some());
    }

    #[test]
    fn test_wind_ti_rose_clone() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 2, 1), 0.08);
        let freq = Array2::from_elem((2, 2), 0.25);

        let original = WindTIRose::new(wd, ws, ti, ti_table, Some(freq), None).unwrap();
        let cloned = original.clone();

        assert_eq!(original.wind_directions.len(), cloned.wind_directions.len());
        assert_eq!(original.wind_speeds.len(), cloned.wind_speeds.len());
        assert_eq!(
            original.turbulence_intensities.len(),
            cloned.turbulence_intensities.len()
        );
        assert_eq!(original.ti_table.shape(), cloned.ti_table.shape());
    }

    // ============================================================================
    // Debug Test
    // ============================================================================

    #[test]
    fn test_wind_ti_rose_debug() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0]);
        let ti = Array1::from_vec(vec![0.06]);
        let ti_table = Array3::from_elem((2, 1, 1), 0.08);
        let wtr = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();

        let debug_str = format!("{:?}", wtr);
        assert!(debug_str.contains("WindTIRose"));
        assert!(debug_str.contains("wind_directions"));
        assert!(debug_str.contains("wind_speeds"));
        assert!(debug_str.contains("turbulence_intensities"));
    }
}
