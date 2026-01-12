//! Wind rose - aggregated wind statistics by direction/speed bins.
//!
//! WindRose represents a wind resource distribution with frequency tables
//! binned by wind direction and wind speed.

use crate::heterogeneous_map::{HeterogeneousInflowConfig, HeterogeneousMap, MultidimConditions};
use crate::types::{Array1, Array2, Float};
use crate::wind_data::traits::{TIParams, WindData};
use crate::Result;
use serde::{Deserialize, Serialize};

/// Wind rose - aggregated wind statistics by direction/speed bins
///
/// WindRose represents a wind resource distribution where wind conditions
/// are binned by wind direction and wind speed, with associated frequency
/// and turbulence intensity tables.
///
/// Corresponds to WindRose in Python FLORIS v4.6
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindRose {
    /// Wind directions for each bin [n_directions]
    pub wind_directions: Array1,
    /// Wind speeds for each bin [n_speeds]
    pub wind_speeds: Array1,
    /// Turbulence intensity table [n_directions, n_speeds]
    pub ti_table: Array2,
    /// Frequency table [n_directions, n_speeds], None = uniform
    pub freq_table: Option<Array2>,
    /// Value table [n_directions, n_speeds], None = unit value
    pub value_table: Option<Array2>,
    /// Flag to compute zero frequency occurrence
    pub compute_zero_freq_occurrence: bool,
    /// Heterogeneous map for spatial variation
    pub heterogeneous_map: Option<HeterogeneousMap>,
    /// Multidimensional conditions (e.g., wave period, significant wave height)
    pub multidim_conditions: Option<MultidimConditions>,
}

impl Default for WindRose {
    fn default() -> Self {
        Self {
            wind_directions: Array1::from_vec(vec![]),
            wind_speeds: Array1::from_vec(vec![]),
            ti_table: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
            freq_table: None,
            value_table: None,
            compute_zero_freq_occurrence: false,
            heterogeneous_map: None,
            multidim_conditions: None,
        }
    }
}

impl WindRose {
    /// Create a new WindRose
    ///
    /// # Arguments
    /// * `wind_directions` - Wind direction bins
    /// * `wind_speeds` - Wind speed bins
    /// * `ti_table` - Turbulence intensity table [n_directions, n_speeds]
    /// * `freq_table` - Optional frequency table [n_directions, n_speeds]
    /// * `value_table` - Optional value table [n_directions, n_speeds]
    pub fn new(
        wind_directions: Array1,
        wind_speeds: Array1,
        ti_table: Array2,
        freq_table: Option<Array2>,
        value_table: Option<Array2>,
    ) -> Result<Self> {
        let n_dir = wind_directions.len();
        let n_ws = wind_speeds.len();

        if ti_table.shape() != &[n_dir, n_ws] {
            anyhow::bail!("ti_table must have shape ({}, {})", n_dir, n_ws);
        }

        if let Some(ref freq) = freq_table {
            if freq.shape() != &[n_dir, n_ws] {
                anyhow::bail!("freq_table must have shape ({}, {})", n_dir, n_ws);
            }
        }

        if let Some(ref value) = value_table {
            if value.shape() != &[n_dir, n_ws] {
                anyhow::bail!("value_table must have shape ({}, {})", n_dir, n_ws);
            }
        }

        Ok(Self {
            wind_directions,
            wind_speeds,
            ti_table,
            freq_table,
            value_table,
            compute_zero_freq_occurrence: false,
            heterogeneous_map: None,
            multidim_conditions: None,
        })
    }

    /// Assign TI using a function of wind direction and wind speed
    pub fn assign_ti_using_wd_ws_function<F>(&mut self, func: F)
    where
        F: Fn(Float, Float) -> Float,
    {
        for i in 0..self.wind_directions.len() {
            for j in 0..self.wind_speeds.len() {
                self.ti_table[[i, j]] =
                    func(self.wind_directions[i], self.wind_speeds[j]).clamp(0.0, 1.0);
            }
        }
    }

    /// Assign TI using IEC method
    pub fn assign_ti_using_iec_method(&mut self, params: Option<TIParams>) {
        let params = params.unwrap_or_default();
        for i in 0..self.wind_directions.len() {
            for j in 0..self.wind_speeds.len() {
                self.ti_table[[i, j]] = params.calculate_ti(self.wind_speeds[j]);
            }
        }
    }

    /// Assign value using a function of wind direction and wind speed
    pub fn assign_value_using_wd_ws_function<F>(&mut self, func: F, normalize: bool)
    where
        F: Fn(Float, Float) -> Float,
    {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let mut value_table = Array2::from_shape_fn((n_dir, n_ws), |(i, j)| {
            func(self.wind_directions[i], self.wind_speeds[j])
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

    /// Downsample (aggregate) wind rose to fewer bins
    pub fn downsample(&self, _wd_step: Float, _ws_step: Float, _inplace: bool) -> Self {
        // Simplified implementation - in real code would aggregate bins
        self.clone()
    }

    /// Upsample (resample) wind rose to more bins
    pub fn upsample(
        &self,
        _wd_step: Float,
        _ws_step: Float,
        _method: &str,
        _inplace: bool,
    ) -> Self {
        // Simplified implementation - in real code would use interpolation
        self.clone()
    }

    /// Aggregate (downsample) wind rose to fewer bins - wrapper for downsample
    ///
    /// This is a convenience method that wraps downsample() for backwards compatibility
    /// with the Python FLORIS API.
    pub fn aggregate(&self, wd_step: Float, ws_step: Float, inplace: bool) -> Self {
        self.downsample(wd_step, ws_step, inplace)
    }

    /// Resample (upsample) wind rose to more bins using interpolation - wrapper for upsample
    ///
    /// This is a convenience method that wraps upsample() for backwards compatibility
    /// with the Python FLORIS API.
    pub fn resample_by_interpolation(
        &self,
        wd_step: Float,
        ws_step: Float,
        method: &str,
        inplace: bool,
    ) -> Self {
        self.upsample(wd_step, ws_step, method, inplace)
    }

    /// Unpack multidimensional inflow conditions
    ///
    /// Returns the multidimensional conditions filtered by non-zero frequency mask.
    pub fn unpack_multidim_conditions(&self) -> Option<MultidimConditions> {
        let conditions = match &self.multidim_conditions {
            Some(c) => c,
            None => return None,
        };

        // Build the non-zero frequency mask
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let n_conditions = n_dir * n_ws;

        let default_freq = Array2::from_elem((n_dir, n_ws), 1.0);
        let freq = self.freq_table.as_ref().unwrap_or(&default_freq);

        // Collect non-zero indices
        let nonzero_indices: Vec<usize> = (0..n_conditions)
            .filter(|&i| {
                let d = i / n_ws;
                let s = i % n_ws;
                freq[[d, s]] > 0.0
            })
            .collect();

        if nonzero_indices.is_empty() {
            return None;
        }

        // Check if all TP values are the same (scalar case)
        let tp_first = conditions.tp[nonzero_indices[0]];
        let all_same_tp = nonzero_indices.iter().all(|&i| {
            let d = i / n_ws;
            let s = i % n_ws;
            (conditions.tp[d * n_ws + s] - tp_first).abs() < 1e-10
        });

        // Check if all HS values are the same (if HS exists)
        let all_same_hs = if let Some(ref hs) = conditions.hs {
            let hs_first = hs[nonzero_indices[0]];
            nonzero_indices.iter().all(|&i| {
                let d = i / n_ws;
                let s = i % n_ws;
                (hs[d * n_ws + s] - hs_first).abs() < 1e-10
            })
        } else {
            true
        };

        // If all values are the same, return scalar conditions
        if all_same_tp && all_same_hs {
            return Some(MultidimConditions {
                tp: Array1::from_vec(vec![tp_first]),
                hs: if conditions.hs.is_some() {
                    Some(Array1::from_vec(vec![conditions.hs.as_ref().unwrap()[0]]))
                } else {
                    None
                },
            });
        }

        // Otherwise, filter by non-zero mask
        let tp_filtered: Vec<Float> = nonzero_indices
            .iter()
            .map(|&i| {
                let d = i / n_ws;
                let s = i % n_ws;
                conditions.tp[d * n_ws + s]
            })
            .collect();

        let hs_filtered: Option<Array1> = if let Some(ref hs) = conditions.hs {
            let hs_vec: Vec<Float> = nonzero_indices
                .iter()
                .map(|&i| {
                    let d = i / n_ws;
                    let s = i % n_ws;
                    hs[d * n_ws + s]
                })
                .collect();
            Some(Array1::from_vec(hs_vec))
        } else {
            None
        };

        Some(MultidimConditions {
            tp: Array1::from_vec(tp_filtered),
            hs: hs_filtered,
        })
    }

    /// Convert to time series
    pub fn to_time_series(&self) -> super::TimeSeries {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();

        let default_freq = Array2::from_elem((n_dir, n_ws), 1.0);
        let freq = self.freq_table.as_ref().unwrap_or(&default_freq);

        let mut wd_flat = Vec::new();
        let mut ws_flat = Vec::new();
        let mut ti_flat = Vec::new();
        let mut count = 0usize;

        for i in 0..n_dir {
            for j in 0..n_ws {
                if freq[[i, j]] > 0.0 {
                    wd_flat.push(self.wind_directions[i]);
                    ws_flat.push(self.wind_speeds[j]);
                    ti_flat.push(self.ti_table[[i, j]]);
                    count += 1;
                }
            }
        }

        super::TimeSeries {
            wind_directions: Array1::from_vec(wd_flat),
            wind_speeds: Array1::from_vec(ws_flat),
            turbulence_intensities: Array1::from_vec(ti_flat),
            values: Array1::from_vec(vec![1.0; count]),
        }
    }
}

impl WindData for WindRose {
    fn wind_speeds(&self) -> Array1 {
        self.wind_speeds.clone()
    }

    fn wind_directions(&self) -> Array1 {
        self.wind_directions.clone()
    }

    fn turbulence_intensities(&self) -> Array1 {
        // For WindRose, we can't return a single TI array since TI depends on both direction and speed
        // Return the TI from the first speed bin for each direction as a representative value
        let n_dir = self.wind_directions.len();
        if n_dir > 0 && self.ti_table.shape().len() >= 2 {
            Array1::from_iter((0..n_dir).map(|i| self.ti_table[[i, 0]]))
        } else {
            Array1::from_vec(vec![])
        }
    }

    fn n_conditions(&self) -> usize {
        self.wind_directions.len() * self.wind_speeds.len()
    }

    fn heterogeneous_inflow_config(&self) -> HeterogeneousInflowConfig {
        let n_conditions = self.n_conditions();
        let n_points = if let Some(ref het_map) = self.heterogeneous_map {
            het_map.x.len()
        } else {
            0
        };

        if n_points > 0 {
            // Use the heterogeneous_map if available
            let het_map = self.heterogeneous_map.as_ref().unwrap();
            het_map
                .get_heterogeneous_inflow_config(
                    self.wind_directions.clone(),
                    self.wind_speeds.clone(),
                )
                .unwrap_or_else(|_| HeterogeneousInflowConfig {
                    x: het_map.x.clone(),
                    y: het_map.y.clone(),
                    z: het_map.z.clone(),
                    wind_speeds: Some(self.wind_speeds.clone()),
                    wind_directions: Some(self.wind_directions.clone()),
                    speed_multipliers: Array2::zeros((n_conditions, n_points)),
                })
        } else {
            // No heterogeneous map, return empty config
            HeterogeneousInflowConfig {
                x: Array1::from_vec(vec![]),
                y: Array1::from_vec(vec![]),
                z: None,
                wind_speeds: Some(self.wind_speeds.clone()),
                wind_directions: Some(self.wind_directions.clone()),
                speed_multipliers: Array2::from_shape_vec((n_conditions, 0), vec![]).unwrap(),
            }
        }
    }

    fn set_layout(&mut self, _layout_x: &Option<Array1>, _layout_y: &Option<Array1>) {
        // WindRose doesn't support layout changes through this interface
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
                if freq[[i, j]] > 0.0 {
                    wind_directions.push(self.wind_directions[i]);
                    wind_speeds.push(self.wind_speeds[j]);
                    turbulence_intensities.push(self.ti_table[[i, j]]);
                    frequencies.push(freq[[i, j]]);
                    values.push(value[[i, j]]);
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
    use crate::types::Array1;

    #[test]
    fn test_wind_rose_creation() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0, 270.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0, 12.0]);
        let ti_table = Array2::from_elem((4, 3), 0.08);
        let freq = Array2::from_shape_vec(
            (4, 3),
            vec![0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        )
        .unwrap();

        let wr = WindRose::new(wd, ws, ti_table, Some(freq), None).unwrap();
        assert_eq!(wr.n_conditions(), 12);
    }

    #[test]
    fn test_wind_rose_assign_ti() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti_table = Array2::from_elem((2, 2), 0.06);

        let mut wr = WindRose::new(wd, ws, ti_table, None, None).unwrap();
        wr.assign_ti_using_iec_method(None);

        // All TI values should now be calculated by IEC
        for val in wr.ti_table.iter() {
            assert!(*val > 0.0 && *val < 1.0);
        }
    }

    #[test]
    fn test_wind_rose_to_time_series() {
        let wd = Array1::from_vec(vec![0.0, 180.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti_table = Array2::from_shape_vec((2, 2), vec![0.06, 0.07, 0.08, 0.09]).unwrap();
        let freq = Array2::from_shape_vec((2, 2), vec![0.25, 0.25, 0.25, 0.25]).unwrap();

        let wr = WindRose::new(wd, ws, ti_table, Some(freq), None).unwrap();
        let ts = wr.to_time_series();

        // Should have 4 conditions (2 dirs × 2 speeds)
        assert_eq!(ts.n_conditions(), 4);
    }
}
