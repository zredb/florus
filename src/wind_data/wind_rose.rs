//! Wind rose - aggregated wind statistics by direction/speed bins.
//!
//! WindRose represents a wind resource distribution with frequency tables
//! binned by wind direction and wind speed.

use crate::core::InterpMethod;
use crate::heterogeneous_map::{HeterogeneousInflowConfig, HeterogeneousMap, MultidimConditions};
use crate::types::{Array1, Array2, Float};
use crate::wind_data::traits::{TIParams, WindData};
use crate::wind_data::{ValidationError, ValidationResult};
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
    pub freq_table: Array2,
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
            freq_table: Array2::ones((0, 0)),
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
        compute_zero_freq_occurrence: bool,
        heterogeneous_map: Option<HeterogeneousMap>,
        multidim_conditions: Option<MultidimConditions>,
    ) -> ValidationResult<Self> {
        let n_dir = wind_directions.len();
        let n_ws = wind_speeds.len();

        crate::wind_data::validate_wind_arrays(&wind_directions.view(), &wind_speeds.view())?;

        if ti_table.shape() != &[n_dir, n_ws] {
            return Err(ValidationError::InvalidShape2(n_dir, n_ws));
        }

        if let Some(ref freq) = freq_table {
            if freq.shape() != &[n_dir, n_ws] {
                return Err(ValidationError::InvalidShape2(n_dir, n_ws));
            }
        }

        if let Some(ref value) = value_table {
            if value.shape() != &[n_dir, n_ws] {
                return Err(ValidationError::InvalidShape2(n_dir, n_ws));
            }
        }
        let freq_table = match freq_table {
            Some(f) => f,
            None => Array2::from_elem((n_dir, n_ws), 1.0),
        };
        let mut freq_table = freq_table;
        // Normalize frequency table to sum to 1.0

        let freq_sum = freq_table.sum();
        if freq_sum > 0.0 {
            freq_table = &freq_table / freq_sum;
        }

        Ok(Self {
            wind_directions,
            wind_speeds,
            ti_table,
            freq_table,
            value_table,
            compute_zero_freq_occurrence,
            heterogeneous_map,
            multidim_conditions,
        })
    }
    fn get_wd_grid(&self) -> Array2 {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        Array2::from_shape_fn((n_dir, n_ws), |(i, _)| self.wind_directions[i])
    }
    fn get_ws_grid(&self) -> Array2 {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        Array2::from_shape_fn((n_dir, n_ws), |(_, j)| self.wind_speeds[j])
    }
    fn get_wd_flat(&self) -> Array1 {
        self.get_wd_grid().flatten().to_owned()
    }
    fn get_ws_flat(&self) -> Array1 {
        self.get_ws_grid().flatten().to_owned()
    }

    fn get_ti_flat(&self) -> Array1 {
        self.ti_table.flatten().to_owned()
    }
    fn get_freq_flat(&self) -> Array1 {
        self.freq_table.flatten().to_owned()
    }
    fn get_value_flat(&self) -> Option<Array1> {
        self.value_table.as_ref().map(|v| v.flatten().to_owned())
    }

    #[allow(dead_code)]
    fn get_multidim_conditions_flat(&self) -> Option<MultidimConditions> {
        self.multidim_conditions.as_ref().map(|mc| {
            let n_dir = self.wind_directions.len();
            let n_ws = self.wind_speeds.len();
            let n_conditions = n_dir * n_ws;

            let tp_flat: Vec<Float> = (0..n_conditions)
                .map(|i| {
                    let d = i / n_ws;
                    let s = i % n_ws;
                    mc.tp[d * n_ws + s]
                })
                .collect();

            let hs_flat: Option<Array1> = if let Some(ref hs) = mc.hs {
                let hs_vec: Vec<Float> = (0..n_conditions)
                    .map(|i| {
                        let d = i / n_ws;
                        let s = i % n_ws;
                        hs[d * n_ws + s]
                    })
                    .collect();
                Some(Array1::from_vec(hs_vec))
            } else {
                None
            };

            MultidimConditions {
                tp: Array1::from_vec(tp_flat),
                hs: hs_flat,
            }
        })
    }
    #[allow(dead_code)]
    fn n_conditions(&self) -> usize {
        self.wind_directions.len() * self.wind_speeds.len()
    }
    #[allow(dead_code)]
    fn get_non_zero_freq_mask(&self) -> Vec<usize> {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let n_conditions = n_dir * n_ws;

        (0..n_conditions)
            .filter(|&i| {
                let d = i / n_ws;
                let s = i % n_ws;
                self.freq_table[[d, s]] > 0.0
            })
            .collect()
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
    pub fn downsample(
        &self,
        wd_step: Option<Float>,
        ws_step: Option<Float>,
        _ti_step: Option<Float>,
    ) -> Self {
        if self.wind_directions.is_empty() || self.wind_speeds.is_empty() {
            return self.clone();
        }

        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();

        let ws_step_current = if n_ws >= 2 {
            self.wind_speeds[1] - self.wind_speeds[0]
        } else {
            1.0
        };

        let wd_step_current = if n_dir >= 2 {
            self.wind_directions[1] - self.wind_directions[0]
        } else {
            360.0
        };

        let wd_step = wd_step.unwrap_or(wd_step_current);
        let ws_step = ws_step.unwrap_or(ws_step_current);

        if wd_step < wd_step_current {
            panic!("wd_step must be >= current step ({}).", wd_step_current);
        }

        if ws_step < ws_step_current {
            panic!("ws_step must be >= current step ({}).", ws_step_current);
        }

        let wd_min = self.wind_directions.iter().fold(Float::INFINITY, |acc, &v| acc.min(v));
        let wd_max = self.wind_directions.iter().fold(Float::NEG_INFINITY, |acc, &v| acc.max(v));
        let ws_min = self.wind_speeds.iter().fold(Float::INFINITY, |acc, &v| acc.min(v));
        let ws_max = self.wind_speeds.iter().fold(Float::NEG_INFINITY, |acc, &v| acc.max(v));

        let wd_range_min = wd_min - wd_step_current / 2.0;
        let wd_range_max = wd_max + wd_step_current / 2.0;
        let ws_range_min = ws_min - ws_step_current / 2.0;
        let ws_range_max = ws_max + ws_step_current / 2.0;

        let mut new_wind_directions = Vec::new();
        let mut wd_cursor = wd_range_min + wd_step / 2.0;
        while wd_cursor <= wd_range_max + 1e-9 {
            let wrapped = ((wd_cursor % 360.0) + 360.0) % 360.0;
            new_wind_directions.push(wrapped);
            wd_cursor += wd_step;
        }

        let mut new_wind_speeds = Vec::new();
        let mut ws_cursor = (ws_range_min + ws_step / 2.0).max(0.0);
        while ws_cursor <= ws_range_max + 1e-9 {
            new_wind_speeds.push(ws_cursor);
            ws_cursor += ws_step;
        }

        let new_shape = (new_wind_directions.len(), new_wind_speeds.len());
        let mut freq_acc = Array2::zeros(new_shape);
        let mut ti_sum = Array2::zeros(new_shape);
        let mut value_sum = self.value_table.as_ref().map(|_| Array2::zeros(new_shape));

        for i in 0..n_dir {
            for j in 0..n_ws {
                let freq = self.freq_table[[i, j]];
                if freq <= 0.0 {
                    continue;
                }

                let wd_idx = (((self.wind_directions[i] - wd_range_min) / wd_step).floor() as isize)
                    .max(0)
                    .min((new_shape.0 - 1) as isize) as usize;
                let ws_idx = (((self.wind_speeds[j] - ws_range_min) / ws_step).floor() as isize)
                    .max(0)
                    .min((new_shape.1 - 1) as isize) as usize;

                freq_acc[[wd_idx, ws_idx]] += freq;
                ti_sum[[wd_idx, ws_idx]] += freq * self.ti_table[[i, j]];
                if let (Some(vs), Some(value_table)) = (&mut value_sum, &self.value_table) {
                    vs[[wd_idx, ws_idx]] += freq * value_table[[i, j]];
                }
            }
        }

        let freq_raw_total: Float = freq_acc.iter().sum();
        let mut freq_table = freq_acc.clone();
        if freq_raw_total > 0.0 {
            for f in freq_table.iter_mut() {
                *f /= freq_raw_total;
            }
        }

        let mut ti_table = Array2::zeros(new_shape);
        for i in 0..new_shape.0 {
            for j in 0..new_shape.1 {
                if freq_acc[[i, j]] > 0.0 {
                    ti_table[[i, j]] = ti_sum[[i, j]] / freq_acc[[i, j]];
                }
            }
        }

        let value_table = if let (Some(sum), Some(_)) = (value_sum, &self.value_table) {
            let mut vt = Array2::zeros(new_shape);
            for i in 0..new_shape.0 {
                for j in 0..new_shape.1 {
                    if freq_acc[[i, j]] > 0.0 {
                        vt[[i, j]] = sum[[i, j]] / freq_acc[[i, j]];
                    }
                }
            }
            Some(vt)
        } else {
            None
        };

        WindRose::new(
            Array1::from_vec(new_wind_directions),
            Array1::from_vec(new_wind_speeds),
            ti_table,
            Some(freq_table),
            value_table,
            self.compute_zero_freq_occurrence,
            self.heterogeneous_map.clone(),
            self.multidim_conditions.clone(),
        )
        .unwrap_or_else(|_| self.clone())
    }

    pub fn downsample_mut(
        &mut self,
        wd_step: Option<Float>,
        ws_step: Option<Float>,
        ti_step: Option<Float>,
    ) {
        let downsampled = self.downsample(wd_step, ws_step, ti_step);
        *self = downsampled;
    }
    /// Upsample (resample) wind rose to more bins
    pub fn upsample(&self, wd_step: Float, ws_step: Float, method: &InterpMethod) -> Self {
        if self.wind_directions.is_empty() || self.wind_speeds.is_empty() {
            return self.clone();
        }

        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();

        let ws_step_current = if n_ws >= 2 {
            self.wind_speeds[1] - self.wind_speeds[0]
        } else {
            1.0
        };

        let wd_step_current = if n_dir >= 2 {
            self.wind_directions[1] - self.wind_directions[0]
        } else {
            360.0
        };

        if wd_step > wd_step_current {
            panic!(
                "wd_step ({}) is larger than current step ({}). Use downsample instead.",
                wd_step, wd_step_current
            );
        }

        if ws_step > ws_step_current {
            panic!(
                "ws_step ({}) is larger than current step ({}). Use downsample instead.",
                ws_step, ws_step_current
            );
        }

        let wd_min = self
            .wind_directions
            .iter()
            .fold(Float::INFINITY, |acc, &v| acc.min(v));
        let wd_max = self
            .wind_directions
            .iter()
            .fold(Float::NEG_INFINITY, |acc, &v| acc.max(v));
        let ws_min = self
            .wind_speeds
            .iter()
            .fold(Float::INFINITY, |acc, &v| acc.min(v));
        let ws_max = self
            .wind_speeds
            .iter()
            .fold(Float::NEG_INFINITY, |acc, &v| acc.max(v));

        let wd_range_min = wd_min - wd_step_current / 2.0;
        let wd_range_max = wd_max + wd_step_current / 2.0;
        let ws_range_min = ws_min - ws_step_current / 2.0;
        let ws_range_max = ws_max + ws_step_current / 2.0;

        let mut new_wind_directions = Vec::new();
        let mut wd_cursor = wd_range_min + wd_step / 2.0;
        while wd_cursor <= wd_range_max + 1e-9 {
            let wrapped = ((wd_cursor % 360.0) + 360.0) % 360.0;
            new_wind_directions.push(wrapped);
            wd_cursor += wd_step;
        }

        let mut new_wind_speeds = Vec::new();
        let mut ws_cursor = (ws_range_min + ws_step / 2.0).max(0.0);
        while ws_cursor <= ws_range_max + 1e-9 {
            new_wind_speeds.push(ws_cursor);
            ws_cursor += ws_step;
        }

        let interp = |table: &Array2, wd: Float, ws: Float| -> Float {
            match method {
                InterpMethod::Nearest => {
                    let wd_idx = Self::nearest_index(&self.wind_directions, wd);
                    let ws_idx = Self::nearest_index(&self.wind_speeds, ws);
                    table[[wd_idx, ws_idx]]
                }
                _ => {
                    let (wd0, wd1, wt) = Self::linear_bounds(&self.wind_directions, wd);
                    let (ws0, ws1, wt_ws) = Self::linear_bounds(&self.wind_speeds, ws);

                    let v00 = table[[wd0, ws0]];
                    let v10 = table[[wd1, ws0]];
                    let v01 = table[[wd0, ws1]];
                    let v11 = table[[wd1, ws1]];

                    let v0 = v00 * (1.0 - wt) + v10 * wt;
                    let v1 = v01 * (1.0 - wt) + v11 * wt;
                    v0 * (1.0 - wt_ws) + v1 * wt_ws
                }
            }
        };

        let new_shape = (new_wind_directions.len(), new_wind_speeds.len());

        let freq_table = Array2::from_shape_fn(new_shape, |(i, j)| {
            let wd = new_wind_directions[i];
            let ws = new_wind_speeds[j];
            interp(&self.freq_table, wd, ws)
        });

        let ti_table = Array2::from_shape_fn(new_shape, |(i, j)| {
            let wd = new_wind_directions[i];
            let ws = new_wind_speeds[j];
            interp(&self.ti_table, wd, ws)
        });

        let value_table = self.value_table.as_ref().map(|values| {
            Array2::from_shape_fn(new_shape, |(i, j)| {
                let wd = new_wind_directions[i];
                let ws = new_wind_speeds[j];
                interp(values, wd, ws)
            })
        });

        let freq_total: Float = freq_table.iter().sum();
        let mut freq_table_norm = freq_table.clone();
        if freq_total > 0.0 {
            for f in freq_table_norm.iter_mut() {
                *f /= freq_total;
            }
        }

        WindRose::new(
            Array1::from_vec(new_wind_directions),
            Array1::from_vec(new_wind_speeds),
            ti_table,
            Some(freq_table_norm),
            value_table,
            self.compute_zero_freq_occurrence,
            self.heterogeneous_map.clone(),
            self.multidim_conditions.clone(),
        )
        .unwrap_or_else(|_| self.clone())
    }
    pub fn upsample_mut(&mut self, wd_step: Float, ws_step: Float, method: &InterpMethod) {
        let upsampled = self.upsample(wd_step, ws_step, method);
        *self = upsampled;
    }

    /// Aggregate (downsample) wind rose to fewer bins - wrapper for downsample
    ///
    /// This is a convenience method that wraps downsample() for backwards compatibility
    /// with the Python FLORIS API.
    pub fn aggregate(&self, wd_step: Float, ws_step: Float, _inplace: bool) -> Self {
        self.downsample(Some(wd_step), Some(ws_step), None)
    }

    /// Resample (upsample) wind rose to more bins using interpolation - wrapper for upsample
    ///
    /// This is a convenience method that wraps upsample() for backwards compatibility
    /// with the Python FLORIS API.
    pub fn resample_by_interpolation(
        &self,
        wd_step: Float,
        ws_step: Float,
        method: &InterpMethod,
        _inplace: bool,
    ) -> Self {
        self.upsample(wd_step, ws_step, method)
    }

    fn nearest_index(arr: &Array1, target: Float) -> usize {
        let mut best_idx = 0usize;
        let mut best_diff = Float::INFINITY;

        for (idx, &val) in arr.iter().enumerate() {
            let diff = (val - target).abs();
            if diff < best_diff {
                best_diff = diff;
                best_idx = idx;
            }
        }

        best_idx
    }

    fn linear_bounds(arr: &Array1, target: Float) -> (usize, usize, Float) {
        if arr.len() == 1 {
            return (0, 0, 0.0);
        }

        let mut left = 0usize;
        while left + 1 < arr.len() && arr[left + 1] <= target {
            left += 1;
        }

        let right = (left + 1).min(arr.len() - 1);
        let span = arr[right] - arr[left];
        let weight = if span.abs() < 1e-9 {
            0.0
        } else {
            ((target - arr[left]) / span).clamp(0.0, 1.0)
        };

        (left, right, weight)
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

        let freq = &self.freq_table;

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

        let freq = &self.freq_table;

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
        self.wind_speeds.len() * self.wind_directions.len()
    }

    fn frequencies(&self) -> Array2 {
        self.freq_table.clone()
    }

    fn heterogeneous_inflow_config(&self) -> HeterogeneousInflowConfig {
        let n_conditions = self.n_conditions();
        let n_points = if let Some(ref het_map) = self.heterogeneous_map {
            het_map.x.len()
        } else {
            0
        };

        if n_points > 0 {
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
        let freq_table_unpack = self.get_freq_flat();
        let wind_directions_unpack =
            filter_by_nonzero_freq(&self.get_wd_flat(), &freq_table_unpack);
        let wind_speeds_unpack = filter_by_nonzero_freq(&self.get_ws_flat(), &freq_table_unpack);
        let ti_table_unpack = filter_by_nonzero_freq(&self.get_ti_flat(), &freq_table_unpack);
        let freq_table_unpack_filtered =
            filter_by_nonzero_freq(&self.get_freq_flat(), &freq_table_unpack);

        // Create 2D arrays from the filtered 1D arrays
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let freq_2d =
            Array2::from_shape_vec((n_dir, n_ws), freq_table_unpack_filtered.to_vec()).unwrap();

        let value_2d = if let Some(value_flat) = self.get_value_flat() {
            let value_table_unpack = filter_by_nonzero_freq(&value_flat, &freq_table_unpack);
            Array2::from_shape_vec((n_dir, n_ws), value_table_unpack.to_vec()).unwrap()
        } else {
            Array2::from_elem((n_dir, n_ws), 1.0)
        };

        (
            wind_directions_unpack,
            wind_speeds_unpack,
            ti_table_unpack,
            freq_2d,
            value_2d,
            self.heterogeneous_inflow_config(),
        )
    }
}

fn filter_by_nonzero_freq(values: &Array1, freq: &Array1) -> Array1 {
    values
        .iter()
        .zip(freq.iter())
        .filter_map(|(&v, &f)| if f > 0.0 { Some(v) } else { None })
        .collect()
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

        let wr = WindRose::new(wd, ws, ti_table, Some(freq), None, false, None, None).unwrap();
        assert_eq!(wr.n_conditions(), 12);
    }

    #[test]
    fn test_wind_rose_assign_ti() {
        let wd = Array1::from_vec(vec![0.0, 90.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti_table = Array2::from_elem((2, 2), 0.06);

        let mut wr = WindRose::new(wd, ws, ti_table, None, None, false, None, None).unwrap();
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

        let wr = WindRose::new(wd, ws, ti_table, Some(freq), None, false, None, None).unwrap();
        let ts = wr.to_time_series();

        // Should have 4 conditions (2 dirs × 2 speeds)
        assert_eq!(ts.n_conditions(), 4);
    }
}
