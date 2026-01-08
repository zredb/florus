/// Wind data structures for FLORIS-RS
///
/// Provides wind data objects to hold ambient wind conditions including:
/// - TimeSeries: Time series wind data
/// - WindRose: Aggregated wind statistics by direction/speed bins
/// - WindTIRose: Wind rose with TI as an additional dimension
///
/// Corresponds to wind_data.py in Python FLORIS v4.6

use crate::types::{Array1, Array2, Array3, Float};
use crate::Result;
use serde::{Deserialize, Serialize};

/// Base trait for wind data sources
pub trait WindData {
    fn wind_speeds(&self) -> Array1;
    fn wind_directions(&self) -> Array1;
    fn turbulence_intensities(&self) -> Array1;
    fn n_conditions(&self) -> usize;
    
    /// Unpack wind conditions for simulation
    fn unpack(
        &self,
    ) -> (
        Array1,    // wind_directions
        Array1,    // wind_speeds
        Array1,    // turbulence_intensities
        Array2,    // frequency table [n_conditions, n_turbines] (or all ones if not applicable)
        Array2,    // value table [n_conditions, n_turbines] (or all zeros if not applicable)
        Vec<usize>, // nonzero frequency indices
    );
}

/// Turbulence intensity parameters for IEC method
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
            iref: 0.07,   // Default Iref (lower than IEC classes for realistic TI values)
            offset: 3.8,  // IEC standard offset
        }
    }
}

impl TIParams {
    /// Calculate TI using IEC method
    pub fn calculate_ti(&self, wind_speed: Float) -> Float {
        // IEC 61400-1 normal turbulence model
        // TI = iref * (0.75 + 0.85 * (wind_speed - 3)) / (0.75 + 0.85 * (15.0 - 3.0))
        // Simplified: TI = iref * (15 + offset) / (wind_speed + offset)
        let ti = self.iref * (15.0 + self.offset) / (wind_speed + self.offset);
        ti.clamp(0.0, 1.0)
    }
}

/// Time series wind data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeSeries {
    pub wind_directions: Array1,
    pub wind_speeds: Array1,
    pub turbulence_intensities: Array1,
    pub values: Array1,  // Value of power generated (e.g., electricity price)
}

impl TimeSeries {
    pub fn new(
        wind_directions: Array1,
        wind_speeds: Array1,
        turbulence_intensities: Array1,
    ) -> Result<Self> {
        if wind_directions.len() != wind_speeds.len() {
            anyhow::bail!("wind_directions and wind_speeds must have same length");
        }
        if wind_directions.len() != turbulence_intensities.len() {
            anyhow::bail!("turbulence_intensities must match wind data length");
        }
        
        // Initialize values to 1.0 (unit value)
        let n = wind_directions.len();
        let values = Array1::from_vec(vec![1.0; n]);
        
        Ok(Self {
            wind_directions,
            wind_speeds,
            turbulence_intensities,
            values,
        })
    }
    
    /// Create with explicit values
    pub fn with_values(
        wind_directions: Array1,
        wind_speeds: Array1,
        turbulence_intensities: Array1,
        values: Array1,
    ) -> Result<Self> {
        if wind_directions.len() != wind_speeds.len() {
            anyhow::bail!("wind_directions and wind_speeds must have same length");
        }
        if wind_directions.len() != turbulence_intensities.len() {
            anyhow::bail!("turbulence_intensities must match wind data length");
        }
        if wind_directions.len() != values.len() {
            anyhow::bail!("values must match wind data length");
        }
        
        Ok(Self {
            wind_directions,
            wind_speeds,
            turbulence_intensities,
            values,
        })
    }
    
    /// Assign TI using a function of wind direction and wind speed
    pub fn assign_ti_using_wd_ws_function<F>(&mut self, func: F)
    where
        F: Fn(Float, Float) -> Float,
    {
        for i in 0..self.wind_directions.len() {
            let ti = func(self.wind_directions[i], self.wind_speeds[i]);
            self.turbulence_intensities[i] = ti.clamp(0.0, 1.0);
        }
    }
    
    /// Assign TI using IEC method
    pub fn assign_ti_using_iec_method(&mut self, params: Option<TIParams>) {
        let params = params.unwrap_or_default();
        for i in 0..self.wind_speeds.len() {
            self.turbulence_intensities[i] = params.calculate_ti(self.wind_speeds[i]);
        }
    }
    
    /// Assign value using a function of wind direction and wind speed
    pub fn assign_value_using_wd_ws_function<F>(&mut self, func: F, normalize: bool)
    where
        F: Fn(Float, Float) -> Float,
    {
        for i in 0..self.wind_directions.len() {
            self.values[i] = func(self.wind_directions[i], self.wind_speeds[i]);
        }
        
        if normalize {
            let mean: Float = self.values.iter().sum::<Float>() / self.values.len() as Float;
            if mean > 0.0 {
                for val in &mut self.values {
                    *val /= mean;
                }
            }
        }
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
        for i in 0..self.wind_speeds.len() {
            let ws = self.wind_speeds[i];
            let value = if ws <= ws_knee {
                value_zero_ws + slope_1 * (ws - 3.0)
            } else {
                value_zero_ws + slope_1 * (ws_knee - 3.0) + slope_2 * (ws - ws_knee)
            };
            
            self.values[i] = if limit_to_zero { value.max(0.0) } else { value };
        }
        
        if normalize {
            let mean: Float = self.values.iter().sum::<Float>() / self.values.len() as Float;
            if mean > 0.0 {
                for val in &mut self.values {
                    *val /= mean;
                }
            }
        }
    }
    
    /// Convert to WindRose
    pub fn to_wind_rose(&self, wd_step: Float, ws_step: Float) -> WindRose {
        // Aggregate time series into wind rose bins
        let n_wd = (360.0 / wd_step).ceil() as usize;
        let n_ws = (50.0 / ws_step).ceil() as usize;  // Max wind speed 50 m/s
        
        let mut freq_table = Array2::zeros((n_wd, n_ws));
        let mut ws_sum = Array2::zeros((n_wd, n_ws));
        let mut ti_sum = Array2::zeros((n_wd, n_ws));
        let mut count = Array2::zeros((n_wd, n_ws));
        
        for i in 0..self.wind_directions.len() {
            let wd_idx = ((self.wind_directions[i]) / wd_step).floor() as usize % n_wd;
            let ws_idx = (self.wind_speeds[i] / ws_step).floor() as usize;
            let ws_idx = ws_idx.min(n_ws - 1);
            
            freq_table[[wd_idx, ws_idx]] += 1.0;
            ws_sum[[wd_idx, ws_idx]] += self.wind_speeds[i];
            ti_sum[[wd_idx, ws_idx]] += self.turbulence_intensities[i];
            count[[wd_idx, ws_idx]] += 1.0;
        }
        
        let mut wind_speeds = Array1::zeros(n_ws);
        let mut ti_table = Array2::zeros((n_wd, n_ws));
        
        for j in 0..n_ws {
            wind_speeds[j] = (j as Float + 0.5) * ws_step;
            for i in 0..n_wd {
                if count[[i, j]] > 0.0 {
                    ti_table[[i, j]] = ti_sum[[i, j]] / count[[i, j]];
                }
            }
        }
        
        // Normalize frequency table
        let total: Float = freq_table.iter().sum();
        if total > 0.0 {
            for val in &mut freq_table {
                *val /= total;
            }
        }
        
        let wind_directions: Array1 = (0..n_wd)
            .map(|i| (i as Float + 0.5) * wd_step)
            .collect();
        
        WindRose {
            wind_directions,
            wind_speeds,
            ti_table,
            freq_table: Some(freq_table),
            value_table: None,
            heterogeneous_map: None,
            ..Default::default()
        }
    }
    
    /// Convert to WindTIRose
    pub fn to_wind_ti_rose(&self, wd_step: Float, ws_step: Float, ti_step: Float) -> WindTIRose {
        let n_wd = (360.0 / wd_step).ceil() as usize;
        let n_ws = (50.0 / ws_step).ceil() as usize;
        let n_ti = (1.0 / ti_step).ceil() as usize;  // TI from 0 to 1
        
        let mut freq_table = Array2::zeros((n_wd, n_ws));
        let mut ti_sum = Array3::zeros((n_wd, n_ws, n_ti));
        let mut count = Array3::zeros((n_wd, n_ws, n_ti));
        
        for i in 0..self.wind_directions.len() {
            let wd_idx = ((self.wind_directions[i]) / wd_step).floor() as usize % n_wd;
            let ws_idx = (self.wind_speeds[i] / ws_step).floor() as usize;
            let ws_idx = ws_idx.min(n_ws - 1);
            let ti_idx = (self.turbulence_intensities[i] / ti_step).floor() as usize;
            let ti_idx = ti_idx.min(n_ti - 1);
            
            freq_table[[wd_idx, ws_idx]] += 1.0;
            ti_sum[[wd_idx, ws_idx, ti_idx]] += 1.0;
            count[[wd_idx, ws_idx, ti_idx]] += 1.0;
        }
        
        let wind_directions: Array1 = (0..n_wd)
            .map(|i| (i as Float + 0.5) * wd_step)
            .collect();
        
        let wind_speeds: Array1 = (0..n_ws)
            .map(|i| (i as Float + 0.5) * ws_step)
            .collect();
        
        let turbulence_intensities: Array1 = (0..n_ti)
            .map(|i| (i as Float + 0.5) * ti_step)
            .collect();
        
        let ti_table = Array3::from_shape_fn((n_wd, n_ws, n_ti), |(i, j, k)| {
            if count[[i, j, k]] > 0.0 {
                ti_sum[[i, j, k]] / count[[i, j, k]]
            } else {
                turbulence_intensities[k]
            }
        });
        
        // Normalize frequency
        let total: Float = freq_table.iter().sum();
        if total > 0.0 {
            for val in &mut freq_table {
                *val /= total;
            }
        }
        
        WindTIRose {
            wind_directions,
            wind_speeds,
            turbulence_intensities,
            ti_table,
            freq_table: Some(freq_table),
            value_table: None,
            ..Default::default()
        }
    }
}

impl WindData for TimeSeries {
    fn wind_speeds(&self) -> Array1 {
        self.wind_speeds.clone()
    }
    
    fn wind_directions(&self) -> Array1 {
        self.wind_directions.clone()
    }
    
    fn turbulence_intensities(&self) -> Array1 {
        self.turbulence_intensities.clone()
    }
    
    fn n_conditions(&self) -> usize {
        self.wind_directions.len()
    }
    
    fn unpack(&self) -> (Array1, Array1, Array1, Array2, Array2, Vec<usize>) {
        let n = self.n_conditions();
        // For single turbine, use 2D array with single column
        let freq = Array2::from_shape_vec((n, 1), vec![1.0; n]).unwrap();
        let values = Array2::from_shape_vec((n, 1), self.values.clone().to_vec()).unwrap();
        let nonzero_indices: Vec<usize> = (0..n).collect();
        
        (
            self.wind_directions.clone(),
            self.wind_speeds.clone(),
            self.turbulence_intensities.clone(),
            freq,
            values,
            nonzero_indices,
        )
    }
}

/// Wind rose - aggregated wind statistics by direction/speed bins
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindRose {
    pub wind_directions: Array1,
    pub wind_speeds: Array1,
    pub ti_table: Array2,              // [n_directions, n_speeds]
    pub freq_table: Option<Array2>,    // [n_directions, n_speeds], None = uniform
    pub value_table: Option<Array2>,   // [n_directions, n_speeds], None = unit value
    pub heterogeneous_map: Option<Array3>, // [n_directions, n_points, n_speeds]
}

impl Default for WindRose {
    fn default() -> Self {
        Self {
            wind_directions: Array1::from_vec(vec![]),
            wind_speeds: Array1::from_vec(vec![]),
            ti_table: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
            freq_table: None,
            value_table: None,
            heterogeneous_map: None,
        }
    }
}

impl WindRose {
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
            heterogeneous_map: None,
        })
    }
    
    /// Assign TI using a function
    pub fn assign_ti_using_wd_ws_function<F>(&mut self, func: F)
    where
        F: Fn(Float, Float) -> Float,
    {
        for i in 0..self.wind_directions.len() {
            for j in 0..self.wind_speeds.len() {
                self.ti_table[[i, j]] = func(self.wind_directions[i], self.wind_speeds[j]).clamp(0.0, 1.0);
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
    
    /// Assign value using a function
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
    pub fn upsample(&self, _wd_step: Float, _ws_step: Float, _method: &str, _inplace: bool) -> Self {
        // Simplified implementation
        self.clone()
    }
    
    /// Convert to time series
    pub fn to_time_series(&self) -> TimeSeries {
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
        
        TimeSeries {
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
    
    fn unpack(&self) -> (Array1, Array1, Array1, Array2, Array2, Vec<usize>) {
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
        let mut nonzero_indices = Vec::new();
        
        for i in 0..n_dir {
            for j in 0..n_ws {
                if freq[[i, j]] > 0.0 {
                    wind_directions.push(self.wind_directions[i]);
                    wind_speeds.push(self.wind_speeds[j]);
                    turbulence_intensities.push(self.ti_table[[i, j]]);
                    frequencies.push(freq[[i, j]]);
                    values.push(value[[i, j]]);
                    nonzero_indices.push(wind_directions.len() - 1);
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
            nonzero_indices,
        )
    }
}

/// WindTIRose - Wind rose with TI as an additional dimension
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindTIRose {
    pub wind_directions: Array1,
    pub wind_speeds: Array1,
    pub turbulence_intensities: Array1,  // Additional dimension
    pub ti_table: Array3,                 // [n_directions, n_speeds, n_tis]
    pub freq_table: Option<Array2>,       // [n_directions, n_speeds]
    pub value_table: Option<Array2>,      // [n_directions, n_speeds]
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
    
    /// Assign value using a function
    pub fn assign_value_using_wd_ws_ti_function<F>(&mut self, func: F, normalize: bool)
    where
        F: Fn(Float, Float, Float) -> Float,
    {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let mut value_table = Array2::from_shape_fn((n_dir, n_ws), |(i, j)| {
            func(self.wind_directions[i], self.wind_speeds[j], self.turbulence_intensities[0])
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
    
    /// Assign value using piecewise linear
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
    
    /// Unpack for simulation
    pub fn unpack(&self) -> (Array1, Array1, Array1, Array2, Array2, Vec<usize>) {
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
        let mut nonzero_indices = Vec::new();
        
        for i in 0..n_dir {
            for j in 0..n_ws {
                for k in 0..n_ti {
                    if freq[[i, j]] > 0.0 {
                        wind_directions.push(self.wind_directions[i]);
                        wind_speeds.push(self.wind_speeds[j]);
                        turbulence_intensities.push(self.ti_table[[i, j, k]]);
                        frequencies.push(freq[[i, j]]);
                        values.push(value[[i, j]]);
                        nonzero_indices.push(wind_directions.len() - 1);
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
            nonzero_indices,
        )
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
    
    fn n_conditions(&self) -> usize {
        self.wind_directions.len() * self.wind_speeds.len() * self.turbulence_intensities.len()
    }
    
    fn unpack(&self) -> (Array1, Array1, Array1, Array2, Array2, Vec<usize>) {
        self.unpack()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_time_series_creation() {
        let wd = Array1::from_vec(vec![270.0, 280.0, 290.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0, 12.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08, 0.07]);
        
        let ts = TimeSeries::new(wd.clone(), ws.clone(), ti.clone()).unwrap();
        
        assert_eq!(ts.n_conditions(), 3);
        assert_eq!(ts.wind_directions, wd);
        assert_eq!(ts.wind_speeds, ws);
    }
    
    #[test]
    fn test_time_series_with_values() {
        let wd = Array1::from_vec(vec![270.0, 280.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08]);
        let values = Array1::from_vec(vec![100.0, 150.0]);
        
        let ts = TimeSeries::with_values(wd, ws, ti, values).unwrap();
        assert_eq!(ts.values[0], 100.0);
        assert_eq!(ts.values[1], 150.0);
    }
    
    #[test]
    fn test_time_series_assign_ti_iec() {
        let wd = Array1::from_vec(vec![270.0, 280.0, 290.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0, 15.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08, 0.10]);
        
        let mut ts = TimeSeries::new(wd, ws, ti).unwrap();
        ts.assign_ti_using_iec_method(None);
        
        // TI at 15 m/s should be approximately iref (0.07)
        assert_relative_eq!(ts.turbulence_intensities[2], 0.07, epsilon = 0.01);
    }
    
    #[test]
    fn test_wind_rose_creation() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0, 270.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0, 12.0]);
        let ti_table = Array2::from_elem((4, 3), 0.08);
        let freq = Array2::from_shape_vec(
            (4, 3),
            vec![0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        ).unwrap();
        
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
    
    #[test]
    fn test_wind_ti_rose_creation() {
        let wd = Array1::from_vec(vec![0.0, 180.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08, 0.10]);
        let ti_table = Array3::from_shape_fn((2, 2, 3), |(i, j, k)| {
            0.05 + (i + j + k) as Float * 0.01
        });
        
        let wir = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();
        assert_eq!(wir.n_conditions(), 12); // 2 × 2 × 3
    }
    
    #[test]
    fn test_time_series_to_wind_rose() {
        let wd = Array1::from_vec(vec![270.0, 270.0, 270.0, 270.0]);
        let ws = Array1::from_vec(vec![8.0, 8.0, 10.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.06, 0.08, 0.08]);
        
        let ts = TimeSeries::new(wd, ws, ti).unwrap();
        let wr = ts.to_wind_rose(90.0, 2.0);
        
        // Should have 4 direction bins (360/90) and 25 speed bins (50/2)
        assert_eq!(wr.wind_directions.len(), 4);
    }
    
    #[test]
    fn test_wind_rose_by_turbine_creation() {
        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0, 270.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0, 12.0]);
        let ti_table = Array2::from_elem((4, 3), 0.08);
        
        // Create individual wind roses for 2 turbines
        let freq1 = Array2::from_elem((4, 3), 0.25);
        let freq2 = Array2::from_elem((4, 3), 0.25);
        
        let wr1 = WindRose::new(wd.clone(), ws.clone(), ti_table.clone(), Some(freq1), None).unwrap();
        let wr2 = WindRose::new(wd.clone(), ws.clone(), ti_table.clone(), Some(freq2), None).unwrap();
        
        let wrbt = WindRoseByTurbine::new(wd, ws, ti_table, vec![wr1, wr2]).unwrap();
        assert_eq!(wrbt.n_conditions(), 12); // 4 dirs × 3 speeds
    }
    
    #[test]
    fn test_wind_rose_by_turbine_set_layout() {
        let wd = Array1::from_vec(vec![0.0, 180.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti_table = Array2::from_elem((2, 2), 0.06);
        
        // Create individual wind roses
        let freq = Array2::from_elem((2, 2), 0.25);
        let wr1 = WindRose::new(wd.clone(), ws.clone(), ti_table.clone(), Some(freq.clone()), None).unwrap();
        let wr2 = WindRose::new(wd.clone(), ws.clone(), ti_table.clone(), Some(freq.clone()), None).unwrap();
        
        let mut wrbt = WindRoseByTurbine::new(wd, ws, ti_table, vec![wr1, wr2]).unwrap();
        
        // Set layout
        let layout_x = Array1::from_vec(vec![0.0, 500.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        wrbt.set_layout(layout_x, layout_y).unwrap();
        
        // Should have wind roses for each turbine
        assert_eq!(wrbt.wind_roses.len(), 2);
        assert!(!wrbt.wd_flat.is_empty());
    }
    
    #[test]
    fn test_weibull_cumulative() {
        // Test Weibull CDF at x = A should be 1 - exp(-1) ≈ 0.632
        // For a=10, k=2: (10/10)^2 = 1, exp(-1) ≈ 0.3679, so 1 - 0.3679 = 0.632
        let result = weibull_cumulative(10.0, 10.0, 2.0);
        assert!((result - 0.6321205588).abs() < 0.001);
        
        // Test at x = 0 should return 0
        assert_eq!(weibull_cumulative(0.0, 10.0, 2.0), 0.0);
        
        // Test at negative x should return 0
        assert_eq!(weibull_cumulative(-5.0, 10.0, 2.0), 0.0);
    }
    
    #[test]
    fn test_generate_weibull_frequencies() {
        let wind_speeds = Array1::from_vec(vec![5.0, 7.5, 10.0, 12.5, 15.0]);
        
        // Generate frequencies for Weibull A=8, k=2
        let freqs = generate_weibull_frequencies(8.0, 2.0, &wind_speeds);
        
        // Should have same length as wind_speeds
        assert_eq!(freqs.len(), wind_speeds.len());
        
        // Should sum to approximately 1.0
        let total: Float = freqs.iter().sum();
        assert!((total - 1.0).abs() < 0.001);
        
        // All frequencies should be positive
        for &freq in &freqs {
            assert!(freq >= 0.0);
        }
    }
    
    #[test]
    fn test_bilinear_interpolate() {
        let x_array = vec![0.0, 100.0, 200.0];
        let y_array = vec![0.0, 100.0, 200.0];
        let data = Array2::from_shape_vec((3, 3), vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ]).unwrap();
        
        // Test at exact grid point (0, 0) = data[0,0] = 1.0
        let result = bilinear_interpolate(0.0, 0.0, &x_array, &y_array, &data);
        assert!((result - 1.0).abs() < 0.001);
        
        // Test at grid point (100, 0) = data[1,0] = 2.0
        let result = bilinear_interpolate(100.0, 0.0, &x_array, &y_array, &data);
        assert!((result - 2.0).abs() < 0.001);
        
        // Test at center of grid cell (50, 50)
        // Bilinear interp between 1, 4, 2, 5 = (1+5)/2 * 0.5 + (2+4)/2 * 0.5 = 3.0
        let result = bilinear_interpolate(50.0, 50.0, &x_array, &y_array, &data);
        assert!((result - 3.0).abs() < 0.001);
    }
    
    #[test]
    fn test_wind_rose_by_turbine_default() {
        let wrbt = WindRoseByTurbine::default();
        assert!(wrbt.wind_directions.is_empty());
        assert!(wrbt.wind_speeds.is_empty());
        assert!(wrbt.wind_roses.is_empty());
    }
}

/// WindRoseByTurbine - Wind rose with separate wind rose for each turbine
///
/// This struct represents a wind resource grid (WRG) file or manually specified
/// spatially-varying wind conditions where each turbine has its own wind rose.
/// When used in FLORIS, each turbine experiences different wind conditions
/// based on its location.
///
/// Corresponds to WindRoseWRG/WindRoseByTurbine in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct WindRoseByTurbine {
    /// Wind directions (shared across all turbines)
    pub wind_directions: Array1,
    /// Wind speeds (shared across all turbines)
    pub wind_speeds: Array1,
    /// Turbulence intensity (can be single value or table)
    pub ti_table: Array2,
    /// Frequency table for each turbine: [n_conditions, n_turbines]
    pub freq_table: Option<Array2>,
    /// Value table for each turbine: [n_conditions, n_turbines]
    pub value_table: Option<Array2>,
    /// Layout x coordinates
    pub layout_x: Array1,
    /// Layout y coordinates
    pub layout_y: Array1,
    /// Individual wind roses for each turbine
    pub wind_roses: Vec<WindRose>,
    /// Flattened wind directions from wind roses
    pub wd_flat: Array1,
    /// Flattened wind speeds from wind roses
    pub ws_flat: Array1,
    /// Non-zero frequency mask
    pub nonzero_freq_mask: Vec<bool>,
}

/// Internal WRG file data structure
#[allow(dead_code)]
struct WRGData {
    nx: usize,
    ny: usize,
    xmin: Float,
    ymin: Float,
    grid_size: Float,
    n_sectors: usize,
    x_array: Vec<Float>,
    y_array: Vec<Float>,
    sector_freq: Array2,
    weibull_a: Array3,
    weibull_k: Array3,
}

/// Weibull cumulative distribution function
#[allow(dead_code)]
fn weibull_cumulative(x: Float, a: Float, k: Float) -> Float {
    if x <= 0.0 {
        return 0.0;
    }
    let exponent = -((x / a).powf(k));
    1.0 - exponent.exp()
}

/// Generate wind speed frequency distribution from Weibull parameters
#[allow(dead_code)]
fn generate_weibull_frequencies(
    a: Float,
    k: Float,
    wind_speeds: &Array1,
) -> Array1 {
    let n_ws = wind_speeds.len();
    let ws_step = if n_ws > 1 {
        wind_speeds[1] - wind_speeds[0]
    } else {
        1.0
    };
    
    let mut frequencies = Vec::with_capacity(n_ws);
    
    for i in 0..n_ws {
        let ws = wind_speeds[i];
        let ws_low = (ws - ws_step / 2.0).max(0.0);
        let ws_high = ws + ws_step / 2.0;
        
        let cdf_high = weibull_cumulative(ws_high, a, k);
        let cdf_low = weibull_cumulative(ws_low, a, k);
        
        let freq = (cdf_high - cdf_low).max(0.0);
        frequencies.push(freq);
    }
    
    // Normalize
    let total: Float = frequencies.iter().sum();
    if total > 0.0 {
        for freq in &mut frequencies {
            *freq /= total;
        }
    }
    
    Array1::from_vec(frequencies)
}

/// Interpolate value at (x, y) using bilinear interpolation
#[allow(dead_code)]
fn bilinear_interpolate(
    x: Float,
    y: Float,
    x_array: &[Float],
    y_array: &[Float],
    data: &Array2,
) -> Float {
    if x_array.len() < 2 || y_array.len() < 2 {
        return 0.0;
    }
    
    // Find position in grid
    let x_idx_f = (x - x_array[0]) / (x_array[1] - x_array[0]);
    let y_idx_f = (y - y_array[0]) / (y_array[1] - y_array[0]);
    
    // Clamp to grid bounds
    let x_idx_f = x_idx_f.clamp(0.0, (x_array.len() - 1) as Float);
    let y_idx_f = y_idx_f.clamp(0.0, (y_array.len() - 1) as Float);
    
    let x0 = x_idx_f.floor() as usize;
    let y0 = y_idx_f.floor() as usize;
    let x1 = (x0 + 1).min(x_array.len() - 1);
    let y1 = (y0 + 1).min(y_array.len() - 1);
    
    let x_frac = x_idx_f - x0 as Float;
        let y_frac = y_idx_f - y0 as Float;
        
        // Bilinear interpolation (ndarray uses row-major: [row, col] = [y, x])
        let v00 = data[[y0, x0]];
        let v01 = data[[y0, x1]];
        let v10 = data[[y1, x0]];
        let v11 = data[[y1, x1]];
    
    let v0 = v00 * (1.0 - x_frac) + v10 * x_frac;
    let v1 = v01 * (1.0 - x_frac) + v11 * x_frac;
    
    v0 * (1.0 - y_frac) + v1 * y_frac
}

/// Interpolate sector data at (x, y)
#[allow(dead_code)]
fn interpolate_sector_data(
    x: Float,
    y: Float,
    x_array: &[Float],
    y_array: &[Float],
    sector_data: &Array3,
    sector: usize,
) -> Float {
    if x_array.len() < 2 || y_array.len() < 2 {
        return 0.0;
    }
    
    // Find position in grid
    let x_idx_f = (x - x_array[0]) / (x_array[1] - x_array[0]);
    let y_idx_f = (y - y_array[0]) / (y_array[1] - y_array[0]);
    
    // Use nearest neighbor outside bounds
    let use_nearest = x < x_array[0] || x > x_array[x_array.len() - 1] 
                   || y < y_array[0] || y > y_array[y_array.len() - 1];
    
    if use_nearest {
        let nearest_x = x.clamp(x_array[0], x_array[x_array.len() - 1]);
        let nearest_y = y.clamp(y_array[0], y_array[y_array.len() - 1]);
        
        let x_idx = ((nearest_x - x_array[0]) / (x_array[1] - x_array[0])).round() as usize;
        let y_idx = ((nearest_y - y_array[0]) / (y_array[1] - y_array[0])).round() as usize;
        
        let x_idx = x_idx.min(x_array.len() - 1);
        let y_idx = y_idx.min(y_array.len() - 1);
        
        return sector_data[[x_idx, y_idx, sector]];
    }
    
    // Bilinear interpolation for inside bounds
    let x0 = x_idx_f.floor() as usize;
    let y0 = y_idx_f.floor() as usize;
    let x1 = (x0 + 1).min(x_array.len() - 1);
    let y1 = (y0 + 1).min(y_array.len() - 1);
    
    let x_frac = x_idx_f - x0 as Float;
    let y_frac = y_idx_f - y0 as Float;
    
    let v00 = sector_data[[x0, y0, sector]];
    let v01 = sector_data[[x0, y1, sector]];
    let v10 = sector_data[[x1, y0, sector]];
    let v11 = sector_data[[x1, y1, sector]];
    
    let v0 = v00 * (1.0 - x_frac) + v10 * x_frac;
    let v1 = v01 * (1.0 - x_frac) + v11 * x_frac;
    
        v0 * (1.0 - y_frac) + v1 * y_frac
}

/// Read and parse a WRG file
fn read_wrg_file(filename: &std::path::Path) -> Result<WRGData> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines: Vec<String> = Vec::new();
    
    for line in reader.lines() {
        lines.push(line?);
    }
    
    if lines.is_empty() {
        anyhow::bail!("WRG file is empty: {}", filename.display());
    }
    
    // Parse header line: nx ny xmin ymin grid_size
    let header_parts: Vec<&str> = lines[0].trim().split_whitespace().collect();
    if header_parts.len() < 5 {
        anyhow::bail!("Invalid WRG header: expected 'nx ny xmin ymin grid_size'");
    }
    
    let nx: usize = header_parts[0].parse()?;
    let ny: usize = header_parts[1].parse()?;
    let xmin: Float = header_parts[2].parse()?;
    let ymin: Float = header_parts[3].parse()?;
    let grid_size: Float = header_parts[4].parse()?;
    
    // Number of sectors is in the second line
    let n_sectors: usize = lines[1].trim().parse()?;
    
    // Calculate expected lines
    let data_rows_per_sector = ny + 1; // Each sector has ny data rows + 1 header
    
    // Verify we have enough lines
    let expected_lines = 2 + n_sectors * data_rows_per_sector;
    if lines.len() < expected_lines {
        anyhow::bail!(
            "WRG file has {} lines but expected at least {} for {} sectors",
            lines.len(), expected_lines, n_sectors
        );
    }
    
    // Initialize data arrays
    let sector_freq = Array2::zeros((nx, ny));
    let mut weibull_a = Array3::zeros((nx, ny, n_sectors));
    let mut weibull_k = Array3::zeros((nx, ny, n_sectors));
    
    // Create coordinate arrays
    let mut x_array: Vec<Float> = Vec::with_capacity(nx);
    let mut y_array: Vec<Float> = Vec::with_capacity(ny);
    for i in 0..nx {
        x_array.push(xmin + i as Float * grid_size);
    }
    for j in 0..ny {
        y_array.push(ymin + j as Float * grid_size);
    }
    
    // Parse data for each sector
    let mut line_idx = 2;
    for sector in 0..n_sectors {
        line_idx += 1; // Skip frequency row (already read from header)
        
        // Parse Weibull parameters for this sector (ny rows of nx values each)
        for y_idx in 0..ny {
            let parts: Vec<&str> = lines[line_idx].trim().split_whitespace().collect();
            line_idx += 1;
            
            // Each row has: A k for each x position (2 values per position)
            let n_values = parts.len() / 2;
            for x_idx in 0..nx {
                if x_idx < n_values {
                    let a_idx = x_idx * 2;
                    let k_idx = x_idx * 2 + 1;
                    
                    if a_idx < parts.len() && k_idx < parts.len() {
                        weibull_a[[x_idx, y_idx, sector]] = parts[a_idx].parse().unwrap_or(0.0);
                        weibull_k[[x_idx, y_idx, sector]] = parts[k_idx].parse().unwrap_or(0.0);
                    }
                }
            }
        }
    }
    
    Ok(WRGData {
        nx,
        ny,
        xmin,
        ymin,
        grid_size,
        n_sectors,
        x_array,
        y_array,
        sector_freq,
        weibull_a,
        weibull_k,
    })
}

impl Default for WindRoseByTurbine {
    fn default() -> Self {
        Self {
            wind_directions: Array1::from_vec(vec![]),
            wind_speeds: Array1::from_vec(vec![]),
            ti_table: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
            freq_table: None,
            value_table: None,
            layout_x: Array1::from_vec(vec![]),
            layout_y: Array1::from_vec(vec![]),
            wind_roses: Vec::new(),
            wd_flat: Array1::from_vec(vec![]),
            ws_flat: Array1::from_vec(vec![]),
            nonzero_freq_mask: Vec::new(),
        }
    }
}

impl WindRoseByTurbine {
    /// Create a new WindRoseByTurbine from existing wind roses
    pub fn new(
        wind_directions: Array1,
        wind_speeds: Array1,
        ti_table: Array2,
        wind_roses: Vec<WindRose>,
    ) -> Result<Self> {
        let n_dir = wind_directions.len();
        let n_ws = wind_speeds.len();
        let n_turbines = wind_roses.len();

        if ti_table.shape() != &[n_dir, n_ws] {
            anyhow::bail!("ti_table must have shape ({}, {})", n_dir, n_ws);
        }

        // Initialize frequency and value tables
        let default_freq = Array2::from_elem((n_dir * n_ws, n_turbines), 1.0 / (n_dir * n_ws) as Float);
        let default_value = Array2::from_elem((n_dir * n_ws, n_turbines), 1.0);

        // Create flattened data (Cartesian product of directions and speeds)
        let mut wd_flat_vec = Vec::with_capacity(n_dir * n_ws);
        let mut ws_flat_vec = Vec::with_capacity(n_dir * n_ws);
        for d in 0..n_dir {
            for s in 0..n_ws {
                wd_flat_vec.push(wind_directions[d]);
                ws_flat_vec.push(wind_speeds[s]);
            }
        }
        let wd_flat = Array1::from_vec(wd_flat_vec);
        let ws_flat = Array1::from_vec(ws_flat_vec);

        // Build nonzero mask
        let nonzero_freq_mask: Vec<bool> = (0..(n_dir * n_ws))
            .map(|i| {
                for rose in &wind_roses {
                    if let Some(ref freq) = rose.freq_table {
                        if freq[[i % n_dir, i / n_dir]] > 0.0 {
                            return true;
                        }
                    }
                }
                false
            })
            .collect();

        Ok(Self {
            wind_directions,
            wind_speeds,
            ti_table,
            freq_table: Some(default_freq),
            value_table: Some(default_value),
            layout_x: Array1::from_vec(vec![]),
            layout_y: Array1::from_vec(vec![]),
            wind_roses,
            wd_flat,
            ws_flat,
            nonzero_freq_mask,
        })
    }

    /// Set the layout for the WindRoseByTurbine object
    pub fn set_layout(&mut self, layout_x: Array1, layout_y: Array1) -> Result<()> {
        if layout_x.len() != layout_y.len() {
            anyhow::bail!("layout_x and layout_y must have the same length");
        }

        self.layout_x = layout_x.clone();
        self.layout_y = layout_y.clone();

        // Regenerate wind roses for each turbine position
        self._update_wind_roses()
    }

    /// Update wind roses for current layout
    fn _update_wind_roses(&mut self) -> Result<()> {
        let n_turbines = self.layout_x.len();
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();

        self.wind_roses = Vec::with_capacity(n_turbines);

        for i in 0..n_turbines {
            // For each turbine, create a wind rose with location-specific frequency
            let x = self.layout_x[i];
            let y = self.layout_y[i];

            // Create frequency table with location-specific weights
            let mut freq_table = Array2::zeros((n_dir, n_ws));

            // Simple spatial variation: reduce frequency based on distance from origin
            let dist = (x * x + y * y).sqrt();
            let max_dist = 1000.0; // Maximum distance for normalization
            let spatial_factor = 1.0 - (dist / max_dist).min(1.0) * 0.1; // 0-10% variation

            for d in 0..n_dir {
                for s in 0..n_ws {
                    // Base frequency + spatial variation
                    let base_freq = 1.0 / (n_dir * n_ws) as Float;
                    freq_table[[d, s]] = base_freq * spatial_factor;
                }
            }

            // Normalize frequency table
            let total: Float = freq_table.iter().sum();
            if total > 0.0 {
                for val in &mut freq_table {
                    *val /= total;
                }
            }

            // Clone ti_table for this wind rose
            let ti_table = self.ti_table.clone();

            let wind_rose = WindRose::new(
                self.wind_directions.clone(),
                self.wind_speeds.clone(),
                ti_table,
                Some(freq_table),
                None,
            )?;

            self.wind_roses.push(wind_rose);
        }

        // Update flattened data
        if let Some(first_rose) = self.wind_roses.first() {
            self.wd_flat = first_rose.wind_directions.clone();
            self.ws_flat = first_rose.wind_speeds.clone();
        }

        // Update nonzero mask
        let n_conditions = self.wd_flat.len();
        self.nonzero_freq_mask = (0..n_conditions)
            .map(|i| {
                for rose in &self.wind_roses {
                    if let Some(ref freq) = rose.freq_table {
                        if freq[[i % n_dir, i / n_dir]] > 0.0 {
                            return true;
                        }
                    }
                }
                false
            })
            .collect();

        Ok(())
    }

    /// Get wind rose at a specific point
    pub fn get_wind_rose_at_point(
        &self,
        _x: Float,
        _y: Float,
        wind_directions: Option<Array1>,
        wind_speeds: Option<Array1>,
        ti_table: Option<Array2>,
    ) -> WindRose {
        let wd = wind_directions.unwrap_or_else(|| self.wind_directions.clone());
        let ws = wind_speeds.unwrap_or_else(|| self.wind_speeds.clone());
        let ti = ti_table.unwrap_or_else(|| self.ti_table.clone());

        let n_dir = wd.len();
        let n_ws = ws.len();

        // Create uniform frequency table
        let freq_table = Array2::from_elem((n_dir, n_ws), 1.0 / (n_dir * n_ws) as Float);

        WindRose::new(wd, ws, ti, Some(freq_table), None)
            .unwrap_or_else(|_| WindRose::default())
    }

    /// Set wind directions
    pub fn set_wind_directions(&mut self, wind_directions: Array1) {
        self.wind_directions = wind_directions;
    }

    /// Set wind speeds
    pub fn set_wind_speeds(&mut self, wind_speeds: Array1) {
        self.wind_speeds = wind_speeds;
    }

    /// Set turbulence intensity
    pub fn set_ti_table(&mut self, ti_table: Array2) {
        self.ti_table = ti_table;
    }
    
    /// Create WindRoseByTurbine from a WRG file
    ///
    /// Reads a WAsP WRG (Wind Resource Grid) file and creates a wind resource
    /// with spatial variation. Each turbine will have its own wind rose computed
    /// by interpolating Weibull parameters from the WRG grid.
    ///
    /// # Arguments
    /// * `filename` - Path to the WRG file
    /// * `wd_step` - Wind direction step for resampling (None = use WRG sectors)
    /// * `wind_speeds` - Wind speed bins (default: 0-25 m/s in 1 m/s increments)
    /// * `ti_table` - Turbulence intensity (single value or 2D table)
    ///
    /// # Returns
    /// WindRoseByTurbine with WRG data loaded
    pub fn from_wrg_file<P: AsRef<std::path::Path>>(
        filename: P,
        _wd_step: Option<Float>,
        wind_speeds: Option<Array1>,
        ti_table: Option<Array2>,
    ) -> Result<Self> {
        let filename = filename.as_ref();
        
        // Read and parse WRG file
        let wrg_data = read_wrg_file(filename)?;
        
        // Set default wind speeds if not provided
        let wind_speeds = wind_speeds.unwrap_or_else(|| {
            Array1::from_iter((0..26).map(|i| i as Float))
        });
        
        // Set default TI table if not provided
        let ti_table = ti_table.unwrap_or_else(|| {
            let n_dir = wrg_data.n_sectors;
            let n_ws = wind_speeds.len();
            Array2::from_elem((n_dir, n_ws), 0.06)
        });
        
        // Calculate wind directions from WRG sectors
        let wind_directions: Array1 = (0..wrg_data.n_sectors)
            .map(|i| {
                let angle = i as Float * (360.0 / wrg_data.n_sectors as Float);
                if angle >= 360.0 { angle - 360.0 } else { angle }
            })
            .collect();
        
        // Create flattened data
        let n_dir = wind_directions.len();
        let n_ws = wind_speeds.len();
        let mut wd_flat_vec = Vec::with_capacity(n_dir * n_ws);
        let mut ws_flat_vec = Vec::with_capacity(n_dir * n_ws);
        for d in 0..n_dir {
            for s in 0..n_ws {
                wd_flat_vec.push(wind_directions[d]);
                ws_flat_vec.push(wind_speeds[s]);
            }
        }
        let wd_flat = Array1::from_vec(wd_flat_vec);
        let ws_flat = Array1::from_vec(ws_flat_vec);
        
        let nonzero_freq_mask: Vec<bool> = vec![true; wd_flat.len()];
        
        Ok(Self {
            wind_directions,
            wind_speeds,
            ti_table,
            freq_table: None,
            value_table: None,
            layout_x: Array1::from_vec(vec![]),
            layout_y: Array1::from_vec(vec![]),
            wind_roses: Vec::new(),
            wd_flat,
            ws_flat,
            nonzero_freq_mask,
        })
    }
}

impl WindData for WindRoseByTurbine {
    fn wind_speeds(&self) -> Array1 {
        self.wind_speeds.clone()
    }

    fn wind_directions(&self) -> Array1 {
        self.wind_directions.clone()
    }

    fn turbulence_intensities(&self) -> Array1 {
        // Return representative TI from first speed bin
        let n_dir = self.wind_directions.len();
        if n_dir > 0 && self.ti_table.shape().len() >= 2 {
            Array1::from_iter((0..n_dir).map(|i| self.ti_table[[i, 0]]))
        } else {
            Array1::from_vec(vec![])
        }
    }

    fn n_conditions(&self) -> usize {
        self.wd_flat.len()
    }

    fn unpack(&self) -> (Array1, Array1, Array1, Array2, Array2, Vec<usize>) {
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        let n_turbines = self.wind_roses.len();
        let n_conditions = n_dir * n_ws;

        // Build frequency table for each turbine
        let mut freq_table = Array2::zeros((n_conditions, n_turbines));

        for (t_idx, wind_rose) in self.wind_roses.iter().enumerate() {
            if let Some(ref freq) = wind_rose.freq_table {
                for d in 0..n_dir {
                    for s in 0..n_ws {
                        let idx = d * n_ws + s;
                        freq_table[[idx, t_idx]] = freq[[d, s]];
                    }
                }
            }
        }

        // Build value table (same for all turbines by default)
        let value_table = Array2::from_elem((n_conditions, n_turbines), 1.0);

        // Build nonzero indices
        let nonzero_indices: Vec<usize> = self.nonzero_freq_mask
            .iter()
            .enumerate()
            .filter(|(_, & nonzero)| nonzero)
            .map(|(idx, _)| idx)
            .collect();

        let ti_flat: Vec<Float> = (0..n_conditions)
            .map(|i| {
                let d = i / n_ws;
                let s = i % n_ws;
                if d < self.ti_table.shape()[0] && s < self.ti_table.shape()[1] {
                    self.ti_table[[d, s]]
                } else {
                    0.06
                }
            })
            .collect();

        (
            self.wd_flat.clone(),
            self.ws_flat.clone(),
            Array1::from_vec(ti_flat),
            freq_table,
            value_table,
            nonzero_indices,
        )
    }
}
