//! Time series wind data.
//!
//! TimeSeries represents sequential wind measurements or synthetic wind
//! conditions over time, with wind direction, speed, and turbulence intensity.

use crate::heterogeneous_map::HeterogeneousInflowConfig;
use crate::types::{Array1, Array2, Array3, Float};
use crate::wind_data::traits::{TIParams, WindData};
use crate::wind_data::WindRose;
use crate::Result;
use serde::{Deserialize, Serialize};

/// Time series wind data
///
/// TimeSeries represents a sequence of wind conditions over time,
/// where each time step has a wind direction, wind speed, and turbulence
/// intensity. It can also include associated values (e.g., electricity prices).
///
/// Corresponds to TimeSeries in Python FLORIS v4.6
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeSeries {
    /// Wind directions for each time step [n_times]
    pub wind_directions: Array1,
    /// Wind speeds for each time step [n_times]
    pub wind_speeds: Array1,
    /// Turbulence intensities for each time step [n_times]
    pub turbulence_intensities: Array1,
    /// Value at each time step (e.g., electricity price) [n_times]
    pub values: Array1,
}

impl TimeSeries {
    /// Create a new TimeSeries
    ///
    /// # Arguments
    /// * `wind_directions` - Wind directions for each time step
    /// * `wind_speeds` - Wind speeds for each time step
    /// * `turbulence_intensities` - TI values for each time step
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
    ///
    /// Aggregates time series into wind rose bins.
    pub fn to_wind_rose(&self, wd_step: Float, ws_step: Float) -> WindRose {
        // Aggregate time series into wind rose bins
        let n_wd = (360.0 / wd_step).ceil() as usize;
        let n_ws = (50.0 / ws_step).ceil() as usize; // Max wind speed 50 m/s

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

        let wind_directions: Array1 = (0..n_wd).map(|i| (i as Float + 0.5) * wd_step).collect();

        WindRose {
            wind_directions,
            wind_speeds,
            ti_table,
            freq_table: Some(freq_table),
            value_table: None,
            heterogeneous_map: None,
            multidim_conditions: None,
            ..Default::default()
        }
    }

    /// Convert to WindTIRose
    ///
    /// Aggregates time series into wind rose bins with TI as a dimension.
    pub fn to_wind_ti_rose(&self, wd_step: Float, ws_step: Float, ti_step: Float) -> super::WindTIRose {
        let n_wd = (360.0 / wd_step).ceil() as usize;
        let n_ws = (50.0 / ws_step).ceil() as usize;
        let n_ti = (1.0 / ti_step).ceil() as usize; // TI from 0 to 1

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

        let wind_directions: Array1 = (0..n_wd).map(|i| (i as Float + 0.5) * wd_step).collect();

        let wind_speeds: Array1 = (0..n_ws).map(|i| (i as Float + 0.5) * ws_step).collect();

        let turbulence_intensities: Array1 =
            (0..n_ti).map(|i| (i as Float + 0.5) * ti_step).collect();

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

        super::WindTIRose {
            wind_directions,
            wind_speeds,
            turbulence_intensities,
            ti_table,
            freq_table: Some(freq_table),
            value_table: None,
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

    fn heterogeneous_inflow_config(&self) -> HeterogeneousInflowConfig {
        HeterogeneousInflowConfig {
            x: Array1::from_vec(vec![]),
            y: Array1::from_vec(vec![]),
            z: None,
            wind_speeds: Some(self.wind_speeds.clone()),
            wind_directions: Some(self.wind_directions.clone()),
            speed_multipliers: Array2::from_shape_vec((self.wind_directions.len(), 0), vec![])
                .unwrap(),
        }
    }

    fn set_layout(&mut self, _layout_x: &Option<Array1>, _layout_y: &Option<Array1>) {
        // TimeSeries doesn't support layout changes
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
        let n = self.n_conditions();
        // For single turbine, use 2D array with single column
        let freq = Array2::from_shape_vec((n, 1), vec![1.0; n]).unwrap();
        let values = Array2::from_shape_vec((n, 1), self.values.clone().to_vec()).unwrap();

        (
            self.wind_directions.clone(),
            self.wind_speeds.clone(),
            self.turbulence_intensities.clone(),
            freq,
            values,
            self.heterogeneous_inflow_config(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;

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
        approx::assert_relative_eq!(ts.turbulence_intensities[2], 0.07, epsilon = 0.01);
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
}
