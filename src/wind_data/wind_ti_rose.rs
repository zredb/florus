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

    fn n_conditions(&self) -> usize {
        self.wind_directions.len() * self.wind_speeds.len() * self.turbulence_intensities.len()
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

    #[test]
    fn test_wind_ti_rose_creation() {
        let wd = Array1::from_vec(vec![0.0, 180.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti = Array1::from_vec(vec![0.06, 0.08, 0.10]);
        let ti_table =
            Array3::from_shape_fn((2, 2, 3), |(i, j, k)| 0.05 + (i + j + k) as Float * 0.01);

        let wir = WindTIRose::new(wd, ws, ti, ti_table, None, None).unwrap();
        assert_eq!(wir.n_conditions(), 12); // 2 × 2 × 3
    }
}
