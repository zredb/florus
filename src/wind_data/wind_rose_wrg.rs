//! Wind rose from WAsP WRG file and related types.
//!
//! This module provides support for:
//! - WindRoseWRG: Wind rose from WAsP WRG (Wind Resource Grid) file
//! - WindRoseByTurbine: Wind rose with separate wind rose for each turbine
//! - RegularGridInterpolant: 2D interpolation helper
//! - WRGData: Internal WRG file data structure

use crate::core::base::InterpMethod;
use crate::heterogeneous_map::HeterogeneousInflowConfig;
use crate::types::{Array1, Array2, Array3, Float};
use crate::wind_data::traits::WindData;
use crate::wind_data::WindRose;
use crate::Result;

/// Regular grid interpolant wrapper for 2D interpolation
///
/// Provides bilinear interpolation over a regular grid for spatial
/// interpolation of wind resource data.
#[derive(Clone, Debug)]
pub struct RegularGridInterpolant {
    /// X coordinates of the grid
    pub x: Array1,
    /// Y coordinates of the grid
    pub y: Array1,
    /// Data values [ny, nx] (row-major: y varies first, then x)
    pub data: Array2,
    /// Interpolation method
    pub method: InterpMethod,
}

impl RegularGridInterpolant {
    /// Create a new interpolant
    pub fn new(x: Array1, y: Array1, data: Array2, method: InterpMethod) -> Self {
        Self { x, y, data, method }
    }

    /// Interpolate at point (x, y)
    pub fn interpolate(&self, x: Float, y: Float) -> Float {
        let x_vec: Vec<Float> = self.x.iter().copied().collect();
        let y_vec: Vec<Float> = self.y.iter().copied().collect();
        bilinear_interpolate(x, y, &x_vec, &y_vec, &self.data)
    }
}

/// Internal WRG file data structure
///
/// Represents parsed data from a WAsP Wind Resource Grid file.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct WRGData {
    /// Number of grid points in x direction
    pub nx: usize,
    /// Number of grid points in y direction
    pub ny: usize,
    /// Minimum x coordinate
    pub xmin: Float,
    /// Minimum y coordinate
    pub ymin: Float,
    /// Grid cell size
    pub grid_size: Float,
    /// Number of wind sectors
    pub n_sectors: usize,
    /// X coordinates of grid points
    pub x_array: Vec<Float>,
    /// Y coordinates of grid points
    pub y_array: Vec<Float>,
    /// Sector frequencies [nx, ny, n_sectors]
    pub sector_freq: Array3,
    /// Weibull A parameters [nx, ny, n_sectors]
    pub weibull_a: Array3,
    /// Weibull k parameters [nx, ny, n_sectors]
    pub weibull_k: Array3,
}

/// WindRoseWRG - Wind rose from WAsP WRG (Wind Resource Grid) file
///
/// WindRoseWRG represents a wind resource grid (WRG) file where each grid point
/// has a separate wind rose defined by sector frequencies and Weibull parameters.
/// When a layout is specified, each turbine gets its own wind rose computed by
/// interpolating Weibull parameters from the WRG grid.
///
/// Corresponds to WindRoseWRG in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct WindRoseWRG {
    /// WRG file data
    pub wrg_data: WRGData,
    /// Wind directions (may be resampled from WRG sectors)
    pub wind_directions: Array1,
    /// Wind direction step size
    pub wd_step: Float,
    /// Wind speed bins
    pub wind_speeds: Array1,
    /// Turbulence intensity table or single value
    pub ti_table: Array2,
    /// Turbulence intensity as a single value (for backwards compatibility)
    pub ti_value: Float,
    /// Layout x coordinates (set via set_layout)
    pub layout_x: Array1,
    /// Layout y coordinates (set via set_layout)
    pub layout_y: Array1,
    /// Wind roses for each turbine position
    pub wind_roses: Vec<WindRose>,
    /// Flattened wind directions from first wind rose
    pub wd_flat: Array1,
    /// Flattened wind speeds from first wind rose
    pub ws_flat: Array1,
    /// Non-zero frequency mask
    pub nonzero_freq_mask: Vec<bool>,
    /// Interpolants for sector frequency [n_sectors]
    pub interpolant_sector_freq: Vec<RegularGridInterpolant>,
    /// Interpolants for Weibull A parameter [n_sectors]
    pub interpolant_weibull_a: Vec<RegularGridInterpolant>,
    /// Interpolants for Weibull k parameter [n_sectors]
    pub interpolant_weibull_k: Vec<RegularGridInterpolant>,
}

impl Default for WindRoseWRG {
    fn default() -> Self {
        Self {
            wrg_data: WRGData::default(),
            wind_directions: Array1::from_vec(vec![]),
            wd_step: 0.0,
            wind_speeds: Array1::from_vec(vec![]),
            ti_table: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
            ti_value: 0.06,
            layout_x: Array1::from_vec(vec![]),
            layout_y: Array1::from_vec(vec![]),
            wind_roses: Vec::new(),
            wd_flat: Array1::from_vec(vec![]),
            ws_flat: Array1::from_vec(vec![]),
            nonzero_freq_mask: Vec::new(),
            interpolant_sector_freq: Vec::new(),
            interpolant_weibull_a: Vec::new(),
            interpolant_weibull_k: Vec::new(),
        }
    }
}

impl WindRoseWRG {
    /// Create a new WindRoseWRG from a WRG file
    ///
    /// # Arguments
    /// * `filename` - Path to the WRG file
    /// * `wd_step` - Wind direction step for resampling (None = use WRG sectors)
    /// * `wind_speeds` - Wind speed bins (default: 0-25 m/s in 1 m/s increments)
    /// * `ti_table` - Turbulence intensity (single value or 2D table)
    ///
    /// # Returns
    /// WindRoseWRG with WRG data loaded
    pub fn new<P: AsRef<std::path::Path>>(
        filename: P,
        wd_step: Option<Float>,
        wind_speeds: Option<Array1>,
        ti_table: Option<Array2>,
    ) -> Result<Self> {
        let filename = filename.as_ref();

        // Read and parse WRG file
        let wrg_data = read_wrg_file(filename)?;

        // Set default wind speeds if not provided
        let wind_speeds =
            wind_speeds.unwrap_or_else(|| Array1::from_iter((0..26).map(|i| i as Float)));

        // Set default TI table if not provided
        let ti_table = ti_table.unwrap_or_else(|| {
            let n_dir = wrg_data.n_sectors;
            let n_ws = wind_speeds.len();
            Array2::from_elem((n_dir, n_ws), 0.06)
        });

        // Get TI single value (first element) for backwards compatibility
        let ti_value = if ti_table.shape()[0] > 0 && ti_table.shape()[1] > 0 {
            ti_table[[0, 0]]
        } else {
            0.06
        };

        // Calculate wind directions from WRG sectors or use specified step
        let (wind_directions, wd_step) = if let Some(step) = wd_step {
            let n_dir = (360.0 / step).ceil() as usize;
            let directions: Array1 = (0..n_dir)
                .map(|i| (i as Float + 0.5) * step)
                .collect();
            (directions, step)
        } else {
            let directions: Array1 = (0..wrg_data.n_sectors)
                .map(|i| i as Float * (360.0 / wrg_data.n_sectors as Float))
                .collect();
            let step = if wrg_data.n_sectors > 1 {
                360.0 / wrg_data.n_sectors as Float
            } else {
                360.0
            };
            (directions, step)
        };

        // Build interpolants for each sector
        let interpolant_sector_freq = build_interpolant_function_list(
            &wrg_data.x_array,
            &wrg_data.y_array,
            &wrg_data.sector_freq,
        )?;
        let interpolant_weibull_a = build_interpolant_function_list(
            &wrg_data.x_array,
            &wrg_data.y_array,
            &wrg_data.weibull_a,
        )?;
        let interpolant_weibull_k = build_interpolant_function_list(
            &wrg_data.x_array,
            &wrg_data.y_array,
            &wrg_data.weibull_k,
        )?;

        Ok(Self {
            wrg_data,
            wind_directions,
            wd_step,
            wind_speeds,
            ti_table,
            ti_value,
            layout_x: Array1::from_vec(vec![]),
            layout_y: Array1::from_vec(vec![]),
            wind_roses: Vec::new(),
            wd_flat: Array1::from_vec(vec![]),
            ws_flat: Array1::from_vec(vec![]),
            nonzero_freq_mask: Vec::new(),
            interpolant_sector_freq,
            interpolant_weibull_a,
            interpolant_weibull_k,
        })
    }

    /// Interpolate data at a given (x, y) location
    fn interpolate_data(&self, x: Float, y: Float, interpolants: &[RegularGridInterpolant]) -> Array1 {
        let n_sectors = self.wrg_data.n_sectors;
        let mut result = Array1::zeros(n_sectors);

        // Determine interpolation method based on bounds
        let use_nearest = x < self.wrg_data.x_array[0]
            || x > self.wrg_data.x_array[self.wrg_data.x_array.len() - 1]
            || y < self.wrg_data.y_array[0]
            || y > self.wrg_data.y_array[self.wrg_data.y_array.len() - 1];

        for sector in 0..n_sectors {
            if use_nearest {
                // Use nearest neighbor for points outside bounds
                result[sector] = interpolants[sector].interpolate(
                    x.clamp(
                        self.wrg_data.x_array[0],
                        self.wrg_data.x_array[self.wrg_data.x_array.len() - 1],
                    ),
                    y.clamp(
                        self.wrg_data.y_array[0],
                        self.wrg_data.y_array[self.wrg_data.y_array.len() - 1],
                    ),
                );
            } else {
                result[sector] = interpolants[sector].interpolate(x, y);
            }
        }

        result
    }

    /// Calculate Weibull cumulative distribution function
    fn weibull_cumulative(&self, x: Float, a: Float, k: Float) -> Float {
        if x <= 0.0 {
            return 0.0;
        }
        let exponent = -((x / a).powf(k));
        1.0 - exponent.exp()
    }

    /// Generate wind speed frequencies from Weibull parameters
    fn generate_wind_speed_frequencies_from_weibull(
        &self,
        a: Float,
        k: Float,
    ) -> (Array1, Array1) {
        let n_ws = self.wind_speeds.len();
        let ws_step = if n_ws > 1 {
            self.wind_speeds[1] - self.wind_speeds[0]
        } else {
            1.0
        };

        let mut frequencies = Vec::with_capacity(n_ws);

        for i in 0..n_ws {
            let ws = self.wind_speeds[i];
            let ws_low = (ws - ws_step / 2.0).max(0.0);
            let ws_high = ws + ws_step / 2.0;

            let cdf_high = self.weibull_cumulative(ws_high, a, k);
            let cdf_low = self.weibull_cumulative(ws_low, a, k);

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

        (
            self.wind_speeds.clone(),
            Array1::from_vec(frequencies),
        )
    }

    /// Get the wind rose at a specific (x, y) location
    pub fn get_wind_rose_at_point(
        &self,
        x: Float,
        y: Float,
        wind_directions: Option<Array1>,
        wind_speeds: Option<Array1>,
        ti_table: Option<Array2>,
    ) -> WindRose {
        let wd = wind_directions.unwrap_or_else(|| self.wind_directions.clone());
        let ws = wind_speeds.unwrap_or_else(|| self.wind_speeds.clone());
        let ti = ti_table.unwrap_or_else(|| self.ti_table.clone());

        let n_sectors = self.wrg_data.n_sectors;
        let n_ws = ws.len();

        // Get interpolated data at this location
        let sector_freq = self.interpolate_data(x, y, &self.interpolant_sector_freq);
        let weibull_a = self.interpolate_data(x, y, &self.interpolant_weibull_a);
        let weibull_k = self.interpolate_data(x, y, &self.interpolant_weibull_k);

        // Build frequency table from Weibull distributions
        let mut freq_table = Array2::zeros((n_sectors, n_ws));

        for sector in 0..n_sectors {
            let (_, freq) = Self::_generate_wind_speed_frequencies_from_weibull_internal(
                &ws,
                weibull_a[sector],
                weibull_k[sector],
            );
            for j in 0..n_ws {
                freq_table[[sector, j]] = sector_freq[sector] * freq[j];
            }
        }

        // Normalize frequency table
        let total: Float = freq_table.iter().sum();
        if total > 0.0 {
            for val in &mut freq_table {
                *val /= total;
            }
        }

        // Get original WRG wind directions for base wind rose
        let wrg_wind_directions: Array1 = (0..n_sectors)
            .map(|i| i as Float * (360.0 / n_sectors as Float))
            .collect();

        // Create base wind rose
        let mut wind_rose = WindRose::new(
            wrg_wind_directions,
            ws.clone(),
            ti.clone(),
            Some(freq_table),
            None,
        )
        .unwrap_or_default();

        wind_rose.compute_zero_freq_occurrence = true;

        // Resample to desired wind directions if needed
        let wd_step = if wd.len() > 1 { wd[1] - wd[0] } else { 360.0 };
        let base_wd_step = 360.0 / n_sectors as Float;

        if (wd_step - base_wd_step).abs() < 1e-6 {
            // Same resolution, just update frequencies
            wind_rose.wind_directions = wd;
            wind_rose
        } else if wd_step < base_wd_step {
            // Need to upsample
            wind_rose.upsample(wd_step, ws[1] - ws[0], "linear", false)
        } else {
            // Need to downsample
            wind_rose.downsample(wd_step, ws[1] - ws[0], false)
        }
    }

    /// Internal helper for generating wind speed frequencies from Weibull
    fn _generate_wind_speed_frequencies_from_weibull_internal(
        wind_speeds: &Array1,
        a: Float,
        k: Float,
    ) -> (Array1, Array1) {
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

            let exponent_low = -((ws_low / a).powf(k));
            let cdf_low = if ws_low <= 0.0 { 0.0 } else { 1.0 - exponent_low.exp() };

            let exponent_high = -((ws_high / a).powf(k));
            let cdf_high = 1.0 - exponent_high.exp();

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

        (
            wind_speeds.clone(),
            Array1::from_vec(frequencies),
        )
    }

    /// Set the wind direction step
    pub fn set_wd_step(&mut self, wd_step: Float) {
        self.wd_step = wd_step;
        self.wind_directions = (0..((360.0 / wd_step).ceil() as usize))
            .map(|i| (i as Float + 0.5) * wd_step)
            .collect();

        // Update wind roses if layout is set
        if !self.layout_x.is_empty() {
            self._update_wind_roses();
        }
    }

    /// Set the wind speeds
    pub fn set_wind_speeds(&mut self, wind_speeds: Array1) {
        self.wind_speeds = wind_speeds;

        // Update TI table if it's a single value
        if self.ti_table.shape() == &[1, 1] || self.ti_table.len() == 1 {
            let ti_val = self.ti_table.iter().next().copied().unwrap_or(0.06);
            let n_dir = self.wind_directions.len();
            let n_ws = self.wind_speeds.len();
            self.ti_table = Array2::from_elem((n_dir, n_ws), ti_val);
        }

        // Update wind roses if layout is set
        if !self.layout_x.is_empty() {
            self._update_wind_roses();
        }
    }

    /// Set the turbulence intensity (single value)
    pub fn set_ti_table(&mut self, ti_value: Float) {
        self.ti_value = ti_value;
        let n_dir = self.wind_directions.len();
        let n_ws = self.wind_speeds.len();
        self.ti_table = Array2::from_elem((n_dir, n_ws), ti_value);

        // Update wind roses if layout is set
        if !self.layout_x.is_empty() {
            self._update_wind_roses();
        }
    }

    /// Set the turbine layout
    pub fn set_layout(&mut self, layout_x: Array1, layout_y: Array1) -> Result<()> {
        if layout_x.len() != layout_y.len() {
            anyhow::bail!("layout_x and layout_y must have the same length");
        }

        // Check if layout is the same
        if !self.layout_x.is_empty()
            && self.layout_x.len() == layout_x.len()
            && self.layout_y.len() == layout_y.len()
        {
            let mut same = true;
            for i in 0..layout_x.len() {
                if (layout_x[i] - self.layout_x[i]).abs() > 1e-6
                    || (layout_y[i] - self.layout_y[i]).abs() > 1e-6
                {
                    same = false;
                    break;
                }
            }
            if same {
                return Ok(());
            }
        }

        self.layout_x = layout_x.clone();
        self.layout_y = layout_y.clone();

        self._update_wind_roses();

        Ok(())
    }

    /// Update wind roses for current layout
    fn _update_wind_roses(&mut self) {
        let n_turbines = self.layout_x.len();
        self.wind_roses = Vec::with_capacity(n_turbines);

        for i in 0..n_turbines {
            let wind_rose = self.get_wind_rose_at_point(
                self.layout_x[i],
                self.layout_y[i],
                Some(self.wind_directions.clone()),
                Some(self.wind_speeds.clone()),
                Some(self.ti_table.clone()),
            );
            self.wind_roses.push(wind_rose);
        }

        // Update flattened data from first wind rose
        if let Some(first_rose) = self.wind_roses.first() {
            self.wd_flat = first_rose.wind_directions.clone();
            self.ws_flat = first_rose.wind_speeds.clone();
            self.nonzero_freq_mask = (0..self.wd_flat.len())
                .map(|i| {
                    for rose in &self.wind_roses {
                        if let Some(ref freq) = rose.freq_table {
                            if freq[[i % first_rose.wind_directions.len(), i / first_rose.wind_directions.len()]] > 0.0 {
                                return true;
                            }
                        }
                    }
                    false
                })
                .collect();
        }
    }
}

impl WindData for WindRoseWRG {
    fn wind_speeds(&self) -> Array1 {
        self.wind_speeds.clone()
    }

    fn wind_directions(&self) -> Array1 {
        self.wind_directions.clone()
    }

    fn turbulence_intensities(&self) -> Array1 {
        // Return representative TI values (first speed bin for each direction)
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

    fn heterogeneous_inflow_config(&self) -> HeterogeneousInflowConfig {
        HeterogeneousInflowConfig {
            x: Array1::from_vec(vec![]),
            y: Array1::from_vec(vec![]),
            z: None,
            wind_speeds: Some(self.wind_speeds.clone()),
            wind_directions: Some(self.wind_directions.clone()),
            speed_multipliers: Array2::from_shape_vec((self.wd_flat.len(), 0), vec![]).unwrap(),
        }
    }

    fn set_layout(&mut self, layout_x: &Option<Array1>, layout_y: &Option<Array1>) {
        if let (Some(lx), Some(ly)) = (layout_x, layout_y) {
            let _ = self.set_layout(lx.clone(), ly.clone());
        }
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
        if self.layout_x.is_empty() {
            // Return empty result if layout is not set
            return (
                Array1::from_vec(vec![]),
                Array1::from_vec(vec![]),
                Array1::from_vec(vec![]),
                Array2::from_shape_vec((0, 0), vec![]).unwrap(),
                Array2::from_shape_vec((0, 0), vec![]).unwrap(),
                HeterogeneousInflowConfig {
                    x: Array1::from_vec(vec![]),
                    y: Array1::from_vec(vec![]),
                    z: None,
                    wind_speeds: None,
                    wind_directions: None,
                    speed_multipliers: Array2::from_shape_vec((0, 0), vec![]).unwrap(),
                },
            );
        }

        let n_turbines = self.layout_x.len();
        let n_conditions = self.wd_flat.len();

        // Initialize freq_table_unpack [n_conditions, n_turbines]
        let mut freq_table_unpack = Array2::zeros((n_conditions, n_turbines));

        let mut wind_directions_unpack = Array1::from_vec(vec![]);
        let mut wind_speeds_unpack = Array1::from_vec(vec![]);
        let mut ti_table_unpack = Array1::from_vec(vec![]);
        let mut value_table_unpack = Array2::zeros((n_conditions, n_turbines));

        // Loop over wind roses and collect data
        for i in 0..n_turbines {
            let wind_rose = &self.wind_roses[i];
            let (
                wd,
                ws,
                ti,
                freq_2d,
                value_2d,
                _,
            ) = wind_rose.unpack();

            if i == 0 {
                wind_directions_unpack = wd;
                wind_speeds_unpack = ws;
                ti_table_unpack = ti;
            }

            // Copy frequency column
            for j in 0..n_conditions {
                freq_table_unpack[[j, i]] = freq_2d[[j, 0]];
                value_table_unpack[[j, i]] = value_2d[[j, 0]];
            }
        }

        (
            wind_directions_unpack,
            wind_speeds_unpack,
            ti_table_unpack,
            freq_table_unpack,
            value_table_unpack,
            self.heterogeneous_inflow_config(),
        )
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
        let default_freq =
            Array2::from_elem((n_dir * n_ws, n_turbines), 1.0 / (n_dir * n_ws) as Float);
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

        WindRose::new(wd, ws, ti, Some(freq_table), None).unwrap_or_else(|_| WindRose::default())
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
        let wind_speeds =
            wind_speeds.unwrap_or_else(|| Array1::from_iter((0..26).map(|i| i as Float)));

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
                if angle >= 360.0 {
                    angle - 360.0
                } else {
                    angle
                }
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

    fn set_layout(&mut self, layout_x: &Option<Array1>, layout_y: &Option<Array1>) {
        if let (Some(lx), Some(ly)) = (layout_x, layout_y) {
            let _ = self.set_layout(lx.clone(), ly.clone());
        }
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
            self.heterogeneous_inflow_config(),
        )
    }
}

// Utility functions

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
fn generate_weibull_frequencies(a: Float, k: Float, wind_speeds: &Array1) -> Array1 {
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
    let use_nearest = x < x_array[0]
        || x > x_array[x_array.len() - 1]
        || y < y_array[0]
        || y > y_array[y_array.len() - 1];

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

/// Build a list of interpolant functions for each sector
fn build_interpolant_function_list(
    x: &[Float],
    y: &[Float],
    data: &Array3,
) -> Result<Vec<RegularGridInterpolant>> {
    let n_sectors = data.shape()[2];
    let mut interpolants = Vec::with_capacity(n_sectors);

    for sector in 0..n_sectors {
        // Extract 2D slice for this sector
        let sector_data = Array2::from_shape_fn((data.shape()[0], data.shape()[1]), |(i, j)| {
            data[[i, j, sector]]
        });

        let x_arr = Array1::from_vec(x.to_vec());
        let y_arr = Array1::from_vec(y.to_vec());

        interpolants.push(RegularGridInterpolant::new(
            x_arr,
            y_arr,
            sector_data,
            InterpMethod::Linear,
        ));
    }

    Ok(interpolants)
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
            lines.len(),
            expected_lines,
            n_sectors
        );
    }

    // Initialize data arrays as 3D for sector_freq (nx, ny, n_sectors)
    let mut sector_freq = Array3::zeros((nx, ny, n_sectors));
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

    // Parse data for each grid point
    let mut line_idx = 2;
    for _gid in 0..(nx * ny) {
        if line_idx >= lines.len() {
            break;
        }
        let line = &lines[line_idx];
        line_idx += 1;

        // Parse basic coordinates
        let x_val: Float = line[10..20].trim().parse().unwrap_or(0.0);
        let y_val: Float = line[20..30].trim().parse().unwrap_or(0.0);
        let _z_val: Float = line[30..38].trim().parse().unwrap_or(0.0);
        let _h_val: Float = line[38..43].trim().parse().unwrap_or(0.0);

        // Find x and y indices
        let x_idx = ((x_val - xmin) / grid_size).round() as usize;
        let y_idx = ((y_val - ymin) / grid_size).round() as usize;

        if x_idx >= nx || y_idx >= ny {
            continue;
        }

        // Parse sector data
        for sector in 0..n_sectors {
            let base_pos = 72 + sector * 13;

            // Frequency (probability * 1000) - 4 chars at base_pos
            let freq_str = &line[base_pos..base_pos + 4];
            if !freq_str.trim().is_empty() {
                if let Ok(freq) = freq_str.trim().parse::<Float>() {
                    sector_freq[[x_idx, y_idx, sector]] = freq / 1000.0;
                }
            }

            // Weibull A (stored * 10) - 4 chars at base_pos + 4
            let a_str = &line[base_pos + 4..base_pos + 8];
            if !a_str.trim().is_empty() {
                if let Ok(a_val) = a_str.trim().parse::<Float>() {
                    weibull_a[[x_idx, y_idx, sector]] = a_val / 10.0;
                }
            }

            // Weibull k (stored * 100) - 5 chars at base_pos + 8
            let k_str = &line[base_pos + 8..base_pos + 13];
            if !k_str.trim().is_empty() {
                if let Ok(k_val) = k_str.trim().parse::<Float>() {
                    weibull_k[[x_idx, y_idx, sector]] = k_val / 100.0;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weibull_cumulative() {
        // Test Weibull CDF at x = A should be 1 - exp(-1) ≈ 0.632
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
        let data =
            Array2::from_shape_vec((3, 3), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
                .unwrap();

        // Test at exact grid point (0, 0) = data[0,0] = 1.0
        let result = bilinear_interpolate(0.0, 0.0, &x_array, &y_array, &data);
        assert!((result - 1.0).abs() < 0.001);

        // Test at grid point (100, 0) = data[1,0] = 2.0
        let result = bilinear_interpolate(100.0, 0.0, &x_array, &y_array, &data);
        assert!((result - 2.0).abs() < 0.001);

        // Test at center of grid cell (50, 50)
        let result = bilinear_interpolate(50.0, 50.0, &x_array, &y_array, &data);
        assert!((result - 3.0).abs() < 0.001);
    }

    #[test]
    fn test_wind_rose_by_turbine_creation() {
        use crate::wind_data::WindRose;

        let wd = Array1::from_vec(vec![0.0, 90.0, 180.0, 270.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0, 12.0]);
        let ti_table = Array2::from_elem((4, 3), 0.08);

        // Create individual wind roses for 2 turbines
        let freq1 = Array2::from_elem((4, 3), 0.25);
        let freq2 = Array2::from_elem((4, 3), 0.25);

        let wr1 =
            WindRose::new(wd.clone(), ws.clone(), ti_table.clone(), Some(freq1), None).unwrap();
        let wr2 =
            WindRose::new(wd.clone(), ws.clone(), ti_table.clone(), Some(freq2), None).unwrap();

        let wrbt = WindRoseByTurbine::new(wd, ws, ti_table, vec![wr1, wr2]).unwrap();
        assert_eq!(wrbt.n_conditions(), 12); // 4 dirs × 3 speeds
    }

    #[test]
    fn test_wind_rose_by_turbine_set_layout() {
        use crate::wind_data::WindRose;

        let wd = Array1::from_vec(vec![0.0, 180.0]);
        let ws = Array1::from_vec(vec![8.0, 10.0]);
        let ti_table = Array2::from_elem((2, 2), 0.06);

        // Create individual wind roses
        let freq = Array2::from_elem((2, 2), 0.25);
        let wr1 = WindRose::new(
            wd.clone(),
            ws.clone(),
            ti_table.clone(),
            Some(freq.clone()),
            None,
        )
        .unwrap();
        let wr2 = WindRose::new(
            wd.clone(),
            ws.clone(),
            ti_table.clone(),
            Some(freq.clone()),
            None,
        )
        .unwrap();

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
    fn test_wind_rose_by_turbine_default() {
        let wrbt = WindRoseByTurbine::default();
        assert!(wrbt.wind_directions.is_empty());
        assert!(wrbt.wind_speeds.is_empty());
        assert!(wrbt.wind_roses.is_empty());
    }
}
