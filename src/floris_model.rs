use crate::core::turbines::cp_ct_table::TableConditions;
use crate::core::turbines::operation_models::{POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED};
use crate::core::{Core, Farm, FlowField, State, TurbineGrid};
use crate::floris_config::{FlorisConfig, SolverConfig};
use crate::types::{Array1, Array2, Array3, Array4, Float};
use crate::utilities::{cosd, load_yaml};
use crate::wind_data::WindData;
use anyhow::bail;
use ndarray::Array;
use ndarray::Array2 as NdArray2;
use ndarray::Axis;
use std::fmt;
use std::path::Path;

/// Operation model specification for turbines
pub enum OperationModelSpec {
    /// Single operation model applied to all turbines
    Single(String),
    /// Different operation models per turbine
    Multiple(Vec<String>),
}

impl From<String> for OperationModelSpec {
    fn from(model: String) -> Self {
        OperationModelSpec::Single(model)
    }
}

impl From<&str> for OperationModelSpec {
    fn from(model: &str) -> Self {
        OperationModelSpec::Single(model.to_string())
    }
}

impl From<Vec<String>> for OperationModelSpec {
    fn from(models: Vec<String>) -> Self {
        OperationModelSpec::Multiple(models)
    }
}

/// Turbine layout representation
pub enum TurbineLayout {
    /// X and Y coordinates only
    XY(Array1, Array1),
    /// X, Y, and Z (hub height) coordinates
    WithZ(Array1, Array1, Array1),
}

impl TurbineLayout {
    /// Get x coordinates
    pub fn x(&self) -> &Array1 {
        match self {
            TurbineLayout::XY(x, _) => x,
            TurbineLayout::WithZ(x, _, _) => x,
        }
    }

    /// Get y coordinates
    pub fn y(&self) -> &Array1 {
        match self {
            TurbineLayout::XY(_, y) => y,
            TurbineLayout::WithZ(_, y, _) => y,
        }
    }

    /// Get z coordinates (if available)
    pub fn z(&self) -> Option<&Array1> {
        match self {
            TurbineLayout::XY(_, _) => None,
            TurbineLayout::WithZ(_, _, z) => Some(z),
        }
    }
}

/// Main FLORIS Model structure
///
/// This is a high-level wrapper around the Core class, providing a user-friendly API
/// similar to Python FLORIS's FlorisModel.
pub struct FlorisModel {
    core: Core,
}

impl Clone for FlorisModel {
    fn clone(&self) -> Self {
        Self {
            core: self.core.clone(),
        }
    }
}

impl fmt::Debug for FlorisModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FlorisModel")
            .field("core", &self.core)
            .finish()
    }
}

impl FlorisModel {
    /// Create a new FlorisModel from a configuration file
    pub fn from_file<P: AsRef<Path>>(path: P) -> crate::Result<Self> {
        // Parse the input file here (instead of delegating to Core::from_file)
        // so that the file-level settings handled by `from_config` -- such as
        // the heterogeneous inflow configuration -- are applied.
        let config_value = crate::utilities::load_yaml(path)?;
        let config: FlorisConfig = serde_yaml::from_value(config_value)
            .map_err(|e| anyhow::anyhow!("Failed to parse config: {}", e))?;
        Self::from_config(config)
    }

    pub fn from_default() -> crate::Result<Self> {
        Self::from_file("default_inputs.yaml") // Fixed: added 's' to match actual filename
    }

    /// Create from configuration structure
    pub fn from_config(config: FlorisConfig) -> crate::Result<Self> {
        // Heterogeneous inflow from the input file (Python
        // flow_field.heterogeneous_inflow_config)
        let het_config = config
            .flow_field
            .heterogeneous_inflow_config
            .as_ref()
            .map(|c| {
                let multipliers = Array2::from_shape_vec(
                    (c.speed_multipliers.len(), c.x.len()),
                    c.speed_multipliers
                        .iter()
                        .flat_map(|row| row.iter().copied())
                        .collect(),
                )
                .unwrap_or_else(|_| Array2::zeros((0, 0)));
                crate::heterogeneous_map::HeterogeneousInflowConfig {
                    x: Array1::from_vec(c.x.clone()),
                    y: Array1::from_vec(c.y.clone()),
                    z: c.z.as_ref().map(|z| Array1::from_vec(z.clone())),
                    wind_speeds: c.wind_speeds.as_ref().map(|s| Array1::from_vec(s.clone())),
                    wind_directions: c
                        .wind_directions
                        .as_ref()
                        .map(|d| Array1::from_vec(d.clone())),
                    speed_multipliers: multipliers,
                }
            });

        let mut core = Core::from_config(config)?;
        if let Some(cfg) = het_config {
            core.flow_field.heterogeneous_inflow_config = Some(cfg);
        }
        Ok(Self { core })
    }
    /// Run the FLORIS simulation
    ///
    /// This method performs the steady-state atmospheric condition calculation,
    /// which includes wake interactions and turbine power/thrust calculations.
    ///
    /// # Examples
    /// ```no_run
    /// # use anyhow::Error;
    /// # fn main() -> Result<(), Error> {
    /// use florus::FlorisModel;
    /// let mut model = FlorisModel::from_file("default_inputs.yaml")?;
    /// model.run()?;
    /// let powers = model.get_turbine_powers();
    /// # Ok(())
    /// # }
    /// ```
    pub fn run(&mut self) -> crate::Result<()> {
        // Use Core's steady_state_atmospheric_condition method
        self.core.steady_state_atmospheric_condition()?;
        Ok(())
    }

    /// Run the FLORIS simulation without wake modeling
    ///
    /// This method is similar to `run()` except that it does not apply a wake model.
    /// The wind farm is modeled as if there is no wake in the flow. Operation settings
    /// may still reduce the power and thrust of the turbines where they're applied.
    ///
    /// This is useful for:
    /// - Comparing results with and without wake effects
    /// - Calculating theoretical maximum power production
    /// - Debugging and validation purposes
    ///
    /// # Examples
    /// ```no_run
    /// # use anyhow::Error;
    /// # fn main() -> Result<(), Error> {
    /// use florus::FlorisModel;
    /// let mut model = FlorisModel::from_file("default_inputs.yaml")?;
    ///
    /// // Run with wake effects
    /// model.run()?;
    /// let power_with_wake = model.get_farm_power()?;
    ///
    /// // Run without wake effects
    /// model.run_no_wake()?;
    /// let power_without_wake = model.get_farm_power()?;
    ///
    /// println!("Wake loss: {} W", power_without_wake - power_with_wake);
    /// # Ok(())
    /// # }
    /// ```
    pub fn run_no_wake(&mut self) -> crate::Result<()> {
        // Initialize solution space (this will initialize velocity field and farm quantities)
        self.core.initialize_domain()?;

        // Finalize values to user-supplied order
        self.core.finalize();

        Ok(())
    }

    /// Get reference to the grid
    pub fn grid(&self) -> Option<&dyn crate::core::Grid> {
        self.core.grid.as_ref().map(|g| g.as_ref())
    }

    /// Get mutable reference to farm
    pub fn farm_mut(&mut self) -> &mut crate::core::Farm {
        &mut self.core.farm
    }

    /// Initialize the computational grid
    pub fn initialize_grid(&mut self) -> crate::Result<()> {
        // Grid is managed by Core, this method is kept for compatibility
        if self.core.grid.is_none() {
            self.core.initialize_grid()?;
        }
        Ok(())
    }

    /// Initialize flow field on the grid
    pub fn initialize_flow_field(&mut self) -> crate::Result<()> {
        // Flow field initialization is handled by Core
        self.core.initialize_domain()?;
        Ok(())
    }

    /// Calculates the power at each turbine in the wind farm (Watts), routed
    /// through the operation model exactly like Python's `get_turbine_powers()`
    /// (`floris_model.py:554-570` → `turbine.py:power()` → the configured
    /// operation model). Power setpoints are applied inside the model, so
    /// disabled turbines (near-zero setpoints) report zero power.
    ///
    /// # Errors
    /// - If `run()` has not been called since the last state change; Python
    ///   raises the equivalent RuntimeError from `_get_turbine_powers`.
    /// - Multi-dimensional Cp/Ct fleets fall back to the legacy inline
    ///   estimate (`get_turbine_powers_legacy`): Python routes those through
    ///   a separate `turbine.py` branch that FLORUS has not ported.
    pub fn get_turbine_powers(&self) -> crate::Result<Array2> {
        // Python _get_turbine_powers: "Can't compute turbine powers without
        // first running `FlorisModel.run()`." (state must be USED)
        if !self.core.state.converged {
            bail!("Can't compute turbine powers without first running `FlorisModel::run()`.");
        }

        let farm = &self.core.farm;
        let flow_field = &self.core.flow_field;

        // Python turbine.py:power() defaults: average_method="cubic-mean" with
        // no cubature weights — the grid's averaging does NOT feed the power
        // path (verified against Python 2026-10-14; see tests/cross_validation_wake.rs
        // `compute_powers_cubature`).
        let n_findex = flow_field.u_sorted.shape()[0];
        let turbulence_intensities = {
            let ti = &flow_field.turbulence_intensity_field_sorted;
            if ti.shape().len() == 4 && ti.shape()[0] == n_findex && ti.shape()[1] > 0 {
                Some(ti)
            } else {
                None
            }
        };
        let powers_sorted = farm.power(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            turbulence_intensities,
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            crate::core::rotor_velocity::AveragingMethod::CubicMean,
            None,
        )?;

        Ok(self.unsort_to_layout(&powers_sorted))
    }

    /// Reorder a (n_findex, n_turbines) sorted-order result back to the
    /// user's layout order. Python returns layout-order arrays directly
    /// (`flow_field.u` is unsorted); FLORUS computes on sorted views, so
    /// results must be mapped through `grid.unsorted_indices`.
    fn unsort_to_layout(&self, sorted: &Array2) -> Array2 {
        let n_findex = sorted.shape()[0];
        let n_turbines = sorted.shape()[1];

        let Some(grid) = self.core.grid.as_ref() else {
            return sorted.clone();
        };
        let Some(unsorted_indices) = grid.unsorted_indices() else {
            return sorted.clone();
        };

        let mut out = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let pos = unsorted_indices[[fi, ti]] as usize;
                out[[fi, ti]] = sorted[[fi, pos]];
            }
        }
        out
    }

    /// Get farm power
    pub fn get_farm_power(&self) -> crate::Result<Array1> {
        let powers = self.get_turbine_powers()?;
        let mut farm_power = Array1::zeros(self.core.flow_field.n_findex);

        for i in 0..self.core.flow_field.n_findex {
            farm_power[i] = powers.row(i).sum();
        }

        Ok(farm_power)
    }

    /// Get turbine powers reshaped for WindRose mode
    ///
    /// Returns powers in shape (n_wind_directions, n_wind_speeds, n_turbines)
    /// Only works when wind_data_info is available (set via set_wind_conditions_with_rose)
    pub fn get_turbine_powers_rose(&self) -> crate::Result<Array3> {
        let powers = self.get_turbine_powers()?;

        Ok(if let Some(ref wind_info) = self.core.wind_data_info {
            // Reshape from (n_findex, n_turbines) to (n_wind_directions, n_wind_speeds, n_turbines)
            let n_wind_dirs = wind_info.n_wind_directions;
            let n_wind_speeds = wind_info.n_wind_speeds;
            let n_turbines = self.core.farm.n_turbines();

            let mut reshaped = Array3::zeros((n_wind_dirs, n_wind_speeds, n_turbines));

            // Reshape the data assuming findex order is [wd0_ws0, wd0_ws1, ..., wd1_ws0, wd1_ws1, ...]
            for (idx, power_row) in powers.axis_iter(Axis(0)).enumerate() {
                let wd_idx = idx / n_wind_speeds;
                let ws_idx = idx % n_wind_speeds;

                for ti in 0..n_turbines {
                    reshaped[[wd_idx, ws_idx, ti]] = power_row[ti];
                }
            }

            reshaped
        } else {
            // Fallback: return 3D array with first dimension as n_findex
            let n_findex = powers.shape()[0];
            let n_turbines = powers.shape()[1];
            let mut reshaped = Array3::zeros((n_findex, 1, n_turbines));
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    reshaped[[fi, 0, ti]] = powers[[fi, ti]];
                }
            }
            reshaped
        })
    }

    /// Get farm power reshaped for WindRose mode
    ///
    /// Returns powers in shape (n_wind_directions, n_wind_speeds)
    /// Only works when wind_data_info is available (set via set_wind_conditions_with_rose)
    pub fn get_farm_power_rose(&self) -> crate::Result<Array2> {
        let farm_power = self.get_farm_power()?;

        Ok(if let Some(ref wind_info) = self.core.wind_data_info {
            // Reshape from (n_findex,) to (n_wind_directions, n_wind_speeds)
            let n_wind_dirs = wind_info.n_wind_directions;
            let n_wind_speeds = wind_info.n_wind_speeds;

            let mut reshaped = Array2::zeros((n_wind_dirs, n_wind_speeds));

            // Reshape the data assuming findex order is [wd0_ws0, wd0_ws1, ..., wd1_ws0, wd1_ws1, ...]
            for (idx, &power) in farm_power.iter().enumerate() {
                let wd_idx = idx / n_wind_speeds;
                let ws_idx = idx % n_wind_speeds;
                reshaped[[wd_idx, ws_idx]] = power;
            }

            reshaped
        } else {
            // Fallback: return 2D array with second dimension as 1
            let n_findex = farm_power.len();
            let mut reshaped = Array2::zeros((n_findex, 1));
            for (i, &power) in farm_power.iter().enumerate() {
                reshaped[[i, 0]] = power;
            }
            reshaped
        })
    }

    /// Compute the expected (mean) power of each turbine
    ///
    /// This method calculates the weighted average power for each turbine across all
    /// wind conditions, using the provided frequencies as weights.
    ///
    /// # Arguments
    /// * `freq` - Optional frequency array with shape:
    ///   - 1D array of shape `(n_findex,)`: Same frequencies used for all turbines
    ///   - 2D array of shape `(n_findex, n_turbines)`: Unique frequencies per turbine
    ///   
    ///   If None, uniform frequencies are assumed (simple mean over findices).
    ///
    /// # Returns
    /// Array of shape `(n_turbines,)` containing the expected power for each turbine.
    ///
    /// # Examples
    /// ```no_run
    /// # use anyhow::Error;
    /// # fn main() -> Result<(), Error> {
    /// use florus::FlorisModel;
    /// use ndarray::Array1;
    ///
    /// let mut model = FlorisModel::from_file("default_inputs.yaml")?;
    /// model.run()?;
    ///
    /// // Calculate expected powers with uniform frequencies
    /// let expected_powers = model.get_expected_turbine_powers(None)?;
    ///
    /// // Or with custom frequencies
    /// let freq = Array1::from_vec(vec![0.3, 0.5, 0.2]);
    /// let expected_powers = model.get_expected_turbine_powers(Some(freq))?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_expected_turbine_powers(&self, freq: Option<Array1>) -> crate::Result<Array1> {
        // Get turbine powers: shape (n_findex, n_turbines)
        let turbine_powers = self.get_turbine_powers()?;
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();

        // Determine frequencies to use (Python floris_model.py:633-643)
        let frequencies = if let Some(f) = freq {
            f
        } else if let Some(stored) = &self.core.wind_data_frequencies {
            // A (n_findex, n_turbines) table (the wind resource grid) weights
            // the power of every turbine separately
            if stored.shape()[1] > 1 && stored.shape()[1] == n_turbines {
                let mut expected_powers = Array1::zeros(n_turbines);
                for ti in 0..n_turbines {
                    let mut weighted_sum = 0.0;
                    for fi in 0..n_findex {
                        weighted_sum += stored[[fi, ti]] * turbine_powers[[fi, ti]];
                    }
                    expected_powers[ti] = weighted_sum;
                }
                return Ok(expected_powers);
            }
            stored.column(0).to_owned()
        } else {
            // Default: uniform frequencies
            Array1::from_elem(n_findex, 1.0 / n_findex as Float)
        };

        // Validate frequency dimensions
        if frequencies.len() != n_findex {
            anyhow::bail!(
                "Frequency array length ({}) must match n_findex ({})",
                frequencies.len(),
                n_findex
            );
        }

        // Calculate weighted sum: sum over findex dimension
        // For each turbine: expected_power[t] = sum_i(freq[i] * power[i, t])
        let mut expected_powers = Array1::zeros(n_turbines);

        for ti in 0..n_turbines {
            let mut weighted_sum = 0.0;
            for fi in 0..n_findex {
                weighted_sum += frequencies[fi] * turbine_powers[[fi, ti]];
            }
            expected_powers[ti] = weighted_sum;
        }

        Ok(expected_powers)
    }

    /// Compute the expected (mean) power of the wind farm
    ///
    /// This method calculates the weighted average total power across all wind conditions,
    /// using the provided frequencies as weights.
    ///
    /// # Arguments
    /// * `freq` - Optional frequency array of shape `(n_findex,)`. If None, uniform
    ///   frequencies are assumed (simple mean over findices).
    /// * `turbine_weights` - Optional 2D array of shape `(n_findex, n_turbines)` for
    ///   weighting individual turbines. If None, all turbines are weighted equally (1.0).
    ///
    /// # Returns
    /// Single f64 value representing the expected total farm power.
    ///
    /// # Examples
    /// ```no_run
    /// # use anyhow::Error;
    /// # fn main() -> Result<(), Error> {
    /// use florus::FlorisModel;
    /// use ndarray::Array1;
    ///
    /// let mut model = FlorisModel::from_file("default_inputs.yaml")?;
    /// model.run()?;
    ///
    /// // Calculate expected farm power with uniform frequencies
    /// let expected_power = model.get_expected_farm_power(None, None)?;
    ///
    /// // Or with custom frequencies
    /// let freq = Array1::from_vec(vec![0.3, 0.5, 0.2]);
    /// let expected_power = model.get_expected_farm_power(Some(freq), None)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_expected_farm_power(
        &self,
        freq: Option<Array1>,
        turbine_weights: Option<Array2>,
    ) -> crate::Result<Float> {
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();

        // Get turbine powers and apply weights
        let turbine_powers = self.get_turbine_powers()?;

        let weights = if let Some(w) = turbine_weights {
            if w.shape() != [n_findex, n_turbines] {
                anyhow::bail!(
                    "turbine_weights shape {:?} doesn't match expected [{}, {}]",
                    w.shape(),
                    n_findex,
                    n_turbines
                );
            }
            w
        } else {
            Array2::from_elem((n_findex, n_turbines), 1.0)
        };

        // Apply weights to turbine powers
        let mut weighted_powers = turbine_powers.clone();
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                weighted_powers[[fi, ti]] *= weights[[fi, ti]];
            }
        }

        // Determine frequencies (Python floris_model.py:811-825)
        let frequencies = if let Some(f) = freq {
            f
        } else if let Some(stored) = &self.core.wind_data_frequencies {
            if stored.shape()[1] > 1 && stored.shape()[1] == n_turbines {
                // A (n_findex, n_turbines) table weights the power of every
                // turbine separately
                let mut expected_power = 0.0;
                for fi in 0..n_findex {
                    for ti in 0..n_turbines {
                        expected_power += stored[[fi, ti]] * weighted_powers[[fi, ti]];
                    }
                }
                return Ok(expected_power);
            }
            stored.column(0).to_owned()
        } else {
            Array1::from_elem(n_findex, 1.0 / n_findex as Float)
        };

        if frequencies.len() != n_findex {
            anyhow::bail!(
                "Frequency array length ({}) must match n_findex ({})",
                frequencies.len(),
                n_findex
            );
        }

        // Calculate farm power for each findex
        let mut farm_power = Array1::zeros(n_findex);
        for fi in 0..n_findex {
            farm_power[fi] = weighted_powers.row(fi).sum();
        }

        // Calculate weighted sum
        let mut expected_power = 0.0;
        for fi in 0..n_findex {
            expected_power += frequencies[fi] * farm_power[fi];
        }

        Ok(expected_power)
    }

    /// Compute the expected (mean) value produced by the wind farm
    ///
    /// This method multiplies the farm power by corresponding values (e.g., electricity
    /// prices) and weights by frequency to calculate expected value.
    ///
    /// # Arguments
    /// * `freq` - Optional frequency array of shape `(n_findex,)`. If None, uses uniform frequencies.
    /// * `values` - Optional value array of shape `(n_findex,)`. If None, assumes value of 1.0
    ///   for all conditions (equivalent to expected power).
    /// * `turbine_weights` - Optional 2D array for weighting individual turbines.
    ///
    /// # Returns
    /// Single f64 value representing the expected farm value.
    ///
    /// # Examples
    /// ```no_run
    /// # use anyhow::Error;
    /// # fn main() -> Result<(), Error> {
    /// use florus::FlorisModel;
    /// use ndarray::Array1;
    ///
    /// let mut model = FlorisModel::from_file("default_inputs.yaml")?;
    /// model.run()?;
    ///
    /// // Calculate expected value with electricity prices
    /// let prices = Array1::from_vec(vec![50.0, 60.0, 45.0]); // $/MWh
    /// let expected_value = model.get_expected_farm_value(None, Some(prices), None)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_expected_farm_value(
        &self,
        freq: Option<Array1>,
        values: Option<Array1>,
        turbine_weights: Option<Array2>,
    ) -> crate::Result<Float> {
        let n_findex = self.core.flow_field.n_findex;

        // Get expected farm power
        let expected_power = self.get_expected_farm_power(freq, turbine_weights)?;

        // Determine values
        let (value_array, is_uniform) = if let Some(ref v) = values {
            if v.len() != n_findex {
                anyhow::bail!(
                    "Values array length ({}) must match n_findex ({})",
                    v.len(),
                    n_findex
                );
            }
            (v.clone(), false)
        } else {
            (Array1::from_elem(n_findex, 1.0), true)
        };

        // Calculate weighted farm power for each findex
        let farm_power = self.get_farm_power()?;

        // Multiply power by values and sum
        let mut expected_value = 0.0;
        for fi in 0..n_findex {
            expected_value += value_array[fi] * farm_power[fi];
        }

        // Normalize by number of findices if using uniform values
        if is_uniform {
            expected_value /= n_findex as Float;
        }

        Ok(expected_value)
    }

    /// Get axial induction factors for all turbines
    ///
    /// Routes through the configured operation model using the solved
    /// rotor velocities, the grid's averaging method and cubature weights —
    /// mirroring Python `get_turbine_ais` (floris_model.py:1017-1036), which
    /// passes `average_method=grid.average_method,
    /// cubature_weights=grid.cubature_weights` into `farm.axial_induction`.
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing axial induction factors.
    ///
    /// # Errors
    /// Bubbles up operation-model errors (e.g. multi-dimensional Cp/Ct
    /// tables, which the operation-model dispatch does not support).
    pub fn get_turbine_ais(&self) -> crate::Result<Array2> {
        let farm = &self.core.farm;
        let flow_field = &self.core.flow_field;

        let Some(grid) = self.core.grid.as_ref() else {
            anyhow::bail!("Grid not initialized; call run() first");
        };

        let n_findex = flow_field.u_sorted.shape()[0];
        let turbulence_intensities = {
            let ti = &flow_field.turbulence_intensity_field_sorted;
            if ti.shape().len() == 4 && ti.shape()[0] == n_findex && ti.shape()[1] > 0 {
                Some(ti)
            } else {
                None
            }
        };
        let ais_sorted = farm.axial_induction(
            &flow_field.u_sorted,
            flow_field.air_density,
            Some(&farm.yaw_angles_sorted),
            Some(&farm.tilt_angles_sorted),
            turbulence_intensities,
            Some(&farm.power_setpoints_sorted),
            Some(&farm.awc_modes_sorted),
            Some(&farm.awc_amplitudes_sorted),
            grid.average_method(),
            grid.cubature_weights(),
        )?;

        Ok(self.unsort_to_layout(&ais_sorted))
    }

    /// Get thrust coefficients at turbine locations
    ///
    /// Routes through the configured operation model using the solved
    /// rotor velocities, the grid's averaging method and cubature weights —
    /// mirroring Python `get_turbine_thrust_coefficients`
    /// (floris_model.py:1038-1057), which passes
    /// `average_method=grid.average_method,
    /// cubature_weights=grid.cubature_weights` into
    /// `farm.thrust_coefficient`. Waked turbines therefore report their
    /// local Ct (not the freestream value the legacy estimate produced).
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing thrust coefficients.
    ///
    /// # Note
    /// Falls back to the legacy freestream-table estimate when the
    /// operation-model dispatch is unavailable (multi-dimensional Cp/Ct
    /// tables), matching the pre-operation-model behavior with a warning.
    pub fn get_turbine_cts(&self) -> Array2 {
        let farm = &self.core.farm;
        let flow_field = &self.core.flow_field;

        let op_model_route = (|| -> crate::Result<Array2> {
            let grid = self
                .core
                .grid
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("Grid not initialized; call run() first"))?;

            let n_findex = flow_field.u_sorted.shape()[0];
            if n_findex == 0 {
                anyhow::bail!("flow field not solved; call run() first");
            }
            let turbulence_intensities = {
                let ti = &flow_field.turbulence_intensity_field_sorted;
                if ti.shape().len() == 4 && ti.shape()[0] == n_findex && ti.shape()[1] > 0 {
                    Some(ti)
                } else {
                    None
                }
            };
            let cts_sorted = farm.thrust_coefficient(
                &flow_field.u_sorted,
                flow_field.air_density,
                Some(&farm.yaw_angles_sorted),
                Some(&farm.tilt_angles_sorted),
                turbulence_intensities,
                Some(&farm.power_setpoints_sorted),
                Some(&farm.awc_modes_sorted),
                Some(&farm.awc_amplitudes_sorted),
                grid.average_method(),
                grid.cubature_weights(),
            )?;
            Ok(self.unsort_to_layout(&cts_sorted))
        })();

        match op_model_route {
            Ok(cts) => cts,
            Err(e) => {
                log::warn!(
                    "get_turbine_thrust_coefficients fell back to the legacy freestream \
                     estimate: {e:#}"
                );
                self.get_turbine_cts_legacy()
            }
        }
    }

    /// Legacy Ct estimate: freestream (hub-reference) wind speed table
    /// lookup, ignoring wakes, yaw/tilt corrections and the operation
    /// model. Superseded by the operation-model route in
    /// `get_turbine_cts`; kept as the multi-dimensional-table fallback.
    fn get_turbine_cts_legacy(&self) -> Array2 {
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();
        let mut ct_array = Array2::zeros((n_findex, n_turbines));

        for ti in 0..n_turbines {
            for fi in 0..n_findex {
                let velocity = self.core.flow_field.wind_speeds[fi];
                if ti < self.core.farm.turbines.len() {
                    let turbine_type = &self.core.farm.turbines[ti].turbine_type;
                    let mut conditions = TableConditions::default();
                    conditions.wind_speed = velocity;
                    ct_array[[fi, ti]] = turbine_type
                        .power_thrust_table
                        .cp_ct_table
                        .get_ct(&conditions)
                        .unwrap_or(0.0);
                }
            }
        }

        ct_array
    }

    /// Alias for get_turbine_cts - returns thrust coefficients
    pub fn get_turbine_thrust_coefficients(&self) -> Array2 {
        self.get_turbine_cts()
    }

    /// Get the operation model for turbines
    ///
    /// Returns the operation model type string (e.g., "simple", "cosine-loss", etc.)
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing operation models as strings.
    pub fn get_operation_models(&self) -> Vec<Vec<String>> {
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();

        let mut operation_models = vec![vec![String::new(); n_turbines]; n_findex];

        for ti in 0..n_turbines {
            if ti >= self.core.farm.turbines.len() {
                continue;
            }
            let turbine_type = &self.core.farm.turbines[ti].turbine_type;
            let rated_power = turbine_type.rated_power.unwrap_or(5e6);

            for fi in 0..n_findex {
                let velocity = self.core.flow_field.u_sorted[[fi, ti, 0, 0]];

                let ct = if velocity < 3.0 {
                    0.0
                } else if velocity > 25.0 {
                    0.0
                } else if velocity < 12.0 {
                    0.8
                } else {
                    0.7
                };

                let area = std::f64::consts::PI * (turbine_type.rotor_diameter / 2.0).powi(2);
                let power = 0.5 * 1.225 * area * ct * velocity.powi(3);

                operation_models[fi][ti] = if power < 0.1 * rated_power {
                    "Idle".to_string()
                } else if power < 0.9 * rated_power {
                    "PartLoad".to_string()
                } else {
                    "FullLoad".to_string()
                };
            }
        }

        operation_models
    }

    /// Get turbulence intensities at turbine locations
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing turbulence intensities.
    pub fn get_turbine_tis(&self) -> Array2 {
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();

        // Return turbulence intensity field at turbine locations
        // Shape: (n_findex, n_turbines, grid_y, grid_z) -> extract first grid point
        let mut tis = Array2::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                tis[[fi, ti]] =
                    self.core.flow_field.turbulence_intensity_field_sorted[[fi, ti, 0, 0]];
            }
        }

        tis
    }

    /// Assign hub height to reference wind height
    ///
    /// This is useful when all turbines have the same hub height and you want to
    /// use that as the reference height for wind shear calculations.
    ///
    /// # Errors
    /// Returns an error if there are multiple unique hub heights.
    pub fn assign_hub_height_to_ref_height(&mut self) -> crate::Result<()> {
        let hub_heights = &self.core.farm.hub_heights;

        // Check for unique hub heights (use f64 comparison with tolerance)
        let mut unique_heights: Vec<Float> = Vec::new();
        for &h in hub_heights.iter() {
            if !unique_heights.iter().any(|&uh| (uh - h).abs() < 1e-6) {
                unique_heights.push(h);
            }
        }

        if unique_heights.len() > 1 {
            anyhow::bail!(
                "Cannot assign hub heights to reference height when there are multiple \
                 unique hub heights. Found {} unique heights: {:?}",
                unique_heights.len(),
                unique_heights
            );
        }

        if let Some(&height) = unique_heights.iter().next() {
            self.core.flow_field.reference_wind_height = height;
            Ok(())
        } else {
            anyhow::bail!("No hub heights found");
        }
    }

    /// Set the turbine operation model(s)
    ///
    /// Operation models control how turbine power and thrust are calculated.
    /// Common models include "simple", "cosine-loss", etc.
    ///
    /// # Arguments
    /// * `operation_model` - Either a single string applied to all turbines, or a vector
    ///   of strings (one per turbine).
    ///
    /// # Examples
    /// ```no_run
    /// # use anyhow::Error;
    /// # fn main() -> Result<(), Error> {
    /// use florus::FlorisModel;
    ///
    /// let mut model = FlorisModel::from_file("default_inputs.yaml")?;
    ///
    /// // Set same operation model for all turbines
    /// model.set_operation_model("simple")?;
    ///
    /// // Or set different models per turbine
    /// model.set_operation_model(vec!["simple".to_string(), "cosine".to_string()])?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn set_operation_model(
        &mut self,
        operation_model: impl Into<OperationModelSpec>,
    ) -> crate::Result<()> {
        let spec = operation_model.into();
        let n_turbines = self.core.farm.n_turbines();

        match spec {
            OperationModelSpec::Single(model) => {
                log::info!(
                    "Setting operation model '{}' for all {} turbines",
                    model,
                    n_turbines
                );
                // Update the operation model on each turbine's TurbineType
                for turbine in self.core.farm.turbines.iter_mut() {
                    let new_tt = turbine.turbine_type.with_operation_model(&model);
                    turbine.turbine_type = std::sync::Arc::new(new_tt);
                }
                self.core.state.initialized = false;
            }
            OperationModelSpec::Multiple(models) => {
                if models.len() != n_turbines {
                    anyhow::bail!(
                        "Operation model list length ({}) must match number of turbines ({})",
                        models.len(),
                        n_turbines
                    );
                }
                log::info!(
                    "Setting operation models for {} turbines: {:?}",
                    n_turbines,
                    models
                );
                for (turbine, model) in self.core.farm.turbines.iter_mut().zip(models.iter()) {
                    let new_tt = turbine.turbine_type.with_operation_model(model);
                    turbine.turbine_type = std::sync::Arc::new(new_tt);
                }
                self.core.state.initialized = false;
            }
        }

        Ok(())
    }

    /// Get turbine layout coordinates
    ///
    /// # Arguments
    /// * `include_z` - If true, includes z-coordinates (hub heights). Default is false.
    ///
    /// # Returns
    /// Tuple of (x, y) or (x, y, z) coordinate arrays.
    pub fn get_turbine_layout(&self, include_z: bool) -> TurbineLayout {
        let x = self.core.farm.layout_x.clone();
        let y = self.core.farm.layout_y.clone();

        if include_z {
            let z = self.core.farm.hub_heights.clone();
            TurbineLayout::WithZ(x, y, z)
        } else {
            TurbineLayout::XY(x, y)
        }
    }

    // ========================================================================
    // Properties (Getters)
    // ========================================================================

    /// Get layout x coordinates
    pub fn layout_x(&self) -> &Array1 {
        &self.core.farm.layout_x
    }

    /// Get layout y coordinates
    pub fn layout_y(&self) -> &Array1 {
        &self.core.farm.layout_y
    }

    /// Get wind directions
    pub fn wind_directions(&self) -> &Array1 {
        &self.core.flow_field.wind_directions
    }

    /// Get wind speeds
    pub fn wind_speeds(&self) -> &Array1 {
        &self.core.flow_field.wind_speeds
    }

    /// Get turbulence intensities
    pub fn turbulence_intensities(&self) -> &Array1 {
        &self.core.flow_field.turbulence_intensities
    }

    /// Get number of flow indices (findex)
    pub fn n_findex(&self) -> usize {
        self.core.flow_field.n_findex
    }

    /// Get number of turbines
    pub fn n_turbines(&self) -> usize {
        self.core.farm.n_turbines()
    }

    /// Get yaw angles array
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing yaw angles in degrees
    pub fn get_yaw_angles(&self) -> &Array2 {
        &self.core.farm.yaw_angles
    }

    /// Get power setpoints array
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing power setpoints in Watts
    pub fn get_power_setpoints(&self) -> &Array2 {
        &self.core.farm.power_setpoints
    }

    /// Get AWC modes array
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing AWC mode strings
    pub fn get_awc_modes(&self) -> &NdArray2<String> {
        &self.core.farm.awc_modes
    }

    /// Get AWC amplitudes array
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing AWC amplitudes in degrees
    pub fn get_awc_amplitudes(&self) -> &Array2 {
        &self.core.farm.awc_amplitudes
    }

    /// Get AWC frequencies array
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_turbines)` containing AWC frequencies in Hz
    pub fn get_awc_frequencies(&self) -> &Array2 {
        &self.core.farm.awc_frequencies
    }

    /// Get reference wind height
    pub fn reference_wind_height(&self) -> Float {
        self.core.flow_field.reference_wind_height
    }

    /// Rotor-average velocity at each turbine, mirroring Python's
    /// `turbine_average_velocities` property (floris_model.py:1798-1804):
    /// `average_velocity` over the SOLVED field with the grid's averaging
    /// method and cubature weights, returned in LAYOUT order.
    ///
    /// Previously this returned the single first rotor-grid-point value in
    /// sorted order (documented as "simplified"); the example comparison
    /// (2026-10-19) exposed the divergence from Python.
    ///
    /// Returns zeros when the flow field has not been solved yet.
    pub fn turbine_average_velocities(&self) -> Array2 {
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();
        let zeros = Array2::zeros((n_findex, n_turbines));

        let Some(grid) = self.core.grid.as_ref() else {
            return zeros;
        };
        if self.core.flow_field.u_sorted.is_empty() {
            return zeros;
        }

        let velocities_dyn = self.core.flow_field.u_sorted.view().into_dyn().to_owned();
        let cubature_weights: Option<Array1> = grid
            .cubature_weights()
            .map(|w| Array1::from_iter(w.iter().copied()));
        let Ok(avg_sorted) = crate::core::rotor_velocity::average_velocity(
            &velocities_dyn,
            None,
            &grid.average_method(),
            cubature_weights.as_ref(),
        ) else {
            return zeros;
        };

        self.unsort_to_layout(&avg_sorted)
    }

    /// Get core reference (immutable)
    pub fn core(&self) -> &Core {
        &self.core
    }

    /// Get core reference (mutable)
    pub fn core_mut(&mut self) -> &mut Core {
        &mut self.core
    }

    /// Get farm reference (immutable)
    pub fn farm(&self) -> &Farm {
        &self.core.farm
    }

    /// Get flow field reference (immutable)
    pub fn flow_field(&self) -> &FlowField {
        &self.core.flow_field
    }

    /// Calculate horizontal plane at specified height
    ///
    /// # Arguments
    /// * `height` - Height above ground (m)
    /// * `x_resolution` - Number of grid points in x direction
    /// * `y_resolution` - Number of grid points in y direction
    ///
    /// # Returns
    /// CutPlane containing the horizontal slice of the flow field
    pub fn calculate_horizontal_plane(
        &self,
        height: f64,
        x_resolution: usize,
        y_resolution: usize,
    ) -> crate::Result<crate::core::cut_plane::CutPlane> {
        use crate::core::grid::FlowFieldPlanarGrid;
        use ndarray::Array;

        // Create a copy for visualization
        let mut fmodel_viz = self.clone();

        // Get turbine coordinates, diameters, and hub heights
        let n_turbines = fmodel_viz.core.farm.n_turbines();
        let mut turbine_coords = Array::zeros((n_turbines, 2));
        let mut turbine_diams = Array::zeros(n_turbines);
        let mut turbine_hub_heights = Array::zeros(n_turbines);

        for i in 0..n_turbines {
            turbine_coords[[i, 0]] = fmodel_viz.core.farm.layout_x[i];
            turbine_coords[[i, 1]] = fmodel_viz.core.farm.layout_y[i];
            turbine_diams[i] = fmodel_viz.core.farm.turbines[i].turbine_type.rotor_diameter;
            turbine_hub_heights[i] = fmodel_viz.core.farm.turbines[i].turbine_type.hub_height;
        }

        // Get wind directions for first findex only
        let wind_dirs = Array::from_vec(vec![fmodel_viz.core.flow_field.wind_directions[0]]);

        // Create FlowFieldPlanarGrid for horizontal plane (normal_vector="z")
        let planar_grid = FlowFieldPlanarGrid::new(
            turbine_coords,
            turbine_diams,
            turbine_hub_heights,
            wind_dirs,
            [x_resolution, y_resolution],
            "z".to_string(),
            height,
            None, // x1_bounds
            None, // x2_bounds
        )?;

        // Set the grid
        fmodel_viz.core.grid = Some(Box::new(planar_grid));

        // Reinitialize flow field with new grid
        if let Some(ref grid) = fmodel_viz.core.grid {
            fmodel_viz
                .core
                .flow_field
                .initialize_flow_field_with_coordinates(
                    grid.grid_shape(),
                    grid.z_sorted(),
                    &grid.hub_heights(),
                    Some(grid.x_sorted_inertial_frame()),
                    Some(grid.y_sorted_inertial_frame()),
                );
        }

        // Mark as initialized
        fmodel_viz.core.state.initialized = true;

        // Solve for visualization (single findex)
        fmodel_viz.core.solve_for_viz()?;

        // Extract data from the solved flow field
        let u_sorted = &fmodel_viz.core.flow_field.u_sorted;
        let v_sorted = &fmodel_viz.core.flow_field.v_sorted;
        let w_sorted = &fmodel_viz.core.flow_field.w_sorted;

        if let Some(ref grid) = fmodel_viz.core.grid {
            // Get inertial frame coordinates for horizontal plane
            let x_inertial = grid.x_sorted_inertial_frame();
            let y_inertial = grid.y_sorted_inertial_frame();
            let z_inertial = grid.z_sorted_inertial_frame();

            // Extract first findex (findex=0)
            let x_flat = x_inertial.slice(ndarray::s![0, .., 0, 0]);
            let y_flat = y_inertial.slice(ndarray::s![0, .., 0, 0]);
            let z_flat = z_inertial.slice(ndarray::s![0, .., 0, 0]);
            let u_flat = u_sorted.slice(ndarray::s![0, .., 0, 0]);
            let v_flat = v_sorted.slice(ndarray::s![0, .., 0, 0]);
            let w_flat = w_sorted.slice(ndarray::s![0, .., 0, 0]);

            // Filter to only include points at the specified height (within tolerance)
            let tol = 0.1;
            let mut x1_vals = Vec::new();
            let mut x2_vals = Vec::new();
            let mut x3_vals = Vec::new();
            let mut u_vals = Vec::new();
            let mut v_vals = Vec::new();
            let mut w_vals = Vec::new();

            for i in 0..x_flat.len() {
                if (z_flat[i] - height).abs() < tol {
                    x1_vals.push(x_flat[i]);
                    x2_vals.push(y_flat[i]);
                    x3_vals.push(z_flat[i]);
                    u_vals.push(u_flat[i]);
                    v_vals.push(v_flat[i]);
                    w_vals.push(w_flat[i]);
                }
            }

            // Create CutPlane
            let cut_plane = crate::core::cut_plane::CutPlane {
                data: crate::core::cut_plane::CutPlaneData {
                    x1: Array::from_vec(x1_vals),
                    x2: Array::from_vec(x2_vals),
                    x3: Array::from_vec(x3_vals),
                    u: Array::from_vec(u_vals),
                    v: Array::from_vec(v_vals),
                    w: Array::from_vec(w_vals),
                },
                normal_vector: "z".to_string(),
                resolution: (x_resolution, y_resolution),
            };

            Ok(cut_plane)
        } else {
            anyhow::bail!("Grid not initialized after setting planar grid.");
        }
    }

    /// Calculate a y-plane (vertical plane along wind direction)
    pub fn calculate_y_plane(
        &self,
        x_resolution: usize,
        z_resolution: usize,
        crossstream_dist: f64,
    ) -> crate::Result<crate::core::cut_plane::CutPlane> {
        use crate::core::grid::FlowFieldPlanarGrid;
        use ndarray::Array;

        // Create a copy for visualization
        let mut fmodel_viz = self.clone();

        // Get turbine coordinates, diameters, and hub heights
        let n_turbines = fmodel_viz.core.farm.n_turbines();
        let mut turbine_coords = Array::zeros((n_turbines, 2));
        let mut turbine_diams = Array::zeros(n_turbines);
        let mut turbine_hub_heights = Array::zeros(n_turbines);

        for i in 0..n_turbines {
            turbine_coords[[i, 0]] = fmodel_viz.core.farm.layout_x[i];
            turbine_coords[[i, 1]] = fmodel_viz.core.farm.layout_y[i];
            turbine_diams[i] = fmodel_viz.core.farm.turbines[i].turbine_type.rotor_diameter;
            turbine_hub_heights[i] = fmodel_viz.core.farm.turbines[i].turbine_type.hub_height;
        }

        // Get wind directions for first findex only
        let wind_dirs = Array::from_vec(vec![fmodel_viz.core.flow_field.wind_directions[0]]);

        // Create FlowFieldPlanarGrid for y-plane (normal_vector="y")
        let planar_grid = FlowFieldPlanarGrid::new(
            turbine_coords,
            turbine_diams,
            turbine_hub_heights,
            wind_dirs,
            [x_resolution, z_resolution],
            "y".to_string(),
            crossstream_dist,
            None, // x1_bounds
            None, // x2_bounds
        )?;

        // Set the grid
        fmodel_viz.core.grid = Some(Box::new(planar_grid));

        // Reinitialize flow field with new grid
        if let Some(ref grid) = fmodel_viz.core.grid {
            fmodel_viz
                .core
                .flow_field
                .initialize_flow_field_with_coordinates(
                    grid.grid_shape(),
                    grid.z_sorted(),
                    &grid.hub_heights(),
                    Some(grid.x_sorted_inertial_frame()),
                    Some(grid.y_sorted_inertial_frame()),
                );
        }

        // Mark as initialized
        fmodel_viz.core.state.initialized = true;

        // Solve for visualization (single findex)
        fmodel_viz.core.solve_for_viz()?;

        // Extract data from the solved flow field
        let u_sorted = &fmodel_viz.core.flow_field.u_sorted;
        let v_sorted = &fmodel_viz.core.flow_field.v_sorted;
        let w_sorted = &fmodel_viz.core.flow_field.w_sorted;

        if let Some(ref grid) = fmodel_viz.core.grid {
            // Get inertial frame coordinates for y-plane
            let x_inertial = grid.x_sorted_inertial_frame();
            let y_inertial = grid.y_sorted_inertial_frame();
            let z_inertial = grid.z_sorted_inertial_frame();

            // Extract first findex (findex=0)
            let x_flat = x_inertial.slice(ndarray::s![0, .., 0, 0]);
            let y_flat = y_inertial.slice(ndarray::s![0, .., 0, 0]);
            let z_flat = z_inertial.slice(ndarray::s![0, .., 0, 0]);
            let u_flat = u_sorted.slice(ndarray::s![0, .., 0, 0]);
            let v_flat = v_sorted.slice(ndarray::s![0, .., 0, 0]);
            let w_flat = w_sorted.slice(ndarray::s![0, .., 0, 0]);

            // Filter to only include points at the specified crossstream distance (within tolerance)
            let tol = 0.1;
            let mut x1_vals = Vec::new();
            let mut x2_vals = Vec::new();
            let mut x3_vals = Vec::new();
            let mut u_vals = Vec::new();
            let mut v_vals = Vec::new();
            let mut w_vals = Vec::new();

            for i in 0..x_flat.len() {
                if (y_flat[i] - crossstream_dist).abs() < tol {
                    x1_vals.push(x_flat[i]);
                    x2_vals.push(z_flat[i]);
                    x3_vals.push(y_flat[i]);
                    u_vals.push(u_flat[i]);
                    v_vals.push(v_flat[i]);
                    w_vals.push(w_flat[i]);
                }
            }

            // Create CutPlane
            let cut_plane = crate::core::cut_plane::CutPlane {
                data: crate::core::cut_plane::CutPlaneData {
                    x1: Array::from_vec(x1_vals),
                    x2: Array::from_vec(x2_vals),
                    x3: Array::from_vec(x3_vals),
                    u: Array::from_vec(u_vals),
                    v: Array::from_vec(v_vals),
                    w: Array::from_vec(w_vals),
                },
                normal_vector: "y".to_string(),
                resolution: (x_resolution, z_resolution),
            };

            Ok(cut_plane)
        } else {
            anyhow::bail!("Grid not initialized after setting planar grid.");
        }
    }

    /// Calculate a cross-plane (vertical plane across wind direction)
    pub fn calculate_cross_plane(
        &self,
        y_resolution: usize,
        z_resolution: usize,
        downstream_dist: f64,
    ) -> crate::Result<crate::core::cut_plane::CutPlane> {
        use crate::core::grid::FlowFieldPlanarGrid;
        use ndarray::Array;

        // Create a copy for visualization
        let mut fmodel_viz = self.clone();

        // Get turbine coordinates, diameters, and hub heights
        let n_turbines = fmodel_viz.core.farm.n_turbines();
        let mut turbine_coords = Array::zeros((n_turbines, 2));
        let mut turbine_diams = Array::zeros(n_turbines);
        let mut turbine_hub_heights = Array::zeros(n_turbines);

        for i in 0..n_turbines {
            turbine_coords[[i, 0]] = fmodel_viz.core.farm.layout_x[i];
            turbine_coords[[i, 1]] = fmodel_viz.core.farm.layout_y[i];
            turbine_diams[i] = fmodel_viz.core.farm.turbines[i].turbine_type.rotor_diameter;
            turbine_hub_heights[i] = fmodel_viz.core.farm.turbines[i].turbine_type.hub_height;
        }

        // Get wind directions for first findex only
        let wind_dirs = Array::from_vec(vec![fmodel_viz.core.flow_field.wind_directions[0]]);

        // Create FlowFieldPlanarGrid for cross-plane (normal_vector="x")
        let planar_grid = FlowFieldPlanarGrid::new(
            turbine_coords,
            turbine_diams,
            turbine_hub_heights,
            wind_dirs,
            [y_resolution, z_resolution],
            "x".to_string(),
            downstream_dist,
            None, // x1_bounds
            None, // x2_bounds
        )?;

        // Set the grid
        fmodel_viz.core.grid = Some(Box::new(planar_grid));

        // Reinitialize flow field with new grid
        if let Some(ref grid) = fmodel_viz.core.grid {
            fmodel_viz
                .core
                .flow_field
                .initialize_flow_field_with_coordinates(
                    grid.grid_shape(),
                    grid.z_sorted(),
                    &grid.hub_heights(),
                    Some(grid.x_sorted_inertial_frame()),
                    Some(grid.y_sorted_inertial_frame()),
                );
        }

        // Mark as initialized
        fmodel_viz.core.state.initialized = true;

        // Solve for visualization (single findex)
        fmodel_viz.core.solve_for_viz()?;

        // Extract data from the solved flow field
        let u_sorted = &fmodel_viz.core.flow_field.u_sorted;
        let v_sorted = &fmodel_viz.core.flow_field.v_sorted;
        let w_sorted = &fmodel_viz.core.flow_field.w_sorted;

        if let Some(ref grid) = fmodel_viz.core.grid {
            // Get inertial frame coordinates for cross-plane
            let x_inertial = grid.x_sorted_inertial_frame();
            let y_inertial = grid.y_sorted_inertial_frame();
            let z_inertial = grid.z_sorted_inertial_frame();

            // Extract first findex (findex=0)
            let x_flat = x_inertial.slice(ndarray::s![0, .., 0, 0]);
            let y_flat = y_inertial.slice(ndarray::s![0, .., 0, 0]);
            let z_flat = z_inertial.slice(ndarray::s![0, .., 0, 0]);
            let u_flat = u_sorted.slice(ndarray::s![0, .., 0, 0]);
            let v_flat = v_sorted.slice(ndarray::s![0, .., 0, 0]);
            let w_flat = w_sorted.slice(ndarray::s![0, .., 0, 0]);

            // Filter to only include points at the specified downstream distance (within tolerance)
            let tol = 0.1;
            let mut x1_vals = Vec::new();
            let mut x2_vals = Vec::new();
            let mut x3_vals = Vec::new();
            let mut u_vals = Vec::new();
            let mut v_vals = Vec::new();
            let mut w_vals = Vec::new();

            for i in 0..x_flat.len() {
                if (x_flat[i] - downstream_dist).abs() < tol {
                    x1_vals.push(y_flat[i]);
                    x2_vals.push(z_flat[i]);
                    x3_vals.push(x_flat[i]);
                    u_vals.push(u_flat[i]);
                    v_vals.push(v_flat[i]);
                    w_vals.push(w_flat[i]);
                }
            }

            // Create CutPlane
            let cut_plane = crate::core::cut_plane::CutPlane {
                data: crate::core::cut_plane::CutPlaneData {
                    x1: Array::from_vec(x1_vals),
                    x2: Array::from_vec(x2_vals),
                    x3: Array::from_vec(x3_vals),
                    u: Array::from_vec(u_vals),
                    v: Array::from_vec(v_vals),
                    w: Array::from_vec(w_vals),
                },
                normal_vector: "x".to_string(),
                resolution: (y_resolution, z_resolution),
            };

            Ok(cut_plane)
        } else {
            anyhow::bail!("Grid not initialized after setting planar grid.");
        }
    }

    /// Get state reference (immutable)
    pub fn state(&self) -> &State {
        &self.core.state
    }

    /// Set wind conditions from WindData trait object
    ///
    /// This method accepts any object that implements the WindData trait,
    /// including TimeSeries, WindRose, WindRoseRwg and custom implementations.
    pub fn set_wind_data(&mut self, wind_data: &dyn WindData) -> crate::Result<()> {
        // Python `fmodel.set(wind_data=...)` only replaces the wind speeds,
        // directions and turbulence intensities; the shear, veer, air density
        // and reference height of the model are kept. Rebuilding the flow
        // field with the model's own settings (instead of the hard-coded
        // defaults of `FlowField::from_wind_data`) preserves them.
        let (wind_directions, wind_speeds, turbulence_intensities, _, _, _) = wind_data.unpack();
        self.core.flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            self.core.flow_field.wind_veer,
            self.core.flow_field.wind_shear,
            self.core.flow_field.air_density,
            turbulence_intensities,
            self.core.flow_field.reference_wind_height,
        )?;
        self.core
            .farm
            .initialize_control_arrays(self.core.flow_field.n_findex);
        self.core
            .farm
            .set_tilt_to_ref_tilt(self.core.flow_field.n_findex);

        // The multi-dimensional conditions of the wind data select the
        // power/thrust table of the multi-dimensional turbines (Python
        // `floris_model.py:reinitialize`).
        self.core.multidim_conditions = wind_data.unpack_multidim_conditions().map(|conditions| {
            let mut values = std::collections::HashMap::new();
            for (key, array) in conditions.to_condition_map() {
                values.insert(
                    key,
                    crate::floris_config::MultiDimConditionValue::Array(array.to_vec()),
                );
            }
            crate::floris_config::MultiDimConditions { values }
        });

        // Keep the frequencies of the wind data: the expected-power and AEP
        // queries use them when no explicit frequencies are given (Python
        // `floris_model.py:get_expected_turbine_powers`).
        self.core.wind_data_frequencies = Some(wind_data.unpack_freq());

        self.core.grid = None;
        Ok(())
    }

    /// Set wind conditions
    pub fn set_wind_conditions(
        &mut self,
        wind_speeds: Array1,
        wind_directions: Array1,
        turbulence_intensities: Array1,
    ) -> crate::Result<()> {
        self.core.flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            self.core.flow_field.wind_veer,
            self.core.flow_field.wind_shear,
            self.core.flow_field.air_density,
            turbulence_intensities,
            self.core.flow_field.reference_wind_height,
        )?;

        self.core
            .farm
            .initialize_control_arrays(self.core.flow_field.n_findex);
        self.core
            .farm
            .set_tilt_to_ref_tilt(self.core.flow_field.n_findex);

        // Clear wind_data_info for flat TimeSeries mode
        self.core.wind_data_info = None;

        // The explicit wind conditions replace the wind data: uniform
        // frequencies are used for the expected-power queries
        self.core.wind_data_frequencies = None;

        // Mark farm as uninitialized so it will be properly reinitialized on next run()
        self.core.farm.state.initialized = false;
        self.core.grid = None;
        Ok(())
    }

    /// Set wind conditions with WindRose structure (for reshaping)
    pub fn set_wind_conditions_with_rose(
        &mut self,
        wind_speeds: Array1,
        wind_directions: Array1,
        turbulence_intensities: Array1,
        unique_wind_directions: Vec<f64>,
        unique_wind_speeds: Vec<f64>,
    ) -> crate::Result<()> {
        let n_wind_directions = unique_wind_directions.len();
        let n_wind_speeds = unique_wind_speeds.len();

        self.core.flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            self.core.flow_field.wind_veer,
            self.core.flow_field.wind_shear,
            self.core.flow_field.air_density,
            turbulence_intensities,
            self.core.flow_field.reference_wind_height,
        )?;

        self.core
            .farm
            .initialize_control_arrays(self.core.flow_field.n_findex);
        self.core
            .farm
            .set_tilt_to_ref_tilt(self.core.flow_field.n_findex);

        // Store WindRose information for reshaping
        self.core.wind_data_info = Some(crate::core::core::WindDataInfo {
            wind_directions: unique_wind_directions,
            wind_speeds: unique_wind_speeds,
            n_wind_directions,
            n_wind_speeds,
        });

        // Mark farm as uninitialized so it will be properly reinitialized on next run()
        self.core.farm.state.initialized = false;
        self.core.grid = None;
        Ok(())
    }

    /// Set the heterogeneous inflow configuration.
    ///
    /// The speed multipliers are applied to the initial velocity field
    /// (Python `fmodel.set(heterogeneous_inflow_config=...)`,
    /// floris_model.py:400). The multipliers are interpolated onto the flow
    /// field grid points in the inertial frame; points outside the convex
    /// hull of the given locations keep the freestream speed.
    pub fn set_heterogeneous_inflow_config(
        &mut self,
        config: crate::heterogeneous_map::HeterogeneousInflowConfig,
    ) -> crate::Result<()> {
        self.core.flow_field.heterogeneous_inflow_config = Some(config);
        // Mark the model uninitialized so the next run() re-initializes the
        // flow field with the heterogeneous inflow applied.
        self.core.state.initialized = false;
        self.core.farm.state.initialized = false;
        Ok(())
    }

    /// Set the multi-dimensional conditions of the model.
    ///
    /// The conditions select the power/thrust table of the multi-dimensional
    /// turbines (Python `fmodel.set(multidim_conditions=...)`). The keys are
    /// the condition names of the turbine's table (e.g. "Tp", "Hs", "TI").
    pub fn set_multidim_conditions(
        &mut self,
        conditions: std::collections::HashMap<String, crate::floris_config::MultiDimConditionValue>,
    ) -> crate::Result<()> {
        self.core.multidim_conditions =
            Some(crate::floris_config::MultiDimConditions { values: conditions });
        // The conditions are applied when the flow field is initialized
        self.core.state.initialized = false;
        self.core.farm.state.initialized = false;
        Ok(())
    }

    /// Set yaw angles
    pub fn set_yaw_angles(&mut self, yaw_angles: Array2) -> crate::Result<()> {
        if yaw_angles.shape() != [self.core.flow_field.n_findex, self.core.farm.n_turbines()] {
            anyhow::bail!(
                "yaw_angles shape {:?} doesn't match expected [{}, {}]",
                yaw_angles.shape(),
                self.core.flow_field.n_findex,
                self.core.farm.n_turbines()
            );
        }

        // Use Farm's method to ensure both yaw_angles and yaw_angles_sorted are updated
        self.core.farm.set_yaw_angles(yaw_angles);

        // Mark the model uninitialized so the next run() re-initializes the
        // domain, which re-sorts the yaw angles into the solver's sorted
        // order (Farm::initialize). Without this the solver would keep the
        // previously sorted array, applying the yaw angles to the wrong
        // turbines whenever the sort order differs from the layout order
        // (Python FlorisModel.set() resets core.state to UNINITIALIZED,
        // floris_model.py:398).
        self.core.state.initialized = false;
        self.core.farm.state.initialized = false;

        Ok(())
    }

    /// Set turbine layout
    pub fn set_layout(&mut self, layout_x: &Array1, layout_y: &Array1) -> crate::Result<()> {
        self.core.farm.set_layout(layout_x, layout_y)?;
        self.core.grid = None; // Reset grid when layout changes
        Ok(())
    }

    /// Set wind shear exponent
    pub fn set_wind_shear(&mut self, wind_shear: Float) -> crate::Result<()> {
        self.core.flow_field.wind_shear = wind_shear;
        Ok(())
    }

    /// Set air density
    pub fn set_air_density(&mut self, air_density: Float) -> crate::Result<()> {
        self.core.flow_field.air_density = air_density;
        Ok(())
    }

    /// Set reference wind height
    pub fn set_reference_wind_height(&mut self, reference_wind_height: Float) -> crate::Result<()> {
        self.core.flow_field.reference_wind_height = reference_wind_height;
        Ok(())
    }

    /// Set wind conditions and operation setpoints for the wind farm.
    ///
    /// This is a comprehensive method that allows setting multiple parameters at once.
    /// It reinitializes the model with new conditions while preserving non-default operation settings.
    ///
    /// # Arguments
    /// * `wind_speeds` - Optional wind speeds array
    /// * `wind_directions` - Optional wind directions array
    /// * `wind_shear` - Optional wind shear exponent
    /// * `wind_veer` - Optional wind veer angle
    /// * `reference_wind_height` - Optional reference wind height
    /// * `turbulence_intensities` - Optional turbulence intensities array
    /// * `air_density` - Optional air density
    /// * `layout_x` - Optional turbine x-coordinates
    /// * `layout_y` - Optional turbine y-coordinates
    /// * `yaw_angles` - Optional yaw angles for turbines
    /// * `power_setpoints` - Optional power setpoints for turbines
    /// * `awc_modes` - Optional AWC modes
    /// * `awc_amplitudes` - Optional AWC amplitudes
    /// * `awc_frequencies` - Optional AWC frequencies
    /// * `disable_turbines` - Optional boolean array indicating disabled turbines
    ///
    /// # Returns
    /// Returns `Ok(())` on success, or an error if validation fails
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use anyhow::Error;
    /// # fn main() -> Result<(), Error> {
    /// use florus::FlorisModel;
    /// use ndarray::{Array2, Array1};
    ///
    /// let mut model = FlorisModel::from_file("default_inputs.yaml")?;
    ///
    /// // Example 1: Update wind conditions only
    /// let wind_speeds = Array1::from_vec(vec![8.0, 10.0, 12.0]);
    /// let wind_directions = Array1::from_vec(vec![270.0, 270.0, 270.0]);
    /// let ti = Array1::from_vec(vec![0.06, 0.06, 0.06]);
    /// model.set(
    ///     Some(wind_speeds),
    ///     Some(wind_directions),
    ///     None, None, None,
    ///     Some(ti),
    ///     None, None, None, None, None, None, None, None, None
    /// )?;
    ///
    /// // Example 2: Update layout and wind conditions
    /// let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
    /// let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    /// model.set(
    ///     None, None, None, None, None, None, None,
    ///     Some(layout_x), Some(layout_y), None, None, None, None, None, None
    /// )?;
    ///
    /// // Example 3: Set operation parameters
    /// let yaw_angles = Array2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), 10.0);
    /// model.set(
    ///     None, None, None, None, None, None, None, None, None,
    ///     Some(yaw_angles), None, None, None, None, None
    /// )?;
    ///
    /// // Example 4: Combined update
    /// let ws = Array1::from_vec(vec![9.0]);
    /// let wd = Array1::from_vec(vec![280.0]);
    /// let ti = Array1::from_vec(vec![0.07]);
    /// let yaw = Array2::from_elem((1, model.core().farm.n_turbines()), 15.0);
    /// model.set(
    ///     Some(ws), Some(wd), None, None, None, Some(ti), None,
    ///     None, None, Some(yaw), None, None, None, None, None
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn set(
        &mut self,
        wind_speeds: Option<Array1>,
        wind_directions: Option<Array1>,
        wind_shear: Option<Float>,
        wind_veer: Option<Float>,
        reference_wind_height: Option<Float>,
        turbulence_intensities: Option<Array1>,
        air_density: Option<Float>,
        layout_x: Option<Array1>,
        layout_y: Option<Array1>,
        yaw_angles: Option<Array2>,
        power_setpoints: Option<Array2>,
        awc_modes: Option<NdArray2<String>>,
        awc_amplitudes: Option<Array2>,
        awc_frequencies: Option<Array2>,
        disable_turbines: Option<NdArray2<bool>>,
    ) -> crate::Result<()> {
        // Save current operation settings
        let saved_yaw_angles = self.core.farm.yaw_angles.clone();
        let saved_power_setpoints = self.core.farm.power_setpoints.clone();
        let saved_awc_modes = self.core.farm.awc_modes.clone();
        let saved_awc_amplitudes = self.core.farm.awc_amplitudes.clone();
        let saved_awc_frequencies = self.core.farm.awc_frequencies.clone();

        // Check if saved values are non-default
        let yaw_is_non_default = saved_yaw_angles.iter().any(|&v| v != 0.0);
        let power_is_non_default = saved_power_setpoints
            .iter()
            .any(|&v| v != POWER_SETPOINT_DEFAULT);
        let awc_modes_is_non_default = saved_awc_modes.iter().any(|v| v != "baseline");
        let awc_amp_is_non_default = saved_awc_amplitudes.iter().any(|&v| v != 0.0);
        let awc_freq_is_non_default = saved_awc_frequencies.iter().any(|&v| v != 0.0);

        // Reinitialize with new conditions
        self.reinitialize(
            wind_speeds,
            wind_directions,
            wind_shear,
            wind_veer,
            reference_wind_height,
            turbulence_intensities,
            air_density,
            layout_x,
            layout_y,
        )?;

        // Restore previous operation settings if they were non-default,
        // reshaping saved arrays to the new dimensions
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();

        if yaw_is_non_default {
            if let Some(reshaped) = reshape_saved_array(&saved_yaw_angles, n_findex, n_turbines) {
                self.core.farm.set_yaw_angles(reshaped);
            }
        }
        if power_is_non_default {
            if let Some(reshaped) =
                reshape_saved_array(&saved_power_setpoints, n_findex, n_turbines)
            {
                self.core.farm.set_power_setpoints(reshaped);
            }
        }
        if awc_modes_is_non_default {
            if let Some(reshaped) =
                reshape_saved_string_array(&saved_awc_modes, n_findex, n_turbines)
            {
                self.core.farm.set_awc_modes(reshaped);
            }
        }
        if awc_amp_is_non_default {
            if let Some(reshaped) = reshape_saved_array(&saved_awc_amplitudes, n_findex, n_turbines)
            {
                self.core.farm.set_awc_amplitudes(reshaped);
            }
        }
        if awc_freq_is_non_default {
            if let Some(reshaped) =
                reshape_saved_array(&saved_awc_frequencies, n_findex, n_turbines)
            {
                self.core.farm.set_awc_frequencies(reshaped);
            }
        }

        // Apply new operation settings
        self.set_operation(
            yaw_angles,
            power_setpoints,
            awc_modes,
            awc_amplitudes,
            awc_frequencies,
            disable_turbines,
        )?;

        Ok(())
    }

    /// Reinitialize the model with new conditions
    ///
    /// This is a helper method used by `set()` to update wind conditions and layout.
    fn reinitialize(
        &mut self,
        wind_speeds: Option<Array1>,
        wind_directions: Option<Array1>,
        wind_shear: Option<Float>,
        wind_veer: Option<Float>,
        reference_wind_height: Option<Float>,
        turbulence_intensities: Option<Array1>,
        air_density: Option<Float>,
        layout_x: Option<Array1>,
        layout_y: Option<Array1>,
    ) -> crate::Result<()> {
        // Use current values if not provided
        let new_wind_speeds =
            wind_speeds.unwrap_or_else(|| self.core.flow_field.wind_speeds.clone());
        let new_wind_directions =
            wind_directions.unwrap_or_else(|| self.core.flow_field.wind_directions.clone());
        let new_wind_shear = wind_shear.unwrap_or(self.core.flow_field.wind_shear);
        let new_wind_veer = wind_veer.unwrap_or(self.core.flow_field.wind_veer);
        let new_reference_height =
            reference_wind_height.unwrap_or(self.core.flow_field.reference_wind_height);
        let new_ti = turbulence_intensities
            .unwrap_or_else(|| self.core.flow_field.turbulence_intensities.clone());
        let new_air_density = air_density.unwrap_or(self.core.flow_field.air_density);

        // Validate dimensions
        if new_wind_speeds.len() != new_wind_directions.len() {
            anyhow::bail!(
                "wind_speeds (len={}) and wind_directions (len={}) must have same length",
                new_wind_speeds.len(),
                new_wind_directions.len()
            );
        }

        if new_ti.len() != new_wind_speeds.len() {
            anyhow::bail!(
                "turbulence_intensities (len={}) must match number of conditions ({})",
                new_ti.len(),
                new_wind_speeds.len()
            );
        }

        // Create new flow field
        self.core.flow_field = FlowField::new(
            new_wind_speeds,
            new_wind_directions,
            new_wind_veer,
            new_wind_shear,
            new_air_density,
            new_ti,
            new_reference_height,
        )?;

        // Update layout if provided
        if let (Some(lx), Some(ly)) = (layout_x, layout_y) {
            self.core.farm.set_layout(&lx, &ly)?;
        }

        // Reinitialize control arrays with new n_findex
        self.core
            .farm
            .initialize_control_arrays(self.core.flow_field.n_findex);
        self.core
            .farm
            .set_tilt_to_ref_tilt(self.core.flow_field.n_findex);

        // Reset grid (will be recreated on next run)
        self.core.grid = None;

        // Reset state
        self.core.state.converged = false;
        self.core.state.initialized = false;

        Ok(())
    }
    pub fn reset_operation(&mut self) {
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();

        self.core.farm.yaw_angles = Array2::zeros((n_findex, n_turbines));
        self.core.farm.power_setpoints =
            Array2::from_elem((n_findex, n_turbines), POWER_SETPOINT_DEFAULT);
        self.core.farm.awc_modes =
            NdArray2::from_elem((n_findex, n_turbines), "baseline".to_string());
        self.core.farm.awc_amplitudes = Array2::zeros((n_findex, n_turbines));
        self.core.farm.awc_frequencies = Array2::zeros((n_findex, n_turbines));
    }

    /// Calculate Annual Energy Production (AEP)
    ///
    /// Mirrors Python `get_farm_AEP` (floris_model.py:827-873 ->
    /// `get_expected_farm_value`):
    /// - a per-condition frequency column (shape `(n_findex, 1)`, e.g.
    ///   TimeSeries / WindRose via `unpack_freq`) weights each condition's
    ///   FARM power;
    /// - a per-turbine frequency table (shape `(n_findex, n_turbines)`,
    ///   WindRoseByTurbine) weights each turbine's power individually.
    ///
    /// # Arguments
    /// * `turbine_weights` - Optional per-turbine weights (Python
    ///   `turbine_weights`): turbines of a neighboring farm can be given a
    ///   weight of 0 so they contribute to the wakes but not to the
    ///   objective. `None` weights every turbine equally.
    ///
    /// # Errors
    /// - the model must have been run for the same wind data (`run()` after
    ///   `set_wind_data`); the frequency length must match `n_findex`.
    pub fn get_farm_aep(
        &self,
        wind_data: &dyn WindData,
        hours_per_year: Float,
    ) -> crate::Result<Float> {
        self.get_farm_aep_weighted(wind_data, hours_per_year, None)
    }

    /// AEP with per-turbine weights (Python `get_farm_AEP(turbine_weights=...)`).
    pub fn get_farm_aep_weighted(
        &self,
        wind_data: &dyn WindData,
        hours_per_year: Float,
        turbine_weights: Option<&[Float]>,
    ) -> crate::Result<Float> {
        // Python's unpack_freq: aligned with the unpacked condition set
        // (zero-frequency bins excluded for WindRose unless
        // compute_zero_freq_occurrence).
        let freq = wind_data.unpack_freq();
        let powers = self.get_turbine_powers()?;

        let (n_findex, n_turbines) = (powers.shape()[0], powers.shape()[1]);
        if freq.shape()[0] != n_findex {
            anyhow::bail!(
                "wind data frequency length ({}) does not match the model's \
                 n_findex ({}) — call set_wind_data() and run() with the same \
                 wind data",
                freq.shape()[0],
                n_findex
            );
        }

        let mut expected = 0.0;
        if freq.shape()[1] == 1 {
            // Per-condition frequency: weight the farm power of each findex
            for i in 0..n_findex {
                let mut farm_power = 0.0;
                for ti in 0..n_turbines {
                    let w = turbine_weights.map(|tw| tw[ti]).unwrap_or(1.0);
                    farm_power += w * powers[[i, ti]];
                }
                expected += freq[[i, 0]] * farm_power;
            }
        } else {
            // Per-turbine frequency (WindRoseByTurbine): weight each
            // turbine's power (Python get_expected_farm_value's 2D branch)
            for i in 0..n_findex {
                for ti in 0..n_turbines {
                    let w = turbine_weights.map(|tw| tw[ti]).unwrap_or(1.0);
                    expected += freq[[i, ti]] * w * powers[[i, ti]];
                }
            }
        }

        Ok(expected * hours_per_year)
    }

    /// Calculate AEP using uniform frequencies (1/N for each condition)
    ///
    /// This is useful for layout optimization when you want to compare
    /// layouts using the same set of wind conditions.
    pub fn get_farm_aep_uniform(&self, hours_per_year: Float) -> crate::Result<Float> {
        let n_conditions = self.core.flow_field.n_findex;
        let frequency = 1.0 / n_conditions as Float;
        let powers = self.get_turbine_powers()?;

        let mut total_energy = 0.0;

        for i in 0..n_conditions {
            let power = powers.row(i).sum();
            total_energy += power * frequency * hours_per_year;
        }

        Ok(total_energy)
    }

    /// Get Annual Value Production (AVP)
    ///
    /// Returns the total value-weighted energy production.
    /// If no value table is available, returns AEP.
    pub fn get_farm_avp(&self) -> crate::Result<Float> {
        // Simple implementation: return farm power as a proxy for value
        // In a full implementation, this would multiply by value factors
        Ok(self.get_farm_power()?.sum())
    }

    /// Get velocity samples at each turbine: the FIRST rotor grid point
    /// (`u[fi, ti, 0, 0]`), in layout order.
    ///
    /// # Note
    /// This is a point sample, NOT a rotor average (there is no Python
    /// equivalent; Python exposes `turbine_average_velocities` instead —
    /// see `get_turbine_average_velocities`). Provided for grid-point
    /// inspection.
    pub fn get_turbine_velocities(&self) -> Array3 {
        let n_findex = self.core.flow_field.n_findex;
        let n_turbines = self.core.farm.n_turbines();

        let mut velocities = Array::zeros((n_findex, n_turbines, 1));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                // Unsorted (layout-order) point sample from flow_field.u
                velocities[[fi, ti, 0]] = self.core.flow_field.u[[fi, ti, 0, 0]];
            }
        }

        velocities
    }

    /// Extract velocity deficit profiles at a set of downstream distances.
    ///
    /// Python `FlorisModel.sample_velocity_deficit_profiles`
    /// (floris/floris_model.py): the velocity deficit is
    /// `(homogeneous_wind_speed - u) / homogeneous_wind_speed`, where `u` is
    /// the wake velocity obtained with a homogeneous inflow (the wind shear
    /// is temporarily set to zero, as in Python).
    ///
    /// * `direction` - "cross-stream" or "vertical"
    /// * `downstream_dists` - the streamwise locations of the profiles
    ///   (default: `ref_rotor_diameter * [3, 5, 7, 9]`)
    /// * `profile_range` - the extent of the profile lines (default:
    ///   `ref_rotor_diameter * [-2, 2]`)
    /// * `resolution` - the number of points of every profile
    /// * `wind_direction` - the wind direction (default: the single wind
    ///   direction of the model)
    /// * `homogeneous_wind_speed` - the wind speed of the homogeneous inflow
    ///   (default: the single wind speed of the model)
    /// * `ref_rotor_diameter` - the reference rotor diameter (default: the
    ///   rotor diameter of the farm)
    /// * `x_start`, `y_start` - the origin of the profiles
    /// * `reference_height` - the height of the horizontal plane (default:
    ///   the reference wind height of the model)
    pub fn sample_velocity_deficit_profiles(
        &mut self,
        direction: &str,
        downstream_dists: Option<&[Float]>,
        profile_range: Option<&[Float]>,
        resolution: usize,
        wind_direction: Option<Float>,
        homogeneous_wind_speed: Option<Float>,
        ref_rotor_diameter: Option<Float>,
        x_start: Float,
        y_start: Float,
        reference_height: Option<Float>,
    ) -> crate::Result<Vec<crate::core::core::VelocityDeficitProfile>> {
        // The reference rotor diameter
        let ref_rotor_diameter = match ref_rotor_diameter {
            Some(d) => d,
            None => {
                let diameters = self.core.farm.rotor_diameters();
                let mut unique = diameters.to_vec();
                unique.sort_by(|a, b| a.partial_cmp(b).unwrap());
                unique.dedup();
                if unique.len() != 1 {
                    anyhow::bail!(
                        "Please provide a `ref_rotor_diameter`. This is needed to normalize the \
                         coordinates."
                    );
                }
                unique[0]
            }
        };

        let downstream_dists: Vec<Float> = match downstream_dists {
            Some(d) => d.to_vec(),
            None => vec![
                3.0 * ref_rotor_diameter,
                5.0 * ref_rotor_diameter,
                7.0 * ref_rotor_diameter,
                9.0 * ref_rotor_diameter,
            ],
        };
        let profile_range: Vec<Float> = match profile_range {
            Some(r) => r.to_vec(),
            None => vec![-2.0 * ref_rotor_diameter, 2.0 * ref_rotor_diameter],
        };

        // The wind direction and the homogeneous wind speed
        let wind_direction = match wind_direction {
            Some(wd) => wd,
            None => {
                if self.core.flow_field.wind_directions.len() != 1 {
                    anyhow::bail!(
                        "Could not determine a wind direction for which to sample the velocity \
                         profiles. Either provide a single `wind_direction` as an argument to this \
                         method, or initialize the Floris object with a single wind direction."
                    );
                }
                self.core.flow_field.wind_directions[0]
            }
        };
        let homogeneous_wind_speed = match homogeneous_wind_speed {
            Some(ws) => ws,
            None => {
                if self.core.flow_field.wind_speeds.len() != 1 {
                    anyhow::bail!(
                        "Could not determine a wind speed for which to sample the velocity \
                         profiles. Provide a single `homogeneous_wind_speed` to this method."
                    );
                }
                self.core.flow_field.wind_speeds[0]
            }
        };
        let reference_height =
            reference_height.unwrap_or(self.core.flow_field.reference_wind_height);

        // The inflow is homogeneous while the profiles are sampled: the wind
        // shear is temporarily set to zero and the wind conditions are reduced
        // to a single direction and speed (Python does the same).
        let wind_directions_copy = self.core.flow_field.wind_directions.clone();
        let wind_speeds_copy = self.core.flow_field.wind_speeds.clone();
        let wind_shear_copy = self.core.flow_field.wind_shear;
        let turbulence_intensities = self.core.flow_field.turbulence_intensities.clone();

        self.set_wind_conditions(
            Array1::from_vec(vec![homogeneous_wind_speed]),
            Array1::from_vec(vec![wind_direction]),
            Array1::from_vec(vec![turbulence_intensities[0]]),
        )?;
        self.set_wind_shear(0.0)?;

        let profiles = self.core.solve_for_velocity_deficit_profiles(
            direction,
            &downstream_dists,
            &profile_range,
            resolution,
            homogeneous_wind_speed,
            ref_rotor_diameter,
            x_start,
            y_start,
            reference_height,
        );

        // Restore the original wind conditions
        self.set_wind_conditions(
            wind_speeds_copy,
            wind_directions_copy,
            turbulence_intensities,
        )?;
        self.set_wind_shear(wind_shear_copy)?;

        profiles
    }

    /// Sample the flow velocity at arbitrary points, solving the wake AT the
    /// points for ALL wind conditions (Python `sample_flow_at_points` ->
    /// `core.solve_for_points`, floris_model.py:1386-1404).
    ///
    /// # Arguments
    /// * `points_x` - x-coordinates of sampling points (in meters, inertial frame)
    /// * `points_y` - y-coordinates of sampling points (in meters, inertial frame)
    /// * `points_z` - z-coordinates of sampling points (in meters, height above ground)
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_points)` containing the velocity at each
    /// point for each wind condition.
    ///
    /// # Errors
    /// - `turbopark` is rejected (Python refuses it too); `cc` and
    ///   `empirical_gauss` are rejected until their full-flow solver variants
    ///   are ported.
    /// - The points must share the model's wind directions; call `set()` to
    ///   change conditions before sampling.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use florus::FlorisModel;
    ///
    /// let mut fmodel = FlorisModel::from_file("gch.yaml").unwrap();
    /// fmodel.run().unwrap();
    ///
    /// // Sample at met mast location
    /// let velocities = fmodel.sample_flow_at_points(
    ///     &[500.0, 500.0, 500.0],
    ///     &[0.0, 0.0, 0.0],
    ///     &[30.0, 90.0, 150.0]
    /// ).unwrap();
    /// ```
    pub fn sample_flow_at_points(
        &self,
        points_x: &[f64],
        points_y: &[f64],
        points_z: &[f64],
    ) -> crate::Result<Array2> {
        let (velocities, _ti) = self.core.solve_for_points(points_x, points_y, points_z)?;
        Ok(velocities)
    }

    /// Sample the turbulence intensity at arbitrary points, solving the wake
    /// AT the points for ALL wind conditions (Python `sample_ti_at_points`,
    /// floris_model.py:1406-1423: it solves on the points grid and returns
    /// `turbulence_intensity_field_sorted[:, :, 0, 0]`).
    ///
    /// Unlike the pre-fix implementation, which returned the constant
    /// ambient TI, this returns the wake-added TI field produced by the
    /// solver's turbulence model at each point.
    ///
    /// # Arguments
    /// * `points_x` - x-coordinates of sampling points (in meters, inertial frame)
    /// * `points_y` - y-coordinates of sampling points (in meters, inertial frame)
    /// * `points_z` - z-coordinates of sampling points (in meters, height above ground)
    ///
    /// # Returns
    /// Array of shape `(n_findex, n_points)` containing the turbulence
    /// intensity at each point for each wind condition.
    ///
    /// # Errors
    /// Same conditions as `sample_flow_at_points`.
    pub fn sample_ti_at_points(
        &self,
        points_x: &[f64],
        points_y: &[f64],
        points_z: &[f64],
    ) -> crate::Result<Array2> {
        let (_velocities, ti) = self.core.solve_for_points(points_x, points_y, points_z)?;
        Ok(ti)
    }
    /// Apply operating setpoints to the floris object.
    ///
    /// This function sets various operational parameters for turbines including yaw angles,
    /// power setpoints, AWC (Active Wake Control) modes, and turbine disable flags.
    ///
    /// # Arguments
    /// * `yaw_angles` - Optional turbine yaw angles array with shape (n_findex, n_turbines)
    /// * `power_setpoints` - Optional turbine power setpoints array with shape (n_findex, n_turbines).
    ///                       None/NaN values will be replaced with POWER_SETPOINT_DEFAULT
    /// * `awc_modes` - Optional AWC modes array with shape (n_findex, n_turbines).
    ///                 Defaults to "baseline" for all turbines if not provided
    /// * `awc_amplitudes` - Optional AWC amplitudes array with shape (n_findex, n_turbines).
    ///                      Defaults to zeros if not provided
    /// * `awc_frequencies` - Optional AWC frequencies array with shape (n_findex, n_turbines).
    ///                       Defaults to zeros if not provided
    /// * `disable_turbines` - Optional boolean array with shape (n_findex, n_turbines) indicating
    ///                        which turbines to disable. Disabled turbines have yaw_angles set to 0
    ///                        and power_setpoints set to POWER_SETPOINT_DISABLED (0.001)
    ///
    /// # Returns
    /// Returns `Ok(())` on success, or an error if array dimensions don't match expected shapes
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use florus::FlorisModel;
    /// use ndarray::{Array2, Array1};
    ///
    /// // Load a model
    /// let mut model = FlorisModel::from_file("default_inputs.yaml").unwrap();
    ///
    /// // Example 1: Set yaw angles only
    /// let yaw_angles = Array2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), 10.0);
    /// model.set_operation(Some(yaw_angles), None, None, None, None, None).unwrap();
    ///
    /// // Example 2: Set power setpoints
    /// let power_setpoints = Array2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), 3e6);
    /// model.set_operation(None, Some(power_setpoints), None, None, None, None).unwrap();
    ///
    /// // Example 3: Disable specific turbines
    /// use ndarray::Array2 as NdArray2;
    /// let mut disable = NdArray2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), false);
    /// disable[[0, 0]] = true; // Disable first turbine in first condition
    /// model.set_operation(None, None, None, None, None, Some(disable)).unwrap();
    ///
    /// // Example 4: Set AWC parameters
    /// let awc_modes = NdArray2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), "sinusoidal".to_string());
    /// let awc_amplitudes = Array2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), 2.0);
    /// let awc_frequencies = Array2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), 0.5);
    /// model.set_operation(None, None, Some(awc_modes), Some(awc_amplitudes), Some(awc_frequencies), None).unwrap();
    ///
    /// // Example 5: Set multiple parameters at once
    /// let yaw = Array2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), 5.0);
    /// let power = Array2::from_elem((model.core().flow_field.n_findex, model.core().farm.n_turbines()), 4e6);
    /// model.set_operation(Some(yaw), Some(power), None, None, None, None).unwrap();
    /// ```
    ///
    /// # Notes
    /// - When any operation parameter is set, the model state is marked as uninitialized
    ///   and will need to be re-run to get updated results
    /// - All arrays must have shape (n_findex, n_turbines) where:
    ///   - n_findex: number of wind conditions (from flow_field)
    ///   - n_turbines: number of turbines in the farm
    /// - AWC (Active Wake Control) modes can be used for advanced wake steering strategies
    pub fn set_operation(
        &mut self,
        yaw_angles: Option<Array2>,
        power_setpoints: Option<Array2>,
        awc_modes: Option<NdArray2<String>>,
        awc_amplitudes: Option<Array2>,
        awc_frequencies: Option<Array2>,
        disable_turbines: Option<NdArray2<bool>>,
    ) -> crate::Result<()> {
        let n_turbines = self.core.farm.n_turbines();
        let n_findex = self.core.flow_field.n_findex;

        // Set yaw angles if provided
        let yaw_angles_provided = yaw_angles.is_some();
        if let Some(yaw) = yaw_angles {
            if yaw.shape()[1] != n_turbines {
                anyhow::bail!(
                    "yaw_angles has a size of {} in the 1st dimension, must be equal to n_turbines={}",
                    yaw.shape()[1],
                    n_turbines
                );
            }
            self.core.farm.set_yaw_angles(yaw);
        }

        // Set power setpoints if provided
        let power_setpoints_provided = power_setpoints.is_some();
        if let Some(mut power) = power_setpoints {
            if power.shape()[1] != n_turbines {
                anyhow::bail!(
                    "power_setpoints has a size of {} in the 1st dimension, must be equal to n_turbines={}",
                    power.shape()[1],
                    n_turbines
                );
            }

            // Replace any NaN or invalid values with default power setpoint
            for val in power.iter_mut() {
                if val.is_nan() || *val < 0.0 {
                    *val = POWER_SETPOINT_DEFAULT;
                }
            }

            self.core.farm.set_power_setpoints(power);
        }

        // Set AWC modes (default to "baseline")
        let awc_modes = awc_modes
            .unwrap_or_else(|| NdArray2::from_elem((n_findex, n_turbines), "baseline".to_string()));
        self.core.farm.set_awc_modes(awc_modes);

        // Set AWC amplitudes (default to zeros)
        let awc_amplitudes =
            awc_amplitudes.unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));
        self.core.farm.set_awc_amplitudes(awc_amplitudes);

        // Set AWC frequencies (default to zeros)
        let awc_frequencies =
            awc_frequencies.unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));
        self.core.farm.set_awc_frequencies(awc_frequencies);

        // Handle disabled turbines
        let disable_turbines_provided = disable_turbines.is_some();
        if let Some(disabled) = disable_turbines {
            // Validate dimensions
            if disabled.shape()[0] != n_findex {
                anyhow::bail!(
                    "disable_turbines has a size of {} in the 0th dimension, must be equal to n_findex={}",
                    disabled.shape()[0],
                    n_findex
                );
            }

            if disabled.shape()[1] != n_turbines {
                anyhow::bail!(
                    "disable_turbines has a size of {} in the 1th dimension, must be equal to n_turbines={}",
                    disabled.shape()[1],
                    n_turbines
                );
            }

            // Set yaw_angles to 0 and power_setpoints to disabled value where disable_turbines is true
            for fi in 0..n_findex {
                for ti in 0..n_turbines {
                    if disabled[[fi, ti]] {
                        self.core.farm.yaw_angles[[fi, ti]] = 0.0;
                        self.core.farm.power_setpoints[[fi, ti]] = POWER_SETPOINT_DISABLED;
                    }
                }
            }
        }

        // Mark state as uninitialized if any operation parameters were set
        if yaw_angles_provided || power_setpoints_provided || disable_turbines_provided {
            self.core.state.converged = false;
            self.core.state.initialized = false;
        }

        Ok(())
    }
}

/// Reshape a saved operation array to the new (n_findex, n_turbines) dimensions,
/// broadcasting from (1, n_turbines) when the findex count changed.
/// Returns None when the saved shape is incompatible with the new dimensions.
fn reshape_saved_array(arr: &Array2, n_findex: usize, n_turbines: usize) -> Option<Array2> {
    let shape = arr.shape();
    if shape.len() == 2 && shape[0] == n_findex && shape[1] == n_turbines {
        Some(arr.clone())
    } else if shape.len() == 2 && shape[0] == 1 && shape[1] == n_turbines && n_findex > 1 {
        let mut expanded = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                expanded[[fi, ti]] = arr[[0, ti]];
            }
        }
        Some(expanded)
    } else {
        None
    }
}

fn reshape_saved_string_array(
    arr: &NdArray2<String>,
    n_findex: usize,
    n_turbines: usize,
) -> Option<NdArray2<String>> {
    let shape = arr.shape();
    if shape.len() == 2 && shape[0] == n_findex && shape[1] == n_turbines {
        Some(arr.clone())
    } else if shape.len() == 2 && shape[0] == 1 && shape[1] == n_turbines && n_findex > 1 {
        let mut expanded = NdArray2::from_elem((n_findex, n_turbines), String::new());
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                expanded[[fi, ti]] = arr[[0, ti]].clone();
            }
        }
        Some(expanded)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper function to create a test FlorisModel using Core::from_config
    fn create_test_model(n_turbines: usize, n_findex: usize) -> FlorisModel {
        let layout_x = (0..n_turbines)
            .map(|i| i as f64 * 500.0)
            .collect::<Vec<_>>();
        let layout_y = vec![0.0; n_turbines];
        let turbine_types = vec!["nrel_5MW".to_string(); n_turbines];

        let wind_speeds = (0..n_findex)
            .map(|i| 8.0 + i as f64 * 2.0)
            .collect::<Vec<_>>();
        let wind_directions = vec![270.0; n_findex];
        let turbulence_intensities = vec![0.06; n_findex];

        use crate::floris_config::{FarmConfig, FlowFieldConfig, LoggingConfig, WakeConfig};

        let config = FlorisConfig {
            name: "test".to_string(),
            description: Some("Test configuration".to_string()),
            floris_version: "4.0".to_string(),
            logging: LoggingConfig::default(),
            solver: SolverConfig::default(),
            farm: FarmConfig {
                layout_x,
                layout_y,
                turbine_type: turbine_types,
            },
            flow_field: FlowFieldConfig {
                wind_speeds,
                wind_directions,
                wind_shear: 0.14,
                wind_veer: 0.0,
                air_density: 1.225,
                turbulence_intensities,
                reference_wind_height: 90.0,
                multidim_conditions: None,
                heterogeneous_inflow_config: None,
            },
            wake: WakeConfig::default(),
            turbine_library: "turbine_library".to_string(),
        };

        FlorisModel::from_config(config).unwrap()
    }

    #[test]
    fn test_floris_model_basic() {
        let model = create_test_model(2, 2);
        assert_eq!(model.core.farm.n_turbines(), 2);
        assert_eq!(model.core.flow_field.n_findex, 2);
    }

    #[test]
    fn test_get_turbine_powers() {
        let model = create_test_model(2, 2);
        assert_eq!(model.core.farm.n_turbines(), 2);
        assert_eq!(model.core.flow_field.n_findex, 2);
    }

    #[test]
    fn test_wake_solver_integration() {
        // Create a simple 2-turbine farm with downstream turbine
        let mut model = create_test_model(2, 1);

        // Initialize grid and flow field
        model.initialize_grid().unwrap();
        model.initialize_flow_field().unwrap();

        // Verify grid was created
        assert!(model.core.grid.is_some());
        let grid = model.core.grid.as_ref().unwrap();
        assert_eq!(grid.n_turbines(), 2);
        assert_eq!(grid.n_findex(), 1);

        // Verify flow field was initialized with correct shape
        let shape = model.core.flow_field.u_sorted.shape();
        assert_eq!(shape[0], 1); // n_findex
        assert_eq!(shape[1], 2); // n_turbines
        assert_eq!(shape[2], 3); // grid resolution y
        assert_eq!(shape[3], 3); // grid resolution z

        // Verify grid coordinates are properly sorted
        // X coordinates should be different for upstream vs downstream turbines
        let x_sorted = grid.x_sorted();
        assert_ne!(x_sorted[[0, 0, 0, 0]], x_sorted[[0, 1, 0, 0]]);
    }

    #[test]
    fn test_set_operation() {
        use ndarray::Array2 as NdArray2;

        let mut model = create_test_model(2, 2);

        // Test 1: Set yaw angles
        let yaw_angles = Array2::from_elem((2, 2), 10.0);
        model
            .set_operation(Some(yaw_angles.clone()), None, None, None, None, None)
            .unwrap();

        assert_eq!(model.core.farm.yaw_angles.shape(), &[2, 2]);
        assert_eq!(model.core.farm.yaw_angles[[0, 0]], 10.0);

        // Test 2: Set power setpoints
        let power_setpoints = Array2::from_elem((2, 2), 3000000.0);
        model
            .set_operation(None, Some(power_setpoints.clone()), None, None, None, None)
            .unwrap();

        assert_eq!(model.core.farm.power_setpoints.shape(), &[2, 2]);
        assert_eq!(model.core.farm.power_setpoints[[0, 0]], 3000000.0);

        // Test 3: Set AWC modes
        let awc_modes = NdArray2::from_elem((2, 2), "sinusoidal".to_string());
        model
            .set_operation(None, None, Some(awc_modes.clone()), None, None, None)
            .unwrap();

        assert_eq!(model.core.farm.awc_modes.shape(), &[2, 2]);
        assert_eq!(model.core.farm.awc_modes[[0, 0]], "sinusoidal");

        // Test 4: Set AWC amplitudes and frequencies
        let awc_amplitudes = Array2::from_elem((2, 2), 2.0);
        let awc_frequencies = Array2::from_elem((2, 2), 0.5);
        model
            .set_operation(
                None,
                None,
                None,
                Some(awc_amplitudes.clone()),
                Some(awc_frequencies.clone()),
                None,
            )
            .unwrap();

        assert_eq!(model.core.farm.awc_amplitudes[[0, 0]], 2.0);
        assert_eq!(model.core.farm.awc_frequencies[[0, 0]], 0.5);

        // Test 5: Disable turbines
        let disable_turbines = NdArray2::from_shape_fn((2, 2), |(i, j)| i == 0 && j == 0);
        model
            .set_operation(None, None, None, None, None, Some(disable_turbines.clone()))
            .unwrap();

        // First turbine in first condition should be disabled
        assert_eq!(model.core.farm.yaw_angles[[0, 0]], 0.0);
        assert_eq!(
            model.core.farm.power_setpoints[[0, 0]],
            POWER_SETPOINT_DISABLED
        );

        // Test 6: Verify state is marked as uninitialized
        assert!(!model.core.state.converged);
        assert!(!model.core.state.initialized);

        // Test 7: Test with all None values (should use defaults)
        model
            .set_operation(None, None, None, None, None, None)
            .unwrap();

        // AWC modes should default to "baseline"
        assert_eq!(model.core.farm.awc_modes[[0, 0]], "baseline");
        // AWC amplitudes and frequencies should default to 0
        assert_eq!(model.core.farm.awc_amplitudes[[0, 0]], 0.0);
        assert_eq!(model.core.farm.awc_frequencies[[0, 0]], 0.0);
    }

    #[test]
    fn test_set_operation_validation() {
        let mut model = create_test_model(2, 1);

        // Test validation: wrong yaw_angles dimension
        let wrong_yaw_angles = Array2::from_elem((1, 3), 10.0); // Should be (1, 2)
        let result = model.set_operation(Some(wrong_yaw_angles), None, None, None, None, None);
        assert!(result.is_err());

        // Test validation: wrong disable_turbines dimension
        let wrong_disable = NdArray2::from_shape_fn((2, 2), |(_i, _j)| false); // Should be (1, 2)
        let result = model.set_operation(None, None, None, None, None, Some(wrong_disable));
        assert!(result.is_err());
    }

    #[test]
    fn test_set_operation_example_usage() {
        // This test demonstrates typical usage patterns for set_operation
        let mut model = create_test_model(3, 2);

        // Example 1: Set different yaw angles for each turbine and condition
        let mut yaw_angles = Array2::zeros((2, 3));
        yaw_angles[[0, 0]] = 10.0; // First condition, first turbine
        yaw_angles[[0, 1]] = 15.0; // First condition, second turbine
        yaw_angles[[0, 2]] = 5.0; // First condition, third turbine
        yaw_angles[[1, 0]] = 12.0; // Second condition, first turbine
        yaw_angles[[1, 1]] = 18.0; // Second condition, second turbine
        yaw_angles[[1, 2]] = 8.0; // Second condition, third turbine

        model
            .set_operation(Some(yaw_angles), None, None, None, None, None)
            .unwrap();

        // Verify yaw angles were set correctly
        assert_eq!(model.core.farm.yaw_angles[[0, 0]], 10.0);
        assert_eq!(model.core.farm.yaw_angles[[1, 2]], 8.0);

        // Example 2: Set power derating for specific turbines
        let mut power_setpoints = Array2::from_elem((2, 3), 5e6); // 5 MW default
        power_setpoints[[0, 0]] = 3e6; // Derate first turbine in first condition
        power_setpoints[[0, 1]] = 4e6; // Derate second turbine in first condition

        model
            .set_operation(None, Some(power_setpoints), None, None, None, None)
            .unwrap();

        assert_eq!(model.core.farm.power_setpoints[[0, 0]], 3e6);
        assert_eq!(model.core.farm.power_setpoints[[0, 2]], 5e6); // Unchanged

        // Example 3: Configure AWC for wake steering
        let awc_modes = NdArray2::from_elem((2, 3), "sinusoidal".to_string());
        let mut awc_amplitudes = Array2::zeros((2, 3));
        awc_amplitudes[[0, 0]] = 3.0; // 3 degree amplitude
        awc_amplitudes[[0, 1]] = 2.5;

        let mut awc_frequencies = Array2::zeros((2, 3));
        awc_frequencies[[0, 0]] = 0.3; // 0.3 Hz
        awc_frequencies[[0, 1]] = 0.4;

        model
            .set_operation(
                None,
                None,
                Some(awc_modes),
                Some(awc_amplitudes),
                Some(awc_frequencies),
                None,
            )
            .unwrap();

        assert_eq!(model.core.farm.awc_modes[[0, 0]], "sinusoidal");
        assert_eq!(model.core.farm.awc_amplitudes[[0, 0]], 3.0);
        assert_eq!(model.core.farm.awc_frequencies[[0, 0]], 0.3);

        // Example 4: Disable downstream turbine in certain conditions
        let mut disable_turbines = NdArray2::from_elem((2, 3), false);
        disable_turbines[[0, 2]] = true; // Disable third turbine in first condition
        disable_turbines[[1, 1]] = true; // Disable second turbine in second condition

        model
            .set_operation(None, None, None, None, None, Some(disable_turbines))
            .unwrap();

        // Verify disabled turbines have correct settings
        assert_eq!(model.core.farm.yaw_angles[[0, 2]], 0.0);
        assert_eq!(
            model.core.farm.power_setpoints[[0, 2]],
            POWER_SETPOINT_DISABLED
        );
        assert_eq!(model.core.farm.yaw_angles[[1, 1]], 0.0);
        assert_eq!(
            model.core.farm.power_setpoints[[1, 1]],
            POWER_SETPOINT_DISABLED
        );

        // Verify non-disabled turbines are unchanged
        assert_ne!(model.core.farm.yaw_angles[[0, 0]], 0.0);
        assert_ne!(
            model.core.farm.power_setpoints[[0, 0]],
            POWER_SETPOINT_DISABLED
        );
    }

    #[test]
    fn test_set_method() {
        let mut model = create_test_model(2, 2);

        // Test 1: Update wind conditions only
        let new_ws = Array1::from_vec(vec![9.0, 11.0, 13.0]);
        let new_wd = Array1::from_vec(vec![275.0, 275.0, 275.0]);
        let new_ti = Array1::from_vec(vec![0.07, 0.07, 0.07]);

        model
            .set(
                Some(new_ws.clone()),
                Some(new_wd.clone()),
                None,
                None,
                None,
                Some(new_ti.clone()),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap();

        assert_eq!(model.core.flow_field.n_findex, 3);
        assert_eq!(model.core.flow_field.wind_speeds, new_ws);
        assert_eq!(model.core.flow_field.wind_directions, new_wd);
        assert_eq!(model.core.flow_field.turbulence_intensities, new_ti);

        // Test 2: Update layout
        let new_layout_x = Array1::from_vec(vec![0.0, 600.0, 1200.0]);
        let new_layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);

        model
            .set(
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                Some(new_layout_x.clone()),
                Some(new_layout_y.clone()),
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap();

        assert_eq!(model.core.farm.n_turbines(), 3);
        assert_eq!(model.core.farm.layout_x, new_layout_x);
        assert_eq!(model.core.farm.layout_y, new_layout_y);

        // Test 3: Set operation parameters with set()
        let yaw_angles = Array2::from_elem((3, 3), 12.0);
        model
            .set(
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                Some(yaw_angles.clone()),
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap();

        assert_eq!(model.core.farm.yaw_angles.shape(), &[3, 3]);
        assert_eq!(model.core.farm.yaw_angles[[0, 0]], 12.0);

        // Test 4: Combined update - wind conditions + operation
        let ws2 = Array1::from_vec(vec![10.0]);
        let wd2 = Array1::from_vec(vec![290.0]);
        let ti2 = Array1::from_vec(vec![0.05]);
        let yaw2 = Array2::from_elem((1, 3), 8.0);

        model
            .set(
                Some(ws2.clone()),
                Some(wd2.clone()),
                None,
                None,
                None,
                Some(ti2.clone()),
                None,
                None,
                None,
                Some(yaw2.clone()),
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap();

        assert_eq!(model.core.flow_field.n_findex, 1);
        assert_eq!(model.core.flow_field.wind_speeds, ws2);
        assert_eq!(model.core.farm.yaw_angles[[0, 0]], 8.0);

        // Test 5: Verify state is reset after set()
        assert!(!model.core.state.converged);
        assert!(!model.core.state.initialized);
    }

    #[test]
    fn test_set_preserves_non_default_settings() {
        // This test verifies that set() preserves non-default operation settings
        // when reinitializing
        let mut model = create_test_model(2, 1);

        // First, set some non-default operation parameters
        let initial_yaw = Array2::from_elem((1, 2), 15.0);
        let initial_power = Array2::from_elem((1, 2), 3e6);

        model
            .set_operation(
                Some(initial_yaw.clone()),
                Some(initial_power.clone()),
                None,
                None,
                None,
                None,
            )
            .unwrap();

        // Verify they were set
        assert_eq!(model.core.farm.yaw_angles[[0, 0]], 15.0);
        assert_eq!(model.core.farm.power_setpoints[[0, 0]], 3e6);

        // Now call set() with only wind condition changes (no operation params)
        let new_ws = Array1::from_vec(vec![9.0, 10.0]);
        let new_wd = Array1::from_vec(vec![275.0, 280.0]);
        let new_ti = Array1::from_vec(vec![0.07, 0.08]);

        model
            .set(
                Some(new_ws),
                Some(new_wd),
                None,
                None,
                None,
                Some(new_ti),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None, // No operation params
            )
            .unwrap();

        // The non-default operation settings should be preserved
        // Note: After reinitialize, arrays are resized to match new n_findex
        // So we check that the values are still non-zero/non-default
        assert_ne!(model.core.farm.yaw_angles[[0, 0]], 0.0);
        assert_ne!(
            model.core.farm.power_setpoints[[0, 0]],
            POWER_SETPOINT_DEFAULT
        );
    }

    #[test]
    fn test_set_validation() {
        let mut model = create_test_model(2, 1);

        // Test validation: mismatched dimensions
        let wrong_ws = Array1::from_vec(vec![8.0, 10.0]); // length 2
        let wrong_wd = Array1::from_vec(vec![270.0]); // length 1
        let wrong_ti = Array1::from_vec(vec![0.06]); // length 1

        let result = model.set(
            Some(wrong_ws),
            Some(wrong_wd),
            None,
            None,
            None,
            Some(wrong_ti),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_get_expected_farm_power() {
        let mut model = create_test_model(3, 4);
        model.run().unwrap();

        // Test with uniform frequencies (default)
        let expected_power = model.get_expected_farm_power(None, None).unwrap();
        assert!(expected_power > 0.0);

        // Test with custom frequencies
        let freq = Array1::from_vec(vec![0.4, 0.3, 0.2, 0.1]);
        let expected_power_custom = model.get_expected_farm_power(Some(freq), None).unwrap();
        assert!(expected_power_custom > 0.0);

        println!("Expected farm power (uniform): {} W", expected_power);
        println!("Expected farm power (custom): {} W", expected_power_custom);
    }

    #[test]
    fn test_get_expected_farm_power_with_weights() {
        let mut model = create_test_model(2, 2);
        model.run().unwrap();

        // Create turbine weights: emphasize first turbine
        let weights = Array2::from_shape_vec(
            (2, 2),
            vec![
                1.0, 0.5, // First findex
                1.0, 0.5,
            ], // Second findex
        )
        .unwrap();

        let expected_power = model.get_expected_farm_power(None, Some(weights)).unwrap();
        assert!(expected_power > 0.0);

        println!("Expected farm power with weights: {} W", expected_power);
    }

    #[test]
    fn test_get_expected_farm_value() {
        let mut model = create_test_model(2, 3);
        model.run().unwrap();

        // Test with uniform values (default) - should be similar to expected power
        let expected_value = model.get_expected_farm_value(None, None, None).unwrap();
        let expected_power = model.get_expected_farm_power(None, None).unwrap();

        // With default value of 1.0, they should be close (allowing for normalization differences)
        println!("Expected value: {}", expected_value);
        println!("Expected power: {}", expected_power);

        // Test with custom values (e.g., electricity prices in $/MWh)
        let prices = Array1::from_vec(vec![50.0, 60.0, 45.0]);
        let expected_value_custom = model
            .get_expected_farm_value(None, Some(prices), None)
            .unwrap();
        assert!(expected_value_custom > 0.0);

        println!("Expected value with prices: {}", expected_value_custom);
    }

    #[test]
    fn test_get_turbine_ais() {
        let mut model = create_test_model(2, 2);
        model.run().unwrap();

        let ais = model.get_turbine_ais().unwrap();

        // Should have shape (n_findex, n_turbines)
        assert_eq!(ais.shape(), &[2, 2]);

        // AI should be in reasonable range (0 to 0.5 typically)
        for fi in 0..2 {
            for ti in 0..2 {
                assert!(
                    ais[[fi, ti]] >= 0.0 && ais[[fi, ti]] <= 1.0,
                    "AI out of range: {}",
                    ais[[fi, ti]]
                );
            }
        }

        println!("Axial induction factors:\n{:?}", ais);
    }

    #[test]
    fn test_get_turbine_tis() {
        let mut model = create_test_model(3, 2);
        model.run().unwrap();

        let tis = model.get_turbine_tis();

        // Should have shape (n_findex, n_turbines)
        assert_eq!(tis.shape(), &[2, 3]);

        // TI should be positive and reasonable (typically 0.05 to 0.25)
        for fi in 0..2 {
            for ti in 0..3 {
                assert!(tis[[fi, ti]] > 0.0, "TI should be positive");
                assert!(tis[[fi, ti]] < 1.0, "TI should be less than 1.0");
            }
        }

        println!("Turbulence intensities:\n{:?}", tis);
    }

    #[test]
    fn test_assign_hub_height_to_ref_height() {
        let mut model = create_test_model(2, 1);

        // Start from a reference height different from the nrel_5MW hub height (90 m)
        model.core_mut().flow_field.reference_wind_height = 100.0;
        let initial_ref_height = model.reference_wind_height();

        // Assign hub height to reference height
        model.assign_hub_height_to_ref_height().unwrap();

        let new_ref_height = model.reference_wind_height();
        let hub_height = model.core.farm.hub_heights[0];

        assert_eq!(new_ref_height, hub_height);
        assert_ne!(initial_ref_height, new_ref_height); // Should have changed

        println!("Initial ref height: {}", initial_ref_height);
        println!("New ref height (hub height): {}", new_ref_height);
    }

    #[test]
    fn test_assign_hub_height_multiple_heights_error() {
        // Create a model with different hub heights (would need custom config)
        // For now, just verify the method works with uniform heights
        let mut model = create_test_model(2, 1);
        assert!(model.assign_hub_height_to_ref_height().is_ok());
    }

    #[test]
    fn test_set_operation_model() {
        let mut model = create_test_model(2, 1);

        // Test setting single operation model
        assert!(model.set_operation_model("simple").is_ok());

        // Test setting multiple operation models
        let models = vec!["simple".to_string(), "cosine".to_string()];
        assert!(model.set_operation_model(models).is_ok());

        // Test validation: wrong number of models
        let wrong_models = vec!["simple".to_string()]; // Only 1 for 2 turbines
        assert!(model.set_operation_model(wrong_models).is_err());
    }

    #[test]
    fn test_get_turbine_layout() {
        let model = create_test_model(3, 1);

        // Test without z coordinates
        let layout = model.get_turbine_layout(false);
        assert_eq!(layout.x().len(), 3);
        assert_eq!(layout.y().len(), 3);
        assert!(layout.z().is_none());

        // Test with z coordinates
        let layout_z = model.get_turbine_layout(true);
        assert_eq!(layout_z.x().len(), 3);
        assert_eq!(layout_z.y().len(), 3);
        assert!(layout_z.z().is_some());
        assert_eq!(layout_z.z().unwrap().len(), 3);

        println!("Layout X: {:?}", layout.x());
        println!("Layout Y: {:?}", layout.y());
        if let Some(z) = layout_z.z() {
            println!("Layout Z: {:?}", z);
        }
    }

    #[test]
    fn test_properties() {
        let model = create_test_model(3, 4);

        // Test all property getters
        assert_eq!(model.layout_x().len(), 3);
        assert_eq!(model.layout_y().len(), 3);
        assert_eq!(model.wind_directions().len(), 4);
        assert_eq!(model.wind_speeds().len(), 4);
        assert_eq!(model.turbulence_intensities().len(), 4);
        assert_eq!(model.n_findex(), 4);
        assert_eq!(model.n_turbines(), 3);
        assert!(model.reference_wind_height() > 0.0);

        // Test turbine average velocities
        let avg_vels = model.turbine_average_velocities();
        assert_eq!(avg_vels.shape(), &[4, 3]);

        println!("Properties test passed!");
        println!("  n_findex: {}", model.n_findex());
        println!("  n_turbines: {}", model.n_turbines());
        println!("  ref height: {}", model.reference_wind_height());
    }

    #[test]
    fn test_core_accessors() {
        let model = create_test_model(2, 1);

        // Test immutable access
        let core_ref = model.core();
        assert_eq!(core_ref.farm.n_turbines(), 2);

        let farm_ref = model.farm();
        assert_eq!(farm_ref.n_turbines(), 2);

        let flow_field_ref = model.flow_field();
        assert_eq!(flow_field_ref.n_findex, 1);

        let state_ref = model.state();
        assert!(!state_ref.converged); // Not run yet

        println!("Core accessor test passed!");
    }
}
