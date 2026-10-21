//! The uncertain FLORIS model: the mean power of a farm for a distribution
//! of wind direction deviations.
//!
//! Python `floris/uncertain_floris_model.py`: the wind conditions of the
//! model are expanded by perturbing the wind direction with a set of sample
//! points, the model is run for the unique expanded conditions and the
//! results are combined with the Gaussian weights of the sample points.

use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2, Array3, Float};
use crate::wind_data::WindData;
use anyhow::Result;
use ndarray::s;
use std::collections::HashMap;

/// Shift a value to the range (-180, 180] (Python `floris.utilities.wrap_180`).
fn wrap_180(angle: Float) -> Float {
    ((angle + 180.0) % 360.0) - 180.0
}

/// A hashable key of a row of inputs (the float values are converted to their
/// bit patterns, so that the rows can be compared exactly).
fn row_key(row: &[Float]) -> Vec<u8> {
    let mut key: Vec<u8> = Vec::with_capacity(row.len() * 8);
    for value in row {
        key.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    key
}

/// Round a value to the nearest multiple of `resolution`.
///
/// The ties are rounded to the even value, exactly as `numpy.round` does
/// (Python's `_get_rounded_inputs` uses `np.round(x / resolution) * resolution`).
fn round_to_resolution(value: Float, resolution: Float) -> Float {
    let scaled = value / resolution;
    let fractional = scaled - scaled.trunc();
    if fractional.abs() == 0.5 {
        // A tie: round to the even multiple
        let floor = scaled.floor();
        let even = if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        };
        return even * resolution;
    }
    scaled.round() * resolution
}

/// The resolutions of the expanded inputs (Python's `_get_rounded_inputs`).
///
/// The inputs of every condition are rounded to these resolutions before the
/// conditions are expanded, so that conditions which only differ below the
/// resolution are solved only once.
#[derive(Debug, Clone, Copy)]
pub struct UncertainResolutions {
    /// The wind direction resolution [deg]
    pub wd: Float,
    /// The wind speed resolution [m/s]
    pub ws: Float,
    /// The turbulence intensity resolution
    pub ti: Float,
    /// The yaw angle resolution [deg]
    pub yaw: Float,
    /// The power setpoint resolution [W]
    pub power_setpoint: Float,
    /// The AWC amplitude resolution [deg]
    pub awc_amplitude: Float,
}

impl Default for UncertainResolutions {
    fn default() -> Self {
        Self {
            wd: 1.0,
            ws: 1.0,
            ti: 0.01,
            yaw: 1.0,
            power_setpoint: 100.0,
            awc_amplitude: 0.1,
        }
    }
}

/// A model with an uncertainty in the wind direction.
///
/// The wind direction of every condition is perturbed with the sample points
/// (`[-2σ, -σ, 0, σ, 2σ]` by default) and the resulting turbine values are
/// averaged with the Gaussian weights of the sample points.
#[derive(Debug, Clone)]
pub struct UncertainFlorisModel {
    /// The model with the original (unexpanded) conditions
    pub model: FlorisModel,

    /// The standard deviation of the wind direction distribution [deg]
    pub wd_std: Float,

    /// The wind direction sample points [deg]
    pub wd_sample_points: Vec<Float>,

    /// The normalized Gaussian weights of the sample points
    pub weights: Array1,

    /// The resolutions of the expanded inputs
    pub resolutions: UncertainResolutions,

    /// The model with the unique expanded conditions
    expanded_model: Option<FlorisModel>,

    /// For every expanded condition, the index of the unique condition
    map_to_expanded: Array1,

    /// The number of the original conditions
    n_unexpanded: usize,

    /// The number of the turbines
    n_turbines: usize,

    /// Whether the expanded model has been run
    pub run_completed: bool,

    /// Print the sizes of the expanded and unique cases
    pub verbose: bool,

    /// Fix the yaw angles to the nominal wind direction, so that they are
    /// shifted together with the wind direction (Python
    /// `fix_yaw_to_nominal_direction`).
    pub fix_yaw_to_nominal_direction: bool,

    /// The turbine values of the unique conditions (set by `run` and
    /// `run_parallel`).
    unique_values: Option<Array2>,
}

impl UncertainFlorisModel {
    /// Create an uncertain model.
    ///
    /// * `model` - the model with the original conditions
    /// * `wd_std` - the standard deviation of the wind direction [deg]
    /// * `wd_sample_points` - the wind direction sample points; the default is
    ///   `[-2σ, -σ, 0, σ, 2σ]`
    pub fn new(
        model: FlorisModel,
        wd_std: Float,
        wd_sample_points: Option<Vec<Float>>,
    ) -> Result<Self> {
        if wd_std <= 0.0 {
            anyhow::bail!("wd_std must be strictly greater than 0.");
        }

        let sample_points = wd_sample_points
            .unwrap_or_else(|| vec![-2.0 * wd_std, -wd_std, 0.0, wd_std, 2.0 * wd_std]);
        if sample_points.len() % 2 != 1 {
            anyhow::bail!("wd_sample_points must have an odd length.");
        }
        if sample_points[sample_points.len() / 2] != 0.0 {
            anyhow::bail!("The middle element of wd_sample_points must be 0.");
        }

        // The Gaussian weights of the sample points
        let gaussian: Vec<Float> = sample_points
            .iter()
            .map(|point| (-(point * point) / (2.0 * wd_std * wd_std)).exp())
            .collect();
        let total: Float = gaussian.iter().sum();
        let weights = Array1::from_vec(gaussian.iter().map(|v| v / total).collect::<Vec<Float>>());

        let mut uncertain = Self {
            model,
            wd_std,
            wd_sample_points: sample_points,
            weights,
            resolutions: UncertainResolutions::default(),
            expanded_model: None,
            map_to_expanded: Array1::zeros(0),
            n_unexpanded: 0,
            n_turbines: 0,
            run_completed: false,
            verbose: false,
            fix_yaw_to_nominal_direction: false,
            unique_values: None,
        };
        uncertain.expand()?;
        Ok(uncertain)
    }

    /// Set the resolutions of the expanded inputs and expand the conditions
    /// again.
    pub fn with_resolutions(mut self, resolutions: UncertainResolutions) -> Self {
        self.resolutions = resolutions;
        self
    }

    /// Set the wind speed resolution of the expanded inputs and expand the
    /// conditions again.
    ///
    /// A coarser wind speed resolution reduces the number of conditions
    /// that have to be solved (the approximate model of Python's
    /// `ApproxFlorisModel`).
    pub fn with_ws_resolution(mut self, ws_resolution: Float) -> Self {
        self.resolutions.ws = ws_resolution;
        self
    }

    /// Expand the conditions again (after changing the parameters).
    pub fn expand_public(&mut self) -> Result<()> {
        self.expand()
    }

    /// The model with the unique expanded conditions
    pub fn expanded_model(&self) -> Option<&FlorisModel> {
        self.expanded_model.as_ref()
    }

    /// The number of the original conditions
    pub fn n_unexpanded(&self) -> usize {
        self.n_unexpanded
    }

    /// The number of the unique expanded conditions
    pub fn n_unique(&self) -> usize {
        self.expanded_model
            .as_ref()
            .map(|m| m.core().flow_field.n_findex)
            .unwrap_or(0)
    }

    /// The number of the expanded conditions
    pub fn n_expanded(&self) -> usize {
        self.map_to_expanded.len()
    }

    /// Expand the wind conditions of the model and set them on a copy of the
    /// model.
    fn expand(&mut self) -> Result<()> {
        let flow_field = self.model.core().flow_field.clone();
        let farm = self.model.core().farm.clone();
        let n_unexpanded = flow_field.wind_directions.len();
        let n_turbines = farm.n_turbines();

        // The inputs of every condition: the wind direction, wind speed and
        // turbulence intensity, followed by the yaw angles, the power
        // setpoints and the AWC amplitudes
        let mut inputs: Vec<Vec<Float>> = Vec::with_capacity(n_unexpanded);
        for fi in 0..n_unexpanded {
            let mut row = vec![
                flow_field.wind_directions[fi],
                flow_field.wind_speeds[fi],
                flow_field.turbulence_intensities[fi],
            ];
            for t in 0..n_turbines {
                row.push(farm.yaw_angles[[fi, t]]);
            }
            for t in 0..n_turbines {
                row.push(farm.power_setpoints[[fi, t]]);
            }
            for t in 0..n_turbines {
                row.push(farm.awc_amplitudes[[fi, t]]);
            }
            inputs.push(row);
        }

        // Round the inputs to the resolutions (Python
        // `_get_rounded_inputs`). The resolutions are per column: the wind
        // direction, the wind speed and the turbulence intensity, then the
        // yaw angles, the power setpoints and the AWC amplitudes.
        for row in inputs.iter_mut() {
            row[0] = round_to_resolution(row[0], self.resolutions.wd);
            row[1] = round_to_resolution(row[1], self.resolutions.ws);
            row[2] = round_to_resolution(row[2], self.resolutions.ti);
            for t in 0..n_turbines {
                row[3 + t] = round_to_resolution(row[3 + t], self.resolutions.yaw);
                row[3 + n_turbines + t] =
                    round_to_resolution(row[3 + n_turbines + t], self.resolutions.power_setpoint);
                row[3 + 2 * n_turbines + t] = round_to_resolution(
                    row[3 + 2 * n_turbines + t],
                    self.resolutions.awc_amplitude,
                );
            }
        }

        // Expand the wind directions: every sample point adds a copy of the
        // inputs with a perturbed wind direction
        let n_sample_points = self.wd_sample_points.len();
        let mut expanded: Vec<Vec<Float>> = Vec::with_capacity(n_unexpanded * n_sample_points);
        for point in &self.wd_sample_points {
            for row in &inputs {
                let mut expanded_row = row.clone();
                expanded_row[0] = (expanded_row[0] + point).rem_euclid(360.0);

                // The yaw angles follow the wind direction when they are fixed
                // to the nominal direction (Python
                // `_expand_wind_directions`, uncertain_floris_model.py:892-897)
                if self.fix_yaw_to_nominal_direction {
                    for t in 0..n_turbines {
                        expanded_row[3 + t] = wrap_180(expanded_row[3 + t] + point);
                    }
                }

                expanded.push(expanded_row);
            }
        }

        // The unique expanded conditions
        let mut unique_rows: Vec<Vec<Float>> = Vec::new();
        let mut row_indices: HashMap<Vec<u8>, usize> = HashMap::new();
        let mut map_to_expanded: Vec<usize> = Vec::with_capacity(expanded.len());
        for row in &expanded {
            let key = row_key(row);
            let index = match row_indices.get(&key) {
                Some(index) => *index,
                None => {
                    let index = unique_rows.len();
                    unique_rows.push(row.clone());
                    row_indices.insert(key, index);
                    index
                }
            };
            map_to_expanded.push(index);
        }

        // The model with the unique conditions
        let mut expanded_model = self.model.clone();
        let n_unique = unique_rows.len();
        let mut wind_directions = Array1::zeros(n_unique);
        let mut wind_speeds = Array1::zeros(n_unique);
        let mut turbulence_intensities = Array1::zeros(n_unique);
        let mut yaw_angles = Array2::zeros((n_unique, n_turbines));
        let mut power_setpoints = Array2::zeros((n_unique, n_turbines));
        let mut awc_amplitudes = Array2::zeros((n_unique, n_turbines));
        for (index, row) in unique_rows.iter().enumerate() {
            wind_directions[index] = row[0];
            wind_speeds[index] = row[1];
            turbulence_intensities[index] = row[2];
            for t in 0..n_turbines {
                yaw_angles[[index, t]] = row[3 + t];
                power_setpoints[[index, t]] = row[3 + n_turbines + t];
                awc_amplitudes[[index, t]] = row[3 + 2 * n_turbines + t];
            }
        }

        expanded_model.set_wind_conditions(wind_speeds, wind_directions, turbulence_intensities)?;
        expanded_model.set_operation(
            Some(yaw_angles),
            Some(power_setpoints),
            None,
            Some(awc_amplitudes),
            None,
            None,
        )?;

        self.expanded_model = Some(expanded_model);
        self.map_to_expanded =
            Array1::from_iter(map_to_expanded.iter().map(|index| *index as Float));
        self.n_unexpanded = n_unexpanded;
        self.n_turbines = n_turbines;
        self.run_completed = false;

        if self.verbose {
            println!("Original num rows: {}", self.n_unexpanded);
            println!("Expanded num rows: {}", self.map_to_expanded.len());
            println!("Unique num rows: {}", unique_rows.len());
        }

        Ok(())
    }

    /// Set the layout of the model and expand the conditions again.
    pub fn set_layout(&mut self, layout_x: &Array1, layout_y: &Array1) -> Result<()> {
        self.model.set_layout(layout_x, layout_y)?;
        self.expand()
    }

    /// Set the wind conditions of the model and expand them.
    pub fn set_wind_conditions(
        &mut self,
        wind_speeds: Array1,
        wind_directions: Array1,
        turbulence_intensities: Array1,
    ) -> Result<()> {
        self.model
            .set_wind_conditions(wind_speeds, wind_directions, turbulence_intensities)?;
        self.expand()
    }

    /// Set the wind data of the model and expand the conditions again.
    ///
    /// Python `UncertainFlorisModel.set(wind_data=...)` replaces the wind
    /// conditions of the unexpanded model and expands them.
    pub fn set_wind_data(&mut self, wind_data: &dyn WindData) -> Result<()> {
        self.model.set_wind_data(wind_data)?;
        self.expand()
    }

    /// Set the yaw angles of the model and expand the conditions again.
    ///
    /// The yaw angles are set on the unexpanded model, so all the expanded
    /// wind directions use the same yaw misalignment (unless
    /// `fix_yaw_to_nominal_direction` is set), exactly as Python's
    /// `UncertainFlorisModel.set(yaw_angles=...)`.
    pub fn set_yaw_angles(&mut self, yaw_angles: Array2) -> Result<()> {
        self.model.set_yaw_angles(yaw_angles)?;
        self.expand()
    }

    /// Run the model with the expanded conditions.
    pub fn run(&mut self) -> Result<()> {
        self.expanded_model.as_mut().unwrap().run()?;
        let powers = self.expanded_model.as_ref().unwrap().get_turbine_powers()?;
        self.unique_values = Some(powers);
        self.run_completed = true;
        Ok(())
    }

    /// Run the expanded conditions, splitting the unique conditions over
    /// several chunks (Python's `ParFlorisModel`).
    ///
    /// The results do not depend on the split: every chunk is solved with a
    /// copy of the expanded model.
    pub fn run_parallel(&mut self, n_chunks: usize) -> Result<()> {
        let n_chunks = n_chunks.max(1);
        let expanded = self.expanded_model.as_ref().unwrap();
        let n_unique = expanded.core().flow_field.n_findex;
        let n_turbines = self.n_turbines;

        let wind_directions = expanded.core().flow_field.wind_directions.clone();
        let wind_speeds = expanded.core().flow_field.wind_speeds.clone();
        let turbulence_intensities = expanded.core().flow_field.turbulence_intensities.clone();
        let yaw_angles = expanded.core().farm.yaw_angles.clone();
        let power_setpoints = expanded.core().farm.power_setpoints.clone();

        let mut values = Array2::zeros((n_unique, n_turbines));
        let chunk_size = ((n_unique + n_chunks - 1) / n_chunks).max(1);
        for start in (0..n_unique).step_by(chunk_size) {
            let end = (start + chunk_size).min(n_unique);

            let mut model = expanded.clone();
            model.set_wind_conditions(
                wind_speeds.slice(s![start..end]).to_owned(),
                wind_directions.slice(s![start..end]).to_owned(),
                turbulence_intensities.slice(s![start..end]).to_owned(),
            )?;
            model.set_operation(
                Some(yaw_angles.slice(s![start..end, ..]).to_owned()),
                Some(power_setpoints.slice(s![start..end, ..]).to_owned()),
                None,
                None,
                None,
                None,
            )?;
            model.run()?;

            let chunk_values = model.get_turbine_powers()?;
            for fi in start..end {
                for t in 0..n_turbines {
                    values[[fi, t]] = chunk_values[[fi - start, t]];
                }
            }
        }

        self.unique_values = Some(values);
        self.run_completed = true;
        Ok(())
    }

    /// Run the model without the wakes.
    pub fn run_no_wake(&mut self) -> Result<()> {
        self.expanded_model.as_mut().unwrap().run_no_wake()?;
        self.run_completed = true;
        Ok(())
    }

    /// The turbine powers, averaged over the wind direction distribution.
    pub fn get_turbine_powers(&self) -> Result<Array2> {
        let unique_powers = match &self.unique_values {
            Some(values) => values.clone(),
            None => self.expanded_model.as_ref().unwrap().get_turbine_powers()?,
        };
        Ok(self.map_turbine_values(&unique_powers))
    }

    /// The farm power for every condition, averaged over the wind direction
    /// distribution.
    pub fn get_farm_power(&self) -> Result<Array1> {
        let powers = self.get_turbine_powers()?;
        Ok(powers.sum_axis(ndarray::Axis(1)))
    }

    /// The expected farm power (the conditions are weighted with the
    /// frequencies of the wind data).
    pub fn get_expected_farm_power(&self, freq: Option<Array1>) -> Result<Float> {
        let farm_power = self.get_farm_power()?;
        let frequencies = match freq {
            Some(f) => f,
            None => {
                if let Some(stored) = &self.model.core().wind_data_frequencies {
                    stored.column(0).to_owned()
                } else {
                    Array1::from_elem(farm_power.len(), 1.0 / farm_power.len() as Float)
                }
            }
        };
        if frequencies.len() != farm_power.len() {
            anyhow::bail!(
                "Frequency array length ({}) must match n_findex ({})",
                frequencies.len(),
                farm_power.len()
            );
        }
        let mut expected = 0.0;
        for i in 0..farm_power.len() {
            expected += frequencies[i] * farm_power[i];
        }
        Ok(expected)
    }

    /// The turbine values of the unique conditions, mapped back to the
    /// expanded conditions and averaged with the weights (Python
    /// `map_turbine_values_uncertain`).
    fn map_turbine_values(&self, unique_values: &Array2) -> Array2 {
        let n_turbines = self.n_turbines;
        let mut result = Array2::zeros((self.n_unexpanded, n_turbines));
        for (fi, mut row) in result.rows_mut().into_iter().enumerate() {
            for p in 0..self.weights.len() {
                let expanded_index = p * self.n_unexpanded + fi;
                let unique_index = self.map_to_expanded[expanded_index] as usize;
                for t in 0..n_turbines {
                    row[t] += self.weights[p] * unique_values[[unique_index, t]];
                }
            }
        }
        result
    }

    /// The turbine values of the unique conditions as a 3-D array
    /// (conditions, sample points, turbines) (used for the diagnostics).
    pub fn expanded_turbine_values(&self) -> Result<Array3> {
        let unique_powers = self.expanded_model.as_ref().unwrap().get_turbine_powers()?;
        let n_sample_points = self.weights.len();
        let mut values = Array3::zeros((self.n_unexpanded, n_sample_points, self.n_turbines));
        for p in 0..n_sample_points {
            for fi in 0..self.n_unexpanded {
                let expanded_index = p * self.n_unexpanded + fi;
                let unique_index = self.map_to_expanded[expanded_index] as usize;
                for t in 0..self.n_turbines {
                    values[[fi, p, t]] = unique_powers[[unique_index, t]];
                }
            }
        }
        Ok(values)
    }
}

/// The approximate FLORIS model (Python `ApproxFlorisModel`).
///
/// A special case of the `UncertainFlorisModel` with a single wind
/// direction sample point: no uncertainty is added, but the resolution of
/// the wind speeds (and of the other inputs) is reduced, so that the
/// conditions can be reused and a faster approximate result is computed.
#[derive(Debug, Clone)]
pub struct ApproxFlorisModel {
    /// The underlying model (an uncertain model with a single sample point)
    pub model: UncertainFlorisModel,
}

impl ApproxFlorisModel {
    /// Create an approximate model.
    ///
    /// * `model` - the model with the original conditions
    /// * `ws_resolution` - the wind speed resolution of the expanded
    ///   conditions [m/s]; a coarser resolution reduces the number of
    ///   conditions that have to be solved
    pub fn new(model: FlorisModel, ws_resolution: Float) -> Result<Self> {
        let mut uncertain = UncertainFlorisModel::new(model, 1.0, Some(vec![0.0]))?;
        uncertain.resolutions.ws = ws_resolution;
        uncertain.expand_public()?;
        Ok(Self { model: uncertain })
    }

    /// Set the layout of the model and expand the conditions again.
    pub fn set_layout(&mut self, layout_x: &Array1, layout_y: &Array1) -> Result<()> {
        self.model.set_layout(layout_x, layout_y)
    }

    /// Set the wind data of the model and expand the conditions again.
    pub fn set_wind_data(&mut self, wind_data: &dyn WindData) -> Result<()> {
        self.model.set_wind_data(wind_data)
    }

    /// Set the wind conditions of the model and expand them.
    pub fn set_wind_conditions(
        &mut self,
        wind_speeds: Array1,
        wind_directions: Array1,
        turbulence_intensities: Array1,
    ) -> Result<()> {
        self.model
            .set_wind_conditions(wind_speeds, wind_directions, turbulence_intensities)
    }

    /// Run the model with the expanded conditions.
    pub fn run(&mut self) -> Result<()> {
        self.model.run()
    }

    /// The turbine powers of the conditions.
    pub fn get_turbine_powers(&self) -> Result<Array2> {
        self.model.get_turbine_powers()
    }

    /// The farm power of every condition.
    pub fn get_farm_power(&self) -> Result<Array1> {
        self.model.get_farm_power()
    }

    /// The expected farm power (the conditions are weighted with the
    /// frequencies of the wind data).
    pub fn get_expected_farm_power(&self, freq: Option<Array1>) -> Result<Float> {
        self.model.get_expected_farm_power(freq)
    }
}
