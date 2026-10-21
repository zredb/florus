//! The uncertain FLORIS model: the mean power of a farm for a distribution
//! of wind direction deviations.
//!
//! Python `floris/uncertain_floris_model.py`: the wind conditions of the
//! model are expanded by perturbing the wind direction with a set of sample
//! points, the model is run for the unique expanded conditions and the
//! results are combined with the Gaussian weights of the sample points.

use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2, Array3, Float};
use anyhow::Result;
use std::collections::HashMap;

/// A hashable key of a row of inputs (the float values are converted to their
/// bit patterns, so that the rows can be compared exactly).
fn row_key(row: &[Float]) -> Vec<u8> {
    let mut key: Vec<u8> = Vec::with_capacity(row.len() * 8);
    for value in row {
        key.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    key
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
            expanded_model: None,
            map_to_expanded: Array1::zeros(0),
            n_unexpanded: 0,
            n_turbines: 0,
            run_completed: false,
            verbose: false,
        };
        uncertain.expand()?;
        Ok(uncertain)
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
        // `_get_rounded_inputs`)
        let resolutions = [
            1.0,   // wind direction [deg]
            1.0,   // wind speed [m/s]
            0.01,  // turbulence intensity
            1.0,   // yaw angle [deg]
            100.0, // power setpoint [W]
            0.1,   // AWC amplitude [deg]
        ];
        for row in inputs.iter_mut() {
            for (i, resolution) in resolutions.iter().enumerate() {
                if i < row.len() {
                    row[i] = (row[i] / resolution).round() * resolution;
                }
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

    /// Run the model with the expanded conditions.
    pub fn run(&mut self) -> Result<()> {
        self.expanded_model.as_mut().unwrap().run()?;
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
        let unique_powers = self.expanded_model.as_ref().unwrap().get_turbine_powers()?;
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
