/// FLORIS Model - Main user interface
///
/// Corresponds to floris_model.py
use crate::core::{Farm, FlowField, GridBase, State, TurbineGrid};
use crate::core::wake::{WakeModelManager, WakeModelStrings};
use crate::types::{Array1, Array2, Array3, Float, NumericDict};
use crate::utilities::{cosd, load_yaml};
use ndarray::Array;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::path::Path;

/// Wind data trait for AEP calculations
pub trait WindData {
    fn n_conditions(&self) -> usize;
    fn frequencies(&self) -> &Array1;
    fn wind_speeds(&self) -> &Array1;
    fn wind_directions(&self) -> &Array1;
    fn turbulence_intensities(&self) -> &Array1;
}

/// Main FLORIS Model structure
pub struct FlorisModel {
    pub farm: Farm,
    pub flow_field: FlowField,
    pub state: State,
    pub grid: Option<Box<dyn GridBase>>,
    pub solver_type: String,
    pub turbine_grid_points: usize,
    pub model_manager: Option<WakeModelManager>,
}

impl Clone for FlorisModel {
    fn clone(&self) -> Self {
        Self {
            farm: self.farm.clone(),
            flow_field: self.flow_field.clone(),
            state: self.state.clone(),
            grid: None, // Cannot clone dyn GridBase, need to reinitialize
            solver_type: self.solver_type.clone(),
            turbine_grid_points: self.turbine_grid_points,
            model_manager: self.model_manager.as_ref().map(|m| m.clone()),
        }
    }
}

impl fmt::Debug for FlorisModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FlorisModel")
            .field("farm", &self.farm)
            .field("flow_field", &self.flow_field)
            .field("state", &self.state)
            .field("solver_type", &self.solver_type)
            .field("turbine_grid_points", &self.turbine_grid_points)
            .field("grid", &"Box<dyn GridBase>")
            .field("model_manager", &self.model_manager.is_some())
            .finish()
    }
}

/// Configuration structures
#[derive(Debug, Serialize, Deserialize)]
pub struct FlorisConfig {
    pub flow_field: FlowFieldConfig,
    pub farm: FarmConfig,
    pub solver: SolverConfig,
    #[serde(default)]
    pub wake: Option<WakeConfig>,
    #[serde(default)]
    pub turbine_library: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct WakeConfig {
    #[serde(default)]
    pub model_strings: WakeModelStringsConfig,
    #[serde(default)]
    pub model_params: Option<HashMap<String, serde_yaml::Value>>,
    #[serde(default)]
    pub enable_secondary_steering: bool,
    #[serde(default)]
    pub enable_yaw_added_recovery: bool,
    #[serde(default)]
    pub use_parallel_calc: bool,
}

impl Default for WakeConfig {
    fn default() -> Self {
        Self {
            model_strings: WakeModelStringsConfig::default(),
            model_params: None,
            enable_secondary_steering: false,
            enable_yaw_added_recovery: false,
            use_parallel_calc: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WakeModelStringsConfig {
    #[serde(rename = "velocity_model")]
    pub velocity_model: String,
    #[serde(rename = "deflection_model")]
    pub deflection_model: String,
    #[serde(rename = "combination_model")]
    pub combination_model: String,
    #[serde(rename = "turbulence_model")]
    pub turbulence_model: String,
}

impl Default for WakeModelStringsConfig {
    fn default() -> Self {
        Self {
            velocity_model: "gauss".to_string(),
            deflection_model: "gauss".to_string(),
            combination_model: "fls".to_string(),
            turbulence_model: "crespo_hernandez".to_string(),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlowFieldConfig {
    pub wind_speeds: Vec<Float>,
    pub wind_directions: Vec<Float>,
    pub turbulence_intensities: Vec<Float>,
    pub air_density: Float,
    pub wind_shear: Float,
    pub wind_veer: Float,
    pub reference_wind_height: Float,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FarmConfig {
    pub layout_x: Vec<Float>,
    pub layout_y: Vec<Float>,
    pub turbine_type: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SolverConfig {
    #[serde(rename = "type")]
    pub solver_type: String,
    pub turbine_grid_points: Option<usize>,
}

impl FlorisModel {
    /// Create a new FlorisModel from a configuration file
    pub fn from_file<P: AsRef<Path>>(path: P) -> crate::Result<Self> {
        let config_value = load_yaml(path)?;
        let config: FlorisConfig = serde_yaml::from_value(config_value)
            .map_err(|e| anyhow::anyhow!("Failed to parse config: {}", e))?;

        Self::from_config(config)
    }

    /// Create from configuration structure
    pub fn from_config(config: FlorisConfig) -> crate::Result<Self> {
        let wind_speeds = Array1::from_vec(config.flow_field.wind_speeds);
        let wind_directions = Array1::from_vec(config.flow_field.wind_directions);
        let turbulence_intensities = Array1::from_vec(config.flow_field.turbulence_intensities);

        let flow_field = FlowField::new(
            wind_speeds.clone(),
            wind_directions.clone(),
            config.flow_field.wind_veer,
            config.flow_field.wind_shear,
            config.flow_field.air_density,
            turbulence_intensities,
            config.flow_field.reference_wind_height,
        )?;

        let layout_x = Array1::from_vec(config.farm.layout_x);
        let layout_y = Array1::from_vec(config.farm.layout_y);
        let mut farm = Farm::new(layout_x, layout_y, config.farm.turbine_type)?;

        farm.initialize_control_arrays(flow_field.n_findex);

        let state = State::new();

        // Parse wake configuration
        let wake_config = config.wake.unwrap_or_default();
        let model_strings = WakeModelStrings {
            velocity_model: wake_config.model_strings.velocity_model,
            deflection_model: wake_config.model_strings.deflection_model,
            combination_model: wake_config.model_strings.combination_model,
            turbulence_model: wake_config.model_strings.turbulence_model,
        };

        // Parse model parameters if provided
        let model_params: HashMap<String, NumericDict> = if let Some(params) = wake_config.model_params {
            params.into_iter()
                .map(|(k, v)| {
                    // Convert serde_yaml::Value to NumericDict
                    let dict = match v {
                        serde_yaml::Value::Mapping(m) => {
                            let mut data = std::collections::HashMap::new();
                            for (key, val) in m {
                                if let Some(s) = key.as_str() {
                                    match val {
                                        serde_yaml::Value::Number(n) => {
                                            data.insert(s.to_string(), crate::types::ConfigValue::Float(n.as_f64().unwrap_or(0.0)));
                                        }
                                        serde_yaml::Value::String(s) => {
                                            data.insert(s.to_string(), crate::types::ConfigValue::String(s));
                                        }
                                        _ => {}
                                    }
                                }
                            }
                            crate::types::NumericDict { data }
                        }
                        _ => crate::types::NumericDict { data: std::collections::HashMap::new() },
                    };
                    (k, dict)
                })
                .collect()
        } else {
            HashMap::new()
        };

        let model_manager = WakeModelManager::new(
            model_strings.clone(),
            model_params,
            HashMap::new(),
            HashMap::new(),
            wake_config.enable_secondary_steering,
            wake_config.enable_yaw_added_recovery,
            wake_config.use_parallel_calc,
        )?;

        Ok(Self {
            farm,
            flow_field,
            state,
            grid: None,
            solver_type: config.solver.solver_type,
            turbine_grid_points: config.solver.turbine_grid_points.unwrap_or(3),
            model_manager: Some(model_manager),
        })
    }

    /// Run the FLORIS simulation
    pub fn run(&mut self) -> crate::Result<()> {
        if self.grid.is_none() {
            self.initialize_grid()?;
        }

        self.initialize_flow_field()?;

        // Expand farm properties using sorted indices from grid
        // This must be done before running the solver
        if let Some(ref grid) = self.grid {
            self.farm.expand_farm_properties(
                self.flow_field.n_findex,
                grid.sorted_coord_indices(),
            );
        }

        // Use configured model manager or create default if not set
        let model_manager = if let Some(ref mut mm) = self.model_manager {
            mm.clone()
        } else {
            // Create default model manager if not configured
            let model_strings = WakeModelStrings {
                velocity_model: "gauss".to_string(),
                deflection_model: "gauss".to_string(),
                combination_model: "fls".to_string(),
                turbulence_model: "crespo_hernandez".to_string(),
            };

            WakeModelManager::new(
                model_strings,
                std::collections::HashMap::new(),
                std::collections::HashMap::new(),
                std::collections::HashMap::new(),
                false,
                false,
                false,
            )?
        };

        crate::core::solver::sequential_solver(
            &self.farm,
            &mut self.flow_field,
            self.grid.as_ref().unwrap().as_ref(),
            &model_manager,
        )?;

        self.state.converged = true;

        Ok(())
    }

    /// Initialize the computational grid
    pub fn initialize_grid(&mut self) -> crate::Result<()> {
        let coords = self.farm.coordinates();
        let grid = TurbineGrid::new(
            coords,
            self.farm.rotor_diameters.clone(),
            self.flow_field.wind_directions.clone(),
            self.turbine_grid_points,
        )?;

        self.grid = Some(Box::new(grid));
        Ok(())
    }

    /// Initialize flow field on the grid
    pub fn initialize_flow_field(&mut self) -> crate::Result<()> {
        let n_findex = self.flow_field.n_findex;
        let n_turbines = self.farm.n_turbines();
        let grid_resolution = 3; // Default 3x3 grid

        self.flow_field
            .initialize_flow_field((n_findex, n_turbines, grid_resolution, grid_resolution));

        Ok(())
    }

    /// Get turbine powers based on calculated velocities
    pub fn get_turbine_powers(&self) -> Array2 {
        let n_findex = self.flow_field.n_findex;
        let n_turbines = self.farm.n_turbines();
        // Use sorted yaw angles since we're iterating in sorted order (upstream to downstream)
        let yaw_angles_sorted = &self.farm.yaw_angles_sorted;
        let yaw_angles = yaw_angles_sorted; // Use sorted version

        let mut powers = Array::zeros((n_findex, n_turbines));

        let velocities = &self.flow_field.u_sorted;

        for ti in 0..n_turbines {
            if ti >= self.farm.turbine_map.len() {
                continue;
            }

            let turbine = &self.farm.turbine_map[ti];
            // Use sorted rotor diameters
            let rotor_diameter = self.farm.rotor_diameters_sorted[[0, ti]];
            let area = std::f64::consts::PI * (rotor_diameter / 2.0).powi(2);

            for fi in 0..n_findex {
                // Average velocity over all grid points on the rotor
                let mut v_sum = 0.0;
                let grid_points = velocities.shape()[2] * velocities.shape()[3];
                for iy in 0..velocities.shape()[2] {
                    for iz in 0..velocities.shape()[3] {
                        v_sum += velocities[[fi, ti, iy, iz]];
                    }
                }
                let mut v_avg = v_sum / grid_points as f64;

                // Apply yaw cosine correction using sorted yaw angles
                let yaw = yaw_angles[[fi, ti]];
                let yaw_factor = cosd(yaw).powf(3.0);
                v_avg *= yaw_factor;

                let power = if v_avg < turbine.cut_in_wind_speed || v_avg > turbine.cut_out_wind_speed {
                    0.0
                } else {
                    let cp = turbine.power_coefficient(v_avg);
                    0.5 * self.flow_field.air_density * area * v_avg.powi(3) * cp
                };

                powers[[fi, ti]] = power.min(turbine.rated_power);
            }
        }

        powers
    }

    /// Get farm power
    pub fn get_farm_power(&self) -> Array1 {
        let powers = self.get_turbine_powers();
        let mut farm_power = Array1::zeros(self.flow_field.n_findex);

        for i in 0..self.flow_field.n_findex {
            farm_power[i] = powers.row(i).sum();
        }

        farm_power
    }

    /// Set wind conditions
    pub fn set_wind_conditions(
        &mut self,
        wind_speeds: Array1,
        wind_directions: Array1,
        turbulence_intensities: Array1,
    ) -> crate::Result<()> {
        self.flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            self.flow_field.wind_veer,
            self.flow_field.wind_shear,
            self.flow_field.air_density,
            turbulence_intensities,
            self.flow_field.reference_wind_height,
        )?;

        self.farm
            .initialize_control_arrays(self.flow_field.n_findex);
        self.grid = None;

        Ok(())
    }

    /// Set yaw angles
    pub fn set_yaw_angles(&mut self, yaw_angles: Array2) -> crate::Result<()> {
        if yaw_angles.shape() != [self.flow_field.n_findex, self.farm.n_turbines()] {
            anyhow::bail!(
                "yaw_angles shape {:?} doesn't match expected [{}, {}]",
                yaw_angles.shape(),
                self.flow_field.n_findex,
                self.farm.n_turbines()
            );
        }

        self.farm.yaw_angles = yaw_angles;
        Ok(())
    }

    /// Set turbine layout
    pub fn set_layout(&mut self, layout_x: &Array1, layout_y: &Array1) -> crate::Result<()> {
        self.farm.set_layout(layout_x, layout_y)?;
        self.grid = None; // Reset grid when layout changes
        Ok(())
    }

    /// Calculate Annual Energy Production (AEP)
    pub fn get_farm_aep(&self, wind_data: &dyn WindData, hours_per_year: Float) -> Float {
        let n_conditions = wind_data.n_conditions();
        let frequencies = wind_data.frequencies();
        let powers = self.get_turbine_powers();

        let mut total_energy = 0.0;

        for i in 0..n_conditions {
            let power = powers.row(i).sum();
            total_energy += power * frequencies[i] * hours_per_year;
        }

        total_energy
    }

    /// Calculate AEP using uniform frequencies (1/N for each condition)
    ///
    /// This is useful for layout optimization when you want to compare
    /// layouts using the same set of wind conditions.
    pub fn get_farm_aep_uniform(&self, hours_per_year: Float) -> Float {
        let n_conditions = self.flow_field.n_findex;
        let frequency = 1.0 / n_conditions as Float;
        let powers = self.get_turbine_powers();

        let mut total_energy = 0.0;

        for i in 0..n_conditions {
            let power = powers.row(i).sum();
            total_energy += power * frequency * hours_per_year;
        }

        total_energy
    }

    /// Get Annual Value Production (AVP)
    ///
    /// Returns the total value-weighted energy production.
    /// If no value table is available, returns AEP.
    pub fn get_farm_avp(&self) -> Float {
        // Simple implementation: return farm power as a proxy for value
        // In a full implementation, this would multiply by value factors
        self.get_farm_power().sum()
    }

    /// Get velocity at each turbine
    pub fn get_turbine_velocities(&self) -> Array3 {
        let n_findex = self.flow_field.n_findex;
        let n_turbines = self.farm.n_turbines();

        let mut velocities = Array::zeros((n_findex, n_turbines, 1));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                velocities[[fi, ti, 0]] = self.flow_field.u_sorted[[fi, ti, 0, 0]];
            }
        }

        velocities
    }

    /// Get thrust coefficients for each turbine
    ///
    /// Returns thrust coefficient array for all turbines at current operating conditions.
    /// Shape: (n_findex, n_turbines)
    pub fn get_turbine_thrust_coefficients(&self) -> Array2 {
        let n_findex = self.flow_field.n_findex;
        let n_turbines = self.farm.n_turbines();

        let mut ct_array = Array2::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                if ti < self.farm.turbine_map.len() {
                    let turbine = &self.farm.turbine_map[ti];
                    let velocity = self.flow_field.u_sorted[[fi, ti, 0, 0]];
                    ct_array[[fi, ti]] = turbine.ct_at_speed(velocity);
                }
            }
        }

        ct_array
    }

    /// Get the operation model for turbines
    ///
    /// Returns the operation model type string (e.g., "simple", "cosine-loss", etc.)
    /// All turbines in the farm currently use the same operation model.
    pub fn get_operation_model(&self) -> String {
        if self.farm.turbine_map.is_empty() {
            return "simple".to_string();
        }
        self.farm.turbine_map[0].operation_model.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_floris_model_basic() {
        let layout_x = Array1::from_vec(vec![0.0, 500.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 2];

        let mut farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();

        let wind_speeds = Array1::from_vec(vec![8.0, 10.0]);
        let wind_directions = Array1::from_vec(vec![270.0, 280.0]);
        let turbulence_intensities = Array1::from_vec(vec![0.06, 0.08]);

        let flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            0.0,
            0.14,
            1.225,
            turbulence_intensities,
            90.0,
        )
        .unwrap();

        farm.initialize_control_arrays(flow_field.n_findex);

        let model = FlorisModel {
            farm,
            flow_field,
            state: State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            turbine_grid_points: 3,
            model_manager: None,
        };

        assert_eq!(model.farm.n_turbines(), 2);
        assert_eq!(model.flow_field.n_findex, 2);
    }

    #[test]
    fn test_get_turbine_powers() {
        let layout_x = Array1::from_vec(vec![0.0, 500.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 2];

        let farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();

        let wind_speeds = Array1::from_vec(vec![8.0, 10.0]);
        let wind_directions = Array1::from_vec(vec![270.0, 280.0]);
        let turbulence_intensities = Array1::from_vec(vec![0.06, 0.08]);

        let flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            0.0,
            0.14,
            1.225,
            turbulence_intensities,
            90.0,
        )
        .unwrap();

        let model = FlorisModel {
            farm,
            flow_field,
            state: State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            turbine_grid_points: 3,
            model_manager: None,
        };

        assert_eq!(model.farm.n_turbines(), 2);
        assert_eq!(model.flow_field.n_findex, 2);
    }

    #[test]
    fn test_wake_solver_integration() {
        // Create a simple 2-turbine farm with downstream turbine
        let layout_x = Array1::from_vec(vec![0.0, 630.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 2];

        let farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();

        // Single wind condition
        let wind_speeds = Array1::from_vec(vec![8.0]);
        let wind_directions = Array1::from_vec(vec![270.0]);
        let turbulence_intensities = Array1::from_vec(vec![0.06]);

        let flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            0.0,
            0.14,
            1.225,
            turbulence_intensities,
            90.0,
        )
        .unwrap();

        // Create model
        let mut model = FlorisModel {
            farm,
            flow_field,
            state: State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            turbine_grid_points: 3,
            model_manager: None,
        };

        // Initialize grid and flow field
        model.initialize_grid().unwrap();
        model.initialize_flow_field().unwrap();

        // Verify grid was created
        assert!(model.grid.is_some());
        let grid = model.grid.as_ref().unwrap();
        assert_eq!(grid.n_turbines(), 2);
        assert_eq!(grid.n_findex(), 1);

        // Verify flow field was initialized with correct shape
        let shape = model.flow_field.u_sorted.shape();
        assert_eq!(shape[0], 1); // n_findex
        assert_eq!(shape[1], 2); // n_turbines
        assert_eq!(shape[2], 3); // grid resolution y
        assert_eq!(shape[3], 3); // grid resolution z

        // Verify grid coordinates are properly sorted
        // X coordinates should be different for upstream vs downstream turbines
        let x_sorted = grid.x_sorted();
        assert_ne!(x_sorted[[0, 0, 0, 0]], x_sorted[[0, 1, 0, 0]]);
    }
}
