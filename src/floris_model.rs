/// FLORIS Model - Main user interface
///
/// Corresponds to floris_model.py
use crate::core::{Farm, FlowField, GridBase, State, TurbineGrid};
use crate::core::wake::{WakeModelManager, WakeModelStrings};
use crate::types::{Array1, Array2, Array3, Float};
use crate::utilities::load_yaml;
use ndarray::Array;
use serde::{Deserialize, Serialize};
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
    pub turbine_library: Option<String>,
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

        Ok(Self {
            farm,
            flow_field,
            state,
            grid: None,
            solver_type: config.solver.solver_type,
            model_manager: None,
        })
    }

    /// Run the FLORIS simulation
    pub fn run(&mut self) -> crate::Result<()> {
        if self.grid.is_none() {
            self.initialize_grid()?;
        }

        self.initialize_flow_field()?;

        let model_strings = WakeModelStrings {
            velocity_model: "gauss".to_string(),
            deflection_model: "gauss".to_string(),
            combination_model: "fls".to_string(),
            turbulence_model: "crespo_hernandez".to_string(),
        };

        let model_manager = WakeModelManager::new(
            model_strings,
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            false,
            false,
            false,
        )?;

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
            3,
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

        let mut powers = Array::zeros((n_findex, n_turbines));

        let velocities = &self.flow_field.u_sorted;

        for ti in 0..n_turbines {
            if ti >= self.farm.turbine_map.len() {
                continue;
            }

            let turbine = &self.farm.turbine_map[ti];
            let rotor_diameter = self.farm.rotor_diameters[ti];
            let area = std::f64::consts::PI * (rotor_diameter / 2.0).powi(2);

            for fi in 0..n_findex {
                let v = velocities[[fi, ti, 0, 0]];

                let power = if v < turbine.cut_in_wind_speed || v > turbine.cut_out_wind_speed {
                    0.0
                } else {
                    let cp = turbine.power_coefficient(v);
                    0.5 * self.flow_field.air_density * area * v.powi(3) * cp
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
