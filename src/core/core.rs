/// Core class - Top-level FLORIS model
use crate::core::farm::Farm;
use crate::core::flow_field::FlowField;
use crate::core::grid::{
    FlowFieldGrid, FlowFieldPlanarGrid, Grid, PointsGrid, TurbineCubatureGrid, TurbineGrid,
};
use crate::core::solver::{cc_solver, empirical_gauss_solver, sequential_solver};
use crate::core::state::State;
use crate::core::wake::WakeModelManager;
use crate::floris_config::{FlorisConfig, LoggingConfig, SolverConfig, SolverType};
use crate::types::{Array1, Array2, Float};
use crate::utilities::load_yaml;
use std::fmt;
use std::path::Path;

/// Core struct - Top-level class that describes a Floris model and initializes the simulation
pub struct Core {
    // Configuration fields (from YAML)
    pub logging: LoggingConfig,
    pub solver: SolverConfig,
    pub wake: WakeModelManager,
    pub farm: Farm,
    pub flow_field: FlowField,

    // Metadata fields (not used in calculations, but kept for completeness)
    pub name: String,
    pub description: String,
    pub floris_version: String,

    // Runtime fields
    pub grid: Option<Box<dyn Grid>>,
    pub state: State,
}

impl fmt::Debug for Core {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Core")
            .field("name", &self.name)
            .field("description", &self.description)
            .field("floris_version", &self.floris_version)
            .field("solver", &self.solver)
            .field("farm", &self.farm)
            .field("flow_field", &self.flow_field)
            .field("wake", &"WakeModelManager")
            .field("grid", &self.grid.as_ref().map(|_| "Box<dyn Grid>"))
            .field("state", &self.state)
            .finish()
    }
}

impl Clone for Core {
    fn clone(&self) -> Self {
        Self {
            logging: self.logging.clone(),
            solver: self.solver.clone(),
            wake: self.wake.clone(),
            farm: self.farm.clone(),
            flow_field: self.flow_field.clone(),
            name: self.name.clone(),
            description: self.description.clone(),
            floris_version: self.floris_version.clone(),
            grid: None, // Grid needs to be recreated
            state: self.state.clone(),
        }
    }
}

impl Core {
    /// Create a Core instance from a YAML configuration file
    ///
    /// # Arguments
    /// * `input_file_path` - Path to the YAML configuration file
    ///
    /// # Returns
    /// Result containing the Core instance or an error
    ///
    /// # Examples
    /// ```no_run
    /// use florus::core::Core;
    /// let core = Core::from_file("default_inputs.yaml").unwrap();
    /// ```
    pub fn from_file<P: AsRef<Path>>(input_file_path: P) -> crate::Result<Self> {
        let config_value = load_yaml(input_file_path)?;

        let config: FlorisConfig = serde_yaml::from_value(config_value)
            .map_err(|e| anyhow::anyhow!("Failed to parse config: {}", e))?;

        Self::from_config(config)
    }

    /// Create a Core instance from default configuration
    pub fn from_default() -> crate::Result<Self> {
        Self::from_file("default_inputs.yaml")
    }

    /// Create a Core instance from a FlorisConfig structure
    ///
    /// This method performs all initialization steps equivalent to Python's __attrs_post_init__
    pub fn from_config(config: FlorisConfig) -> crate::Result<Self> {
        dbg!("Initializing Core from config: {:?}", &config);
        // Create flow field
        let wind_speeds = Array1::from_vec(config.flow_field.wind_speeds.clone());
        let wind_directions = Array1::from_vec(config.flow_field.wind_directions.clone());
        let turbulence_intensities =
            Array1::from_vec(config.flow_field.turbulence_intensities.clone());

        let mut flow_field = FlowField::new(
            wind_speeds,
            wind_directions,
            config.flow_field.wind_veer,
            config.flow_field.wind_shear,
            config.flow_field.air_density,
            turbulence_intensities,
            config.flow_field.reference_wind_height,
        )?;

        // Create farm
        let layout_x = Array1::from_vec(config.farm.layout_x.clone());
        let layout_y = Array1::from_vec(config.farm.layout_y.clone());
        let mut farm = Farm::new(layout_x, layout_y, &config.farm.turbine_type)?;

        // Initialize control arrays
        farm.initialize_control_arrays(flow_field.n_findex);

        // Set reference values
        farm.set_yaw_angles_to_ref_yaw(flow_field.n_findex);
        farm.set_tilt_to_ref_tilt(flow_field.n_findex);
        farm.set_power_setpoints_to_ref_power(flow_field.n_findex);
        farm.set_awc_modes_to_ref_mode(flow_field.n_findex);
        farm.set_awc_amplitudes_to_ref_amp(flow_field.n_findex);
        farm.set_awc_frequencies_to_ref_freq(flow_field.n_findex);

        // Create wake model manager
        let wake = WakeModelManager::from_config(&config.wake)?;

        // Create state
        let state = State::new();

        let mut core = Self {
            logging: config.logging.clone(),
            solver: config.solver.clone(),
            wake,
            farm,
            flow_field,
            name: config.name.clone(),
            description: config.description.unwrap_or_default(),
            floris_version: config.floris_version.clone(),
            grid: None,
            state,
        };

        // Initialize grid based on solver type
        core.initialize_grid()?;

        // Expand farm properties if using turbine grid
        if let Some(ref grid) = core.grid {
            if grid.is_turbine_grid() {
                core.farm
                    .expand_farm_properties(core.flow_field.n_findex, grid.sorted_coord_indices());
            }
        }

        Ok(core)
    }

    /// Initialize the computational grid based on solver configuration
    pub fn initialize_grid(&mut self) -> crate::Result<()> {
        let coords = self.farm.coordinates();
        let rotor_diameters = self.farm.rotor_diameters.clone();
        let wind_directions = self.flow_field.wind_directions.clone();

        match self.solver.solver_type {
            SolverType::TurbineGrid => {
                let grid = TurbineGrid::new(
                    coords,
                    rotor_diameters,
                    wind_directions,
                    self.solver.turbine_grid_points,
                )?;
                self.grid = Some(Box::new(grid));
                if let Some(ref g) = self.grid {
                    self.farm
                        .expand_farm_properties(self.flow_field.n_findex, g.sorted_coord_indices());
                }
            }
            SolverType::TurbineCubatureGrid => {
                let grid = TurbineCubatureGrid::new(
                    coords,
                    rotor_diameters,
                    wind_directions,
                    self.solver.turbine_grid_points,
                )?;
                self.grid = Some(Box::new(grid));
                if let Some(ref g) = self.grid {
                    self.farm
                        .expand_farm_properties(self.flow_field.n_findex, g.sorted_coord_indices());
                }
            }
            SolverType::FlowFieldGrid => {
                // For these solvers, use TurbineGrid by default
                let grid = FlowFieldGrid::new(
                    coords,
                    rotor_diameters,
                    wind_directions,
                    self.solver.turbine_grid_points,
                )?;
                self.grid = Some(Box::new(grid));
            }
            SolverType::FlowFieldPlanarGrid => {
                let grid_resolution = self.solver.grid_resolution.unwrap_or([100, 100]);
                let normal_vector = self
                    .solver
                    .normal_vector
                    .clone()
                    .unwrap_or_else(|| "z".to_string());
                let planar_coordinate = self.solver.planar_coordinate.unwrap_or(90.0);
                let x1_bounds = self.solver.x1_bounds;
                let x2_bounds = self.solver.x2_bounds;

                let grid = FlowFieldPlanarGrid::new(
                    coords,
                    rotor_diameters,
                    wind_directions,
                    grid_resolution,
                    normal_vector,
                    planar_coordinate,
                    x1_bounds,
                    x2_bounds,
                )?;
                self.grid = Some(Box::new(grid));
            }
            _ => {
                anyhow::bail!(
                    "Solver type {:?} not yet fully implemented. Supported types: TurbineGrid, CC, Sequential",
                    self.solver.solver_type
                );
            }
        }

        // Mark as uninitialized since grid changed - will be reinitialized on next run()
        self.state.initialized = false;

        Ok(())
    }

    /// Initialize solution space prior to wake calculations
    ///
    /// This method should be called before performing wake calculations.
    pub fn initialize_domain(&mut self) -> crate::Result<()> {
        if self.grid.is_none() {
            anyhow::bail!("Grid must be initialized before calling initialize_domain()");
        }

        let grid = self.grid.as_ref().unwrap();

        self.flow_field.initialize_flow_field(
            grid.grid_shape(),
            grid.z_sorted(),
            &grid.hub_heights(),
        );

        self.farm.initialize(grid.sorted_coord_indices());

        self.state.initialized = true;

        Ok(())
    }

    /// Perform steady-state wind farm wake calculations
    ///
    /// Note: initialize_domain() must be called before this function.
    pub fn steady_state_atmospheric_condition(&mut self) -> crate::Result<()> {
        if self.grid.is_none() {
            self.initialize_grid()?;
        }

        if !self.state.initialized {
            self.initialize_domain()?;
        }

        let vel_model = self.wake.velocity_model_name.clone();

        // Warning for tilt corrections without vertical wake deflection
        if vel_model != "empirical_gauss" && self.farm.correct_cp_ct_for_tilt.iter().any(|&v| v) {
            log::warn!(
                "The current model does not account for vertical wake deflection due to tilt. \
                 Corrections to power and thrust coefficient can be included, but no vertical \
                 wake deflection will occur."
            );
        }

        // Check for AWC operation model
        let operation_model_awc = self
            .farm
            .turbines
            .iter()
            .any(|t| t.turbine_type.operation_model.model_name() == "awc");

        if vel_model != "empirical_gauss" && operation_model_awc {
            log::warn!(
                "The current model `{}` does not account for additional wake mixing due to \
                 active wake control. Corrections to power and thrust coefficient can be \
                 included, but no enhanced wake recovery will occur.",
                vel_model
            );
        }

        // Get grid reference
        let grid_ref = self.grid.as_ref().unwrap();

        // Try to downcast to specific grid types
        if let Some(turbine_grid) = grid_ref.as_any().downcast_ref::<TurbineGrid>() {
            // Use TurbineGrid for most solvers
            match vel_model.as_str() {
                "cc" => {
                    // For CC solver, we need TurbineCubatureGrid
                    // This is a limitation - CC requires cubature grid
                    anyhow::bail!(
                        "CC solver requires TurbineCubatureGrid, but TurbineGrid was provided"
                    );
                }
                "turbopark" => {
                    log::warn!(
                        "The turbopark model has been superseded by the turboparkgauss model. \
                         We recommend using `velocity_model: turboparkgauss` instead."
                    );
                    // turbopark_solver also needs specific grid type
                    anyhow::bail!("turbopark solver not yet implemented for TurbineGrid");
                }
                "empirical_gauss" => {
                    empirical_gauss_solver(
                        &self.farm,
                        &mut self.flow_field,
                        turbine_grid,
                        &self.wake,
                    )?;
                }
                _ => {
                    // Default to sequential solver
                    sequential_solver(&self.farm, &mut self.flow_field, turbine_grid, &self.wake)?;
                }
            }
        } else if let Some(cubature_grid) = grid_ref.as_any().downcast_ref::<TurbineCubatureGrid>()
        {
            // Use TurbineCubatureGrid for CC and other advanced solvers
            match vel_model.as_str() {
                "cc" => {
                    cc_solver(&self.farm, &mut self.flow_field, cubature_grid, &self.wake)?;
                }
                "turbopark" => {
                    log::warn!(
                        "The turbopark model has been superseded by the turboparkgauss model. \
                         We recommend using `velocity_model: turboparkgauss` instead."
                    );
                    anyhow::bail!("turbopark solver not yet implemented for TurbineCubatureGrid");
                }
                "empirical_gauss" => {
                    empirical_gauss_solver(
                        &self.farm,
                        &mut self.flow_field,
                        cubature_grid,
                        &self.wake,
                    )?;
                }
                _ => {
                    sequential_solver(&self.farm, &mut self.flow_field, cubature_grid, &self.wake)?;
                }
            }
        } else {
            anyhow::bail!("Unsupported grid type for solver");
        }

        self.finalize();

        Ok(())
    }

    /// Finalize the calculation - unsort values to match user-supplied order
    pub fn finalize(&mut self) {
        if let Some(grid) = self.grid.as_ref() {
            if let Some(unsorted_indices) = grid.unsorted_indices() {
                self.flow_field.finalize(unsorted_indices);
                self.farm.finalize(unsorted_indices);
            }
        }
        self.state.converged = true;
    }

    // /// Get turbine powers after calculation
    // pub fn get_turbine_powers(&self) -> crate::Result<Array1> {
    //     if !self.state.converged {
    //         anyhow::bail!("Must run steady_state_atmospheric_condition() before getting powers");
    //     }

    //     let n_turbines = self.farm.n_turbines();
    //     let powers = Array1::zeros(n_turbines);

    //     // TODO: Implement proper power calculation through operation model
    //     // For now, return zeros - the actual power is calculated in the solver

    //     Ok(powers)
    // }

    // /// Get farm total power
    // pub fn get_farm_power(&self) -> crate::Result<Float> {
    //     let powers = self.get_turbine_powers()?;
    //     Ok(powers.sum())
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_core_from_file() {
        // This test requires the default_inputs.yaml file to exist
        // Skip if file doesn't exist
        if Path::new("default_inputs.yaml").exists() {
            let core = Core::from_file("default_inputs.yaml");
            assert!(core.is_ok());
        }
    }

    #[test]
    fn test_core_clone() {
        if Path::new("default_inputs.yaml").exists() {
            let core = Core::from_file("default_inputs.yaml").unwrap();
            let cloned = core.clone();

            assert_eq!(core.name, cloned.name);
            assert_eq!(core.farm.n_turbines(), cloned.farm.n_turbines());
        }
    }
}
