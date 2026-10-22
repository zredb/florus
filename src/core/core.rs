/// Core class - Top-level FLORIS model
use crate::core::farm::Farm;
use crate::core::flow_field::FlowField;
use crate::core::grid::{
    FlowFieldGrid, FlowFieldPlanarGrid, Grid, PointsGrid, TurbineCubatureGrid, TurbineGrid,
};
use crate::core::solver::{
    cc_solver, empirical_gauss_solver, full_flow_sequential_solver, sequential_solver,
    turbopark_solver,
};
use crate::core::state::State;
use crate::core::wake::WakeModelManager;
use crate::floris_config::{FlorisConfig, LoggingConfig, SolverConfig, SolverType};
use crate::types::{Array1, Array2, Float};
use crate::utilities::load_yaml;
use crate::wind_data::traits::WindData;
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

    // Wind data metadata for reshaping (WindRose information)
    pub wind_data_info: Option<WindDataInfo>,

    /// The multi-dimensional conditions of the input file (Python
    /// `flow_field.multidim_conditions`). They select the power/thrust table
    /// of the multi-dimensional turbines.
    pub multidim_conditions: Option<crate::floris_config::MultiDimConditions>,

    /// The wake configuration of the model. It can be changed at runtime
    /// (Python rebuilds the model from a modified dictionary), which is
    /// applied with `FlorisModel::set_wake_config`.
    pub wake_config: crate::floris_config::WakeConfig,

    /// The frequencies of the wind data (Python `wind_data.unpack_freq()`).
    ///
    /// They are used by the expected-power and AEP queries when no explicit
    /// frequencies are given: a single column weights every wind condition,
    /// while a `(n_findex, n_turbines)` table (the wind resource grid) weights
    /// the power of every turbine separately.
    pub wind_data_frequencies: Option<Array2>,

    /// The values (e.g. the electricity price) of the wind data
    /// (Python wind_data.unpack_value()).
    pub wind_data_values: Option<Array2>,

    /// The wind data of the model (Python `FlorisModel._wind_data`).
    ///
    /// It is kept so that it can be updated when the layout of the farm
    /// changes: the wind resource grid generates the wind roses of the
    /// turbines from their positions (Python `floris_model.py:204-207`).
    pub wind_data: Option<Box<dyn WindData>>,
}

/// Information about wind data structure for reshaping results
#[derive(Debug, Clone)]
pub struct WindDataInfo {
    pub wind_directions: Vec<f64>,
    pub wind_speeds: Vec<f64>,
    pub n_wind_directions: usize,
    pub n_wind_speeds: usize,
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
            wind_data_info: self.wind_data_info.clone(),
            multidim_conditions: self.multidim_conditions.clone(),
            wind_data_frequencies: self.wind_data_frequencies.clone(),
            wind_data_values: self.wind_data_values.clone(),
            wind_data: self.wind_data.as_ref().map(|wd| wd.clone_box()),
            wake_config: self.wake_config.clone(),
        }
    }
}

/// One velocity deficit profile (Python: the DataFrame returned by
/// Core.solve_for_velocity_deficit_profiles).
#[derive(Debug, Clone)]
pub struct VelocityDeficitProfile {
    /// The inertial x coordinates of the sample points
    pub x: Array1,
    /// The inertial y coordinates of the sample points
    pub y: Array1,
    /// The heights of the sample points
    pub z: Array1,
    /// The streamwise coordinate, normalized by the reference rotor diameter
    pub x1: Array1,
    /// The cross-stream coordinate, normalized by the reference rotor diameter
    pub x2: Array1,
    /// The vertical coordinate, normalized by the reference rotor diameter
    pub x3: Array1,
    /// The velocity deficit of every sample point
    pub velocity_deficit: Array1,
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
            wind_data_info: None, // Default to None for YAML-based initialization
            multidim_conditions: config.flow_field.multidim_conditions.clone(),
            wind_data_frequencies: None,
            wind_data_values: None,
            wind_data: None,
            wake_config: config.wake.clone(),
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

        // Get hub heights from farm turbines
        use ndarray::Array;
        let n_turbines = self.farm.n_turbines();
        let mut hub_heights = Array::zeros(n_turbines);
        for i in 0..n_turbines {
            hub_heights[i] = self.farm.turbines[i].turbine_type.hub_height;
        }

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
                    hub_heights,
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

        self.flow_field.initialize_flow_field_with_coordinates(
            grid.grid_shape(),
            grid.z_sorted(),
            &grid.hub_heights(),
            Some(grid.x_sorted_inertial_frame()),
            Some(grid.y_sorted_inertial_frame()),
        );

        self.farm.initialize(grid.sorted_coord_indices());

        // Apply the multi-dimensional conditions (Python
        // `flow_field.multidim_conditions`): they select the power/thrust
        // table of the multi-dimensional turbines.
        if let Some(conditions) = &self.multidim_conditions {
            let n_findex = self.flow_field.n_findex;
            self.farm.multidim_conditions = Some(conditions.expand(n_findex)?);
        }

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
                    // Per Python, the cc model runs through cc_solver with the
                    // turbine grid (the cumulative state architecture)
                    cc_solver(&self.farm, &mut self.flow_field, turbine_grid, &self.wake)?;
                }
                "turbopark" => {
                    log::warn!(
                        "The turbopark model has been superseded by the turboparkgauss model. \
                         We recommend using `velocity_model: turboparkgauss` instead."
                    );
                    turbopark_solver(&self.farm, &mut self.flow_field, turbine_grid, &self.wake)?;
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
                    turbopark_solver(&self.farm, &mut self.flow_field, cubature_grid, &self.wake)?;
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
        } else if let Some(planar_grid) = grid_ref.as_any().downcast_ref::<FlowFieldPlanarGrid>() {
            // Use FlowFieldPlanarGrid for visualization
            // This is similar to solve_for_viz in Python
            // For planar grids, we use the generic solver dispatch
            match vel_model.as_str() {
                "empirical_gauss" => {
                    empirical_gauss_solver(
                        &self.farm,
                        &mut self.flow_field,
                        planar_grid,
                        &self.wake,
                    )?;
                }
                _ => {
                    // For other models, use sequential solver with dyn Grid
                    sequential_solver(&self.farm, &mut self.flow_field, planar_grid, &self.wake)?;
                }
            }
        } else {
            anyhow::bail!("Unsupported grid type for solver");
        }

        self.finalize();

        Ok(())
    }

    /// Solve for visualization with FlowFieldPlanarGrid
    ///
    /// This method is specifically designed for visualization purposes.
    /// It should be used with FlowFieldPlanarGrid after setting up the grid.
    pub fn solve_for_viz(&mut self) -> crate::Result<()> {
        if self.grid.is_none() {
            anyhow::bail!("Grid must be initialized before calling solve_for_viz()");
        }

        let grid = self.grid.as_ref().unwrap();

        self.flow_field.initialize_flow_field_with_coordinates(
            grid.grid_shape(),
            grid.z_sorted(),
            &grid.hub_heights(),
            Some(grid.x_sorted_inertial_frame()),
            Some(grid.y_sorted_inertial_frame()),
        );

        let grid_ref = self.grid.as_ref().unwrap();
        if let Some(planar_grid) = grid_ref.as_any().downcast_ref::<FlowFieldPlanarGrid>() {
            full_flow_sequential_solver(&self.farm, &mut self.flow_field, planar_grid, &self.wake)?;
        } else {
            anyhow::bail!("solve_for_viz requires a FlowFieldPlanarGrid");
        }

        self.finalize();
        self.state.initialized = true;

        Ok(())
    }

    /// Solve the wake at arbitrary points for ALL wind conditions, returning
    /// `(velocities, turbulence_intensities)` each shaped
    /// `(n_findex, n_points)` — the `[:, :, 0, 0]` slices of the solved
    /// fields (Python `core.solve_for_points`, core.py:230-266, and the
    /// fields `sample_flow_at_points` / `sample_ti_at_points` consume).
    ///
    /// Differences from Python: the solve runs on a scratch copy of the flow
    /// field and the model state is left intact (Python re-initializes
    /// `self.flow_field` on the points grid, clobbering the main solve
    /// state) — numerically identical, strictly safer. The `cc` and
    /// `empirical_gauss` models are rejected until their full-flow solver
    /// variants are ported; `turbopark` matches Python's explicit refusal.
    /// Extract velocity deficit profiles at a set of downstream distances.
    ///
    /// Python `core.solve_for_velocity_deficit_profiles`
    /// (floris/core/core.py): for every downstream distance a line of sample
    /// points is created along the cross-stream (or vertical) direction, the
    /// points are transformed into the inertial frame and the velocity is
    /// sampled. The deficit is `(homogeneous_wind_speed - u) /
    /// homogeneous_wind_speed`.
    pub fn solve_for_velocity_deficit_profiles(
        &self,
        direction: &str,
        downstream_dists: &[Float],
        profile_range: &[Float],
        resolution: usize,
        homogeneous_wind_speed: Float,
        ref_rotor_diameter: Float,
        x_start: Float,
        y_start: Float,
        reference_height: Float,
    ) -> crate::Result<Vec<VelocityDeficitProfile>> {
        if direction != "cross-stream" && direction != "vertical" {
            anyhow::bail!("`direction` must be either `cross-stream` or `vertical`.");
        }
        if resolution < 2 {
            anyhow::bail!("`resolution` must be at least 2.");
        }

        // The rotation of the wind-direction relative frame into the inertial
        // frame (Python utilities.py:reverse_rotate_coordinates_rel_west):
        // the angle is the deviation from west.
        let wind_direction = self.flow_field.wind_directions[0];
        let angle = (-(wind_direction - 270.0)).to_radians();
        let (sin_angle, cos_angle) = angle.sin_cos();

        let n_lines = downstream_dists.len();
        let mut points_x: Vec<Float> = Vec::with_capacity(n_lines * resolution);
        let mut points_y: Vec<Float> = Vec::with_capacity(n_lines * resolution);
        let mut points_z: Vec<Float> = Vec::with_capacity(n_lines * resolution);
        let mut x1_values: Vec<Float> = Vec::with_capacity(n_lines * resolution);
        let mut x2_values: Vec<Float> = Vec::with_capacity(n_lines * resolution);
        let mut x3_values: Vec<Float> = Vec::with_capacity(n_lines * resolution);

        for downstream in downstream_dists {
            for j in 0..resolution {
                let single_line = profile_range[0]
                    + (profile_range[1] - profile_range[0]) * j as Float
                        / (resolution - 1) as Float;
                let (x1, x2, x3) = if direction == "cross-stream" {
                    (*downstream, single_line, 0.0)
                } else {
                    (*downstream, 0.0, single_line)
                };

                // Rotate about the origin and translate
                points_x.push(x1 * cos_angle - x2 * sin_angle + x_start);
                points_y.push(x1 * sin_angle + x2 * cos_angle + y_start);
                points_z.push(x3 + reference_height);

                x1_values.push(x1 / ref_rotor_diameter);
                x2_values.push(x2 / ref_rotor_diameter);
                x3_values.push(x3 / ref_rotor_diameter);
            }
        }

        let (velocities, _) = self.solve_for_points(&points_x, &points_y, &points_z)?;

        let mut profiles = Vec::with_capacity(n_lines);
        for i in 0..n_lines {
            let start = i * resolution;
            let end = start + resolution;
            let deficit: Array1 = (start..end)
                .map(|k| (homogeneous_wind_speed - velocities[[0, k]]) / homogeneous_wind_speed)
                .collect();
            profiles.push(VelocityDeficitProfile {
                x: Array1::from_vec(points_x[start..end].to_vec()),
                y: Array1::from_vec(points_y[start..end].to_vec()),
                z: Array1::from_vec(points_z[start..end].to_vec()),
                x1: Array1::from_vec(x1_values[start..end].to_vec()),
                x2: Array1::from_vec(x2_values[start..end].to_vec()),
                x3: Array1::from_vec(x3_values[start..end].to_vec()),
                velocity_deficit: deficit,
            });
        }

        Ok(profiles)
    }

    pub fn solve_for_points(
        &self,
        points_x: &[f64],
        points_y: &[f64],
        points_z: &[f64],
    ) -> crate::Result<(Array2, Array2)> {
        let vel_model = self.wake.velocity_model_name.clone();

        // Python core.py:254-258
        if vel_model == "turbopark" {
            anyhow::bail!(
                "solve_for_points is not available for the legacy 'turbopark' model. \
                 However, it is available for 'turboparkgauss'."
            );
        }
        // Python routes these through full_flow_cc_solver /
        // full_flow_empirical_gauss_solver (core.py:259-262); those variants
        // are not yet ported to FLORUS.
        if vel_model == "cc" || vel_model == "empirical_gauss" {
            anyhow::bail!(
                "solve_for_points is not yet supported for the '{}' model in FLORUS \
                 (full_flow solver not ported); Python FLORIS supports it.",
                vel_model
            );
        }

        // Rotation center: the bounding-box center of the farm — the same
        // center the internal TurbineGrid (built by full_flow_sequential_solver
        // from farm.coordinates()) rotates about, and the same one the main
        // turbine grids store. Passing it explicitly keeps the points in the
        // farm's rotated frame (Python passes the turbine grid's centers,
        // core.py:246-247). No prior run() is required, matching Python,
        // because the solve happens on a scratch flow field.
        let coords = self.farm.coordinates();
        let x_min = coords
            .column(0)
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let x_max = coords
            .column(0)
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let y_min = coords
            .column(1)
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let y_max = coords
            .column(1)
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let x_center = Some((x_min + x_max) / 2.0);
        let y_center = Some((y_min + y_max) / 2.0);

        let points_grid = PointsGrid::new(
            Array1::from_iter(points_x.iter().copied()),
            Array1::from_iter(points_y.iter().copied()),
            Array1::from_iter(points_z.iter().copied()),
            self.flow_field.wind_directions.clone(),
            x_center,
            y_center,
        )?;

        // Solve on a scratch copy; the main solve state stays intact.
        let mut points_flow_field = self.flow_field.clone();
        // The wind shear is referenced to the hub height of the farm, not to
        // the heights of the sampling points (Python uses the turbine grid's
        // hub heights, core.py:246-247).
        let farm_hub_heights = self.farm.hub_heights().clone();
        points_flow_field.initialize_flow_field_with_coordinates(
            points_grid.grid_shape(),
            points_grid.z_sorted(),
            &farm_hub_heights,
            Some(points_grid.x_sorted_inertial_frame()),
            Some(points_grid.y_sorted_inertial_frame()),
        );

        full_flow_sequential_solver(&self.farm, &mut points_flow_field, &points_grid, &self.wake)?;

        let n_findex = points_flow_field.u_sorted.shape()[0];
        let n_points = points_grid.n_turbines();
        let mut velocities = Array2::zeros((n_findex, n_points));
        let mut turbulence_intensities = Array2::zeros((n_findex, n_points));
        for fi in 0..n_findex {
            for pi in 0..n_points {
                velocities[[fi, pi]] = points_flow_field.u_sorted[[fi, pi, 0, 0]];
                turbulence_intensities[[fi, pi]] =
                    points_flow_field.turbulence_intensity_field_sorted[[fi, pi, 0, 0]];
            }
        }
        Ok((velocities, turbulence_intensities))
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
