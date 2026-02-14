use crate::core::wake::WakeModelManager;
use crate::core::{Farm, FlowField, GridBase, State, TurbineGrid};
use crate::floris_config::{FlorisConfig, SolverConfig};
use crate::types::{Array1, Array2, Array3, Float};
use crate::utilities::{cosd, load_yaml};
use crate::wind_data::WindData;
use ndarray::Array;
use std::fmt;
use std::path::Path;

/// Main FLORIS Model structure
pub struct FlorisModel {
    pub farm: Farm,
    pub flow_field: FlowField,
    pub state: State,
    pub grid: Option<Box<dyn GridBase>>,
    pub solver: SolverConfig,
    pub model_manager: Option<WakeModelManager>,
}

impl Clone for FlorisModel {
    fn clone(&self) -> Self {
        Self {
            farm: self.farm.clone(),
            flow_field: self.flow_field.clone(),
            state: self.state.clone(),
            grid: None, // Cannot clone dyn GridBase, need to reinitialize
            solver: self.solver.clone(),
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
            .field("solver", &self.solver)
            .field("grid", &"Box<dyn GridBase>")
            .field("model_manager", &self.model_manager.is_some())
            .finish()
    }
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
        let wake_config = &config.wake;

        let model_manager = WakeModelManager::from_config(wake_config)?;

        Ok(Self {
            farm,
            flow_field,
            state,
            grid: None,
            solver: config.solver,
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
            self.farm
                .expand_farm_properties(self.flow_field.n_findex, grid.sorted_coord_indices());
        }

        // Use configured model manager or create default if not set
        let model_manager = if let Some(ref mut mm) = self.model_manager {
            mm.clone()
        } else {
            // Create default model manager if not configured
            WakeModelManager::default_gauss()?
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
            self.solver.turbine_grid_points,
        )?;

        self.grid = Some(Box::new(grid));
        Ok(())
    }

    /// Initialize flow field on the grid
    pub fn initialize_flow_field(&mut self) -> crate::Result<()> {
        let n_findex = self.flow_field.n_findex;
        let n_turbines = self.farm.n_turbines();
        let grid_resolution = 3; // Default 3x3 grid

        // Get z coordinates and hub heights from grid
        if let Some(ref grid) = self.grid {
            let z_grid = grid.z_sorted();
            let hub_heights = self.farm.hub_heights.clone();
            
            self.flow_field.initialize_flow_field(
                (n_findex, n_turbines, grid_resolution, grid_resolution),
                z_grid,
                &hub_heights,
            );
        } else {
            // Fallback if grid not initialized
            self.flow_field.initialize_flow_field(
                (n_findex, n_turbines, grid_resolution, grid_resolution),
                &ndarray::Array::zeros((n_findex, n_turbines, grid_resolution, grid_resolution)),
                &self.farm.hub_heights,
            );
        }

        Ok(())
    }

    /// Get turbine powers based on calculated velocities
    /// 
    /// This implementation matches Python FLORIS's CosineLossTurbine.power() method:
    /// 1. Uses cubic-mean for rotor average velocity
    /// 2. Applies air density correction
    /// 3. Applies yaw cosine correction
    /// 4. Uses power curve interpolation directly
    pub fn get_turbine_powers(&self) -> Array2 {
        let n_findex = self.flow_field.n_findex;
        let n_turbines = self.farm.n_turbines();
        let yaw_angles = &self.farm.yaw_angles_sorted;
        let tilt_angles = &self.farm.tilt_angles_sorted;

        let mut powers = Array::zeros((n_findex, n_turbines));
        let velocities = &self.flow_field.u_sorted;

        // Default parameters matching Python FLORIS
        let ref_air_density = 1.225;  // Standard air density
        let cosine_loss_exponent_yaw = 1.88;  // Default from Python FLORIS
        let cosine_loss_exponent_tilt = 1.88;
        let ref_tilt = 5.0;  // Default reference tilt angle for nrel_5MW

        for ti in 0..n_turbines {
            if ti >= self.farm.turbine_map.len() {
                continue;
            }

            let turbine = &self.farm.turbine_map[ti];
            let power_curve = turbine.turbine_type.power_curve();
            let rated_power = power_curve.values.iter().cloned().fold(0.0, f64::max) * 1000.0;

            for fi in 0..n_findex {
                // Step 1: Calculate rotor average velocity using cubic mean
                // cubic_mean = cbrt(mean(v^3))
                let mut v_cubed_sum = 0.0;
                let grid_points = velocities.shape()[2] * velocities.shape()[3];
                for iy in 0..velocities.shape()[2] {
                    for iz in 0..velocities.shape()[3] {
                        v_cubed_sum += velocities[[fi, ti, iy, iz]].powi(3);
                    }
                }
                let rotor_avg_velocity = (v_cubed_sum / grid_points as f64).powf(1.0 / 3.0);

                // Step 2: Apply air density correction
                // rotor_effective_velocity = (air_density/ref_air_density)^(1/3) * rotor_avg_velocity
                let density_factor = (self.flow_field.air_density / ref_air_density).powf(1.0 / 3.0);
                let mut rotor_effective_velocity = rotor_avg_velocity * density_factor;

                // Step 3: Apply yaw cosine correction
                // rotor_effective_velocity *= cosd(yaw)^(cosine_loss_exponent_yaw/3)
                let yaw = yaw_angles[[fi, ti]];
                let yaw_correction = cosd(yaw).powf(cosine_loss_exponent_yaw / 3.0);
                rotor_effective_velocity *= yaw_correction;

                // Step 4: Apply tilt cosine correction
                let tilt = tilt_angles[[fi, ti]];
                let tilt_correction = (cosd(tilt) / cosd(ref_tilt)).powf(cosine_loss_exponent_tilt / 3.0);
                rotor_effective_velocity *= tilt_correction;

                // Step 5: Interpolate power from power curve (convert kW to W)
                let power = power_curve.interpolate(rotor_effective_velocity) * 1000.0;

                // Clip to rated power
                powers[[fi, ti]] = power.min(rated_power).max(0.0);
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
    /// Set wind conditions from WindData trait object
    ///
    /// This method accepts any object that implements the WindData trait,
    /// including TimeSeries, WindRose, WindRoseRwg and custom implementations.
    pub fn set_wind_data(&mut self, wind_data: &dyn WindData) -> crate::Result<()> {
        self.flow_field = FlowField::from_wind_data(wind_data);

        self.farm
            .initialize_control_arrays(self.flow_field.n_findex);
        self.grid = None;

        Ok(())
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
            let freq = frequencies.row(i).sum();
            total_energy += power * freq * hours_per_year;
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
                    ct_array[[fi, ti]] = turbine.turbine_type.get_ct(velocity);
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
            solver: SolverConfig::default(),
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
            solver: SolverConfig::default(),
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
            solver: SolverConfig::default(),
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
