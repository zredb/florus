//! Layout Optimization Scipy Module
//!
//! Provides scipy-style layout optimization using coordinate descent algorithm.
//! This module corresponds to layout_optimization_scipy.py in Python FLORIS v4.6.

use super::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer,
};
use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2, Float};
use crate::Result;

/// The yaw angles of the geometric yaw optimization for the current layout
/// (Python's `_get_geoyaw_angles`).
pub fn geometric_yaw_angles(fmodel: &mut FlorisModel) -> Result<Array2> {
    use crate::optimization::yaw_optimization::{YawOptimization, YawOptimizationGeometric};

    let mut optimizer = YawOptimizationGeometric::new();
    let result = optimizer.optimize(fmodel, None)?;
    Ok(result.yaw_angles)
}

/// LayoutOptimizationScipy - Layout optimization using scipy-style optimization
///
/// Corresponds to floris.optimization.LayoutOptimizationScipy in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct LayoutOptimizationScipy {
    /// FlorisModel for simulation
    fmodel: FlorisModel,
    /// Configuration
    config: LayoutOptimizationConfig,
    /// Boundary bounds
    xmin: Float,
    xmax: Float,
    ymin: Float,
    ymax: Float,
    /// Number of turbines
    n_turbines: usize,
    /// Initial AEP/AVP
    initial_value: Float,
}

impl LayoutOptimizationScipy {
    /// Create new scipy-style layout optimizer
    pub fn new(fmodel: &FlorisModel, boundaries: Boundary) -> Result<Self> {
        let config = LayoutOptimizationConfig::default();
        let n_turbines = fmodel.farm().n_turbines();
        let (xmin, xmax, ymin, ymax) = boundaries.bounds();

        // Use 0.0 as initial value if model hasn't been run yet
        // The actual value will be computed during optimization
        let initial_value = 0.0;

        Ok(Self {
            fmodel: fmodel.clone(),
            config,
            xmin,
            xmax,
            ymin,
            ymax,
            n_turbines,
            initial_value,
        })
    }

    /// Create with custom configuration
    pub fn with_config(mut self, config: LayoutOptimizationConfig) -> Self {
        self.config = config;
        self
    }

    /// Set minimum distance between turbines
    pub fn with_min_dist(mut self, min_dist: Float) -> Self {
        self.config.min_dist = Some(min_dist);
        self
    }

    /// Set solver
    pub fn with_solver(mut self, solver: String) -> Self {
        self.config.solver = solver;
        self
    }

    /// Get farm layout x coordinates
    fn farm_layout_x(&self) -> Array1 {
        let coords = self.fmodel.farm().coordinates();
        let n_turbines = coords.nrows();
        let mut x = Array1::zeros(n_turbines);

        for (i, coord) in coords.outer_iter().enumerate() {
            x[i] = coord[0];
        }

        x
    }

    /// Get farm layout y coordinates
    fn farm_layout_y(&self) -> Array1 {
        let coords = self.fmodel.farm().coordinates();
        let n_turbines = coords.nrows();
        let mut y = Array1::zeros(n_turbines);

        for (i, coord) in coords.outer_iter().enumerate() {
            y[i] = coord[1];
        }

        y
    }
}

impl LayoutOptimizer for LayoutOptimizationScipy {
    fn n_turbines(&self) -> usize {
        self.n_turbines
    }

    fn min_dist(&self) -> Float {
        // Python: "If no minimum distance is provided, assume a value of 2
        // rotor diameters" (layout_optimization_base.py:72)
        self.config.min_dist.unwrap_or_else(|| {
            2.0 * self
                .fmodel
                .farm()
                .rotor_diameters
                .iter()
                .copied()
                .fold(0.0, Float::max)
        })
    }

    fn boundaries(&self) -> &Boundary {
        &self.config.boundaries
    }

    fn calculate_objective(&self, x: &Array1, y: &Array1) -> Float {
        // Clone fmodel for this calculation to avoid borrowing issues
        let mut fmodel = self.fmodel.clone();

        // Set the new layout
        if let Err(e) = fmodel.set_layout(x, y) {
            eprintln!("Warning: Failed to set layout: {}", e);
            return 0.0;
        }

        // Run the model
        if let Err(e) = fmodel.run() {
            eprintln!("Warning: Failed to run model: {}", e);
            return 0.0;
        }

        // The yaw angles of the geometric yaw optimization (Python's
        // `_get_geoyaw_angles`), which couple the yaw angles to the layout
        if self.config.enable_geometric_yaw {
            match geometric_yaw_angles(&mut fmodel) {
                Ok(yaw_angles) => {
                    if let Err(e) = fmodel.set_yaw_angles(yaw_angles) {
                        eprintln!("Warning: Failed to set the yaw angles: {}", e);
                        return 0.0;
                    }
                    if let Err(e) = fmodel.run() {
                        eprintln!("Warning: Failed to run model: {}", e);
                        return 0.0;
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Failed to compute the geometric yaw angles: {}", e);
                    return 0.0;
                }
            }
        }

        // The objective is the annual energy (or value) production, weighted
        // with the frequencies of the wind data (Python's `get_farm_AEP` and
        // `get_farm_AVP`).
        if self.config.use_value {
            fmodel
                .get_expected_farm_value(None, None, None)
                .unwrap_or(0.0)
                * 8760.0
        } else {
            fmodel.get_expected_farm_power(None, None).unwrap_or(0.0) * 8760.0
        }
    }

    fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        let min_dist = self.min_dist();

        // Initial layout
        let initial_x = self.farm_layout_x();
        let initial_y = self.farm_layout_y();

        // Create optimization variables in normalized space [0, 1]
        let initial_value = self.calculate_objective(&initial_x, &initial_y);
        let mut x_norm: Vec<Float> = (0..self.n_turbines)
            .map(|i| self.normalize(initial_x[i], self.xmin, self.xmax))
            .collect();
        let mut y_norm: Vec<Float> = (0..self.n_turbines)
            .map(|i| self.normalize(initial_y[i], self.ymin, self.ymax))
            .collect();

        // A Hooke-Jeeves pattern search: the turbines are moved one at a time
        // along the eight compass directions (the first move that improves
        // the objective is taken), and the whole layout is then moved along
        // the direction of the last successful move (the "pattern move",
        // which allows the search to leave a local optimum).
        //
        // Python optimizes the same problem with scipy's SLSQP, which uses
        // finite-difference gradients and converges to the same optimum.
        let max_iter = self.config.max_iterations;
        let tolerance = self.config.tolerance;
        // The search is run twice: once with the pattern move and once
        // without it (the pattern move can leave a local optimum that the
        // plain pattern search does not reach); the best layout of the two
        // runs is kept.
        let mut total_iterations = 0;
        let mut best_x = x_norm.clone();
        let mut best_y = y_norm.clone();
        let mut best_value = Float::NEG_INFINITY;
        for use_pattern_move in [false, true] {
            let mut x_norm = x_norm.clone();
            let mut y_norm = y_norm.clone();
            let mut iteration = 0;
            let mut step = 0.1;
            let mut current_value = self.calculate_objective(
                &Array1::from_vec(x_norm.clone()),
                &Array1::from_vec(y_norm.clone()),
            );
            // The direction of the last successful move
            let mut direction = (0.0, 0.0);

            while iteration < max_iter && step > tolerance {
                iteration += 1;

                // The exploratory phase: the turbines are moved one at a time
                // along the eight compass directions
                let mut improved = false;
                'turbine: for ti in 0..self.n_turbines {
                    let original_x = x_norm[ti];
                    let original_y = y_norm[ti];

                    for k in 0..8 {
                        let angle = k as Float * std::f64::consts::FRAC_PI_4;
                        let (dx, dy) = (angle.cos() * step, angle.sin() * step);

                        x_norm[ti] = (original_x + dx).clamp(0.0, 1.0);
                        y_norm[ti] = (original_y + dy).clamp(0.0, 1.0);

                        let new_x: Array1 = x_norm
                            .iter()
                            .map(|&v| self.denormalize(v, self.xmin, self.xmax))
                            .collect();
                        let new_y: Array1 = y_norm
                            .iter()
                            .map(|&v| self.denormalize(v, self.ymin, self.ymax))
                            .collect();

                        // The turbines have to stay inside the boundary (the
                        // points on the boundary are allowed) and keep their
                        // minimum distance from each other
                        if self.check_boundary_distance(&new_x, &new_y)
                            && self.calculate_min_distance(&new_x, &new_y) >= min_dist
                        {
                            let new_value = self.calculate_objective(&new_x, &new_y);

                            if new_value > current_value {
                                current_value = new_value;
                                direction = (dx, dy);
                                improved = true;
                                break 'turbine;
                            }
                        }

                        // Revert
                        x_norm[ti] = original_x;
                        y_norm[ti] = original_y;
                    }
                }

                if improved {
                    continue;
                }

                if !use_pattern_move {
                    // No improvement: reduce the step size
                    step *= 0.5;
                    continue;
                }

                // The pattern move: the whole layout is moved along the
                // direction of the last successful move
                if direction.0 != 0.0 || direction.1 != 0.0 {
                    let original_x = x_norm.clone();
                    let original_y = y_norm.clone();
                    for ti in 0..self.n_turbines {
                        x_norm[ti] = (x_norm[ti] + direction.0 * step).clamp(0.0, 1.0);
                        y_norm[ti] = (y_norm[ti] + direction.1 * step).clamp(0.0, 1.0);
                    }

                    let new_x: Array1 = x_norm
                        .iter()
                        .map(|&v| self.denormalize(v, self.xmin, self.xmax))
                        .collect();
                    let new_y: Array1 = y_norm
                        .iter()
                        .map(|&v| self.denormalize(v, self.ymin, self.ymax))
                        .collect();

                    if self.check_boundary_distance(&new_x, &new_y)
                        && self.calculate_min_distance(&new_x, &new_y) >= min_dist
                    {
                        let new_value = self.calculate_objective(&new_x, &new_y);
                        if new_value > current_value {
                            current_value = new_value;
                            continue;
                        }
                    }

                    // Revert the pattern move
                    x_norm = original_x;
                    y_norm = original_y;
                }

                // No improvement: reduce the step size
                step *= 0.5;
                direction = (0.0, 0.0);
            }
            total_iterations += iteration;
            // Keep the best layout of the two searches
            if current_value > best_value {
                best_value = current_value;
                best_x = x_norm;
                best_y = y_norm;
            }
        }

        // The value of the initial layout (used for the improvement)
        self.initial_value = initial_value;
        x_norm = best_x;
        y_norm = best_y;

        // Convert back to physical coordinates
        let opt_x: Array1 = x_norm
            .iter()
            .map(|&v| self.denormalize(v, self.xmin, self.xmax))
            .collect();
        let opt_y: Array1 = y_norm
            .iter()
            .map(|&v| self.denormalize(v, self.ymin, self.ymax))
            .collect();

        let final_value = self.calculate_objective(&opt_x, &opt_y);
        let improvement_pct = if self.initial_value > 0.0 {
            100.0 * (final_value - self.initial_value) / self.initial_value
        } else {
            0.0
        };
        Ok(LayoutOptimizationResult {
            x: opt_x,
            y: opt_y,
            value: final_value,
            iterations: total_iterations,
            improvement_pct,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Farm;
    use crate::core::FlowField;
    use crate::floris_config::SolverConfig;
    use crate::floris_model::FlorisModel;
    use crate::types::Array1;

    #[test]
    fn test_scipy_optimizer_creation() {
        use crate::floris_config::{
            FarmConfig, FlorisConfig, FlowFieldConfig, SolverConfig, WakeConfig,
        };

        let farm_config = FarmConfig {
            layout_x: vec![0.0, 500.0],
            layout_y: vec![0.0, 0.0],
            turbine_type: vec!["nrel_5MW".to_string(); 2],
        };

        let flow_field_config = FlowFieldConfig {
            wind_speeds: vec![8.0],
            wind_directions: vec![270.0],
            turbulence_intensities: vec![0.06],
            wind_shear: 0.14,
            wind_veer: 0.0,
            air_density: 1.225,
            reference_wind_height: 90.0,
            multidim_conditions: None,
            heterogeneous_inflow_config: None,
        };

        let solver_config = SolverConfig::default();
        let wake_config = WakeConfig::default();

        let config = FlorisConfig {
            name: "test".to_string(),
            description: Some("test".to_string()),
            floris_version: "v4".to_string(),
            logging: Default::default(),
            solver: solver_config,
            farm: farm_config,
            flow_field: flow_field_config,
            wake: wake_config,
            turbine_library: "turbine_library".to_string(),
        };

        let model = FlorisModel::from_config(config).unwrap();

        let boundary = Boundary::Rectangle {
            min_x: 0.0,
            max_x: 1000.0,
            min_y: 0.0,
            max_y: 1000.0,
        };

        let optimizer = LayoutOptimizationScipy::new(&model, boundary).unwrap();
        assert_eq!(optimizer.n_turbines(), 2);
    }
}
