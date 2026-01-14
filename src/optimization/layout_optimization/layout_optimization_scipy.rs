//! Layout Optimization Scipy Module
//!
//! Provides scipy-style layout optimization using coordinate descent algorithm.
//! This module corresponds to layout_optimization_scipy.py in Python FLORIS v4.6.

use crate::floris_model::FlorisModel;
use super::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer,
};
use crate::types::{Array1, Float};
use crate::Result;

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
        let n_turbines = fmodel.farm.n_turbines();
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
        let coords = self.fmodel.farm.coordinates();
        let n_turbines = coords.nrows();
        let mut x = Array1::zeros(n_turbines);

        for (i, coord) in coords.outer_iter().enumerate() {
            x[i] = coord[0];
        }

        x
    }

    /// Get farm layout y coordinates
    fn farm_layout_y(&self) -> Array1 {
        let coords = self.fmodel.farm.coordinates();
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
        self.config.min_dist.unwrap_or(500.0)
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

        // Return AEP or AVP using uniform frequencies
        if self.config.use_value {
            fmodel.get_farm_avp()
        } else {
            fmodel.get_farm_aep_uniform(8760.0)
        }
    }

    fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        let min_dist = self.min_dist();
        
        // Initial layout
        let initial_x = self.farm_layout_x();
        let initial_y = self.farm_layout_y();
        
        // Create optimization variables in normalized space [0, 1]
        let mut x_norm: Vec<Float> = (0..self.n_turbines)
            .map(|i| self.normalize(initial_x[i], self.xmin, self.xmax))
            .collect();
        let mut y_norm: Vec<Float> = (0..self.n_turbines)
            .map(|i| self.normalize(initial_y[i], self.ymin, self.ymax))
            .collect();
        
        // Simple coordinate descent optimization
        let max_iter = self.config.max_iterations;
        let tolerance = self.config.tolerance;
        let mut iteration = 0;
        
        for _ in 0..max_iter {
            iteration += 1;
            let mut improved = false;
            let step_size = 0.05 / (1.0 + iteration as Float / 20.0);
            
            // Try moving each turbine in each direction
            for ti in 0..self.n_turbines {
                let original_x = x_norm[ti];
                let original_y = y_norm[ti];
                
                // Try small steps in x and y directions
                let directions = [
                    (step_size, 0.0),
                    (-step_size, 0.0),
                    (0.0, step_size),
                    (0.0, -step_size),
                ];
                
                let current_value = self.calculate_objective(
                    &Array1::from_vec(x_norm.clone()),
                    &Array1::from_vec(y_norm.clone()),
                );
                
                for &(dx, dy) in &directions {
                    x_norm[ti] = (original_x + dx).clamp(0.0, 1.0);
                    y_norm[ti] = (original_y + dy).clamp(0.0, 1.0);
                    
                    let new_x: Array1 = x_norm.iter().map(|&v| self.denormalize(v, self.xmin, self.xmax)).collect();
                    let new_y: Array1 = y_norm.iter().map(|&v| self.denormalize(v, self.ymin, self.ymax)).collect();
                    
                    // Check minimum distance constraint
                    let min_dist_actual = self.calculate_min_distance(&new_x, &new_y);
                    let in_bounds = self.check_boundary_constraints(&new_x, &new_y);
                    
                    if min_dist_actual >= min_dist && in_bounds {
                        let new_value = self.calculate_objective(&new_x, &new_y);
                        
                        // For maximization, check if new value is better
                        if new_value > current_value {
                            improved = true;
                            break;
                        }
                    }
                    
                    // Revert
                    x_norm[ti] = original_x;
                    y_norm[ti] = original_y;
                }
            }
            
            if !improved {
                break;
            }
        }
        
        // Convert back to physical coordinates
        let opt_x: Array1 = x_norm.iter().map(|&v| self.denormalize(v, self.xmin, self.xmax)).collect();
        let opt_y: Array1 = y_norm.iter().map(|&v| self.denormalize(v, self.ymin, self.ymax)).collect();
        
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
            iterations: iteration,
            improvement_pct,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;
    use crate::core::Farm;
    use crate::floris_model::FlorisModel;
    use crate::core::FlowField;

    #[test]
    fn test_scipy_optimizer_creation() {
        // Create a simple farm for testing
        let layout_x = Array1::from_vec(vec![0.0, 500.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 2];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();
        
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
        ).unwrap();
        
        let model = FlorisModel {
            farm,
            flow_field,
            state: crate::core::State::new(),
            grid: None,
            solver_type: "turbine_grid".to_string(),
            model_manager: None,
        };
        
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
