//! Layout Optimization Boundary Grid Module
//!
//! Provides grid-based layout optimization within specified boundaries.
//! This module corresponds to layout_optimization_boundary_grid.py in Python FLORIS v4.6.

use crate::floris_model::FlorisModel;
use super::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer,
};
use crate::types::{Array1, Float};
use crate::Result;
use rand::prelude::SliceRandom;
use rand::Rng;
use serde::{Deserialize, Serialize};

/// Grid-based layout optimization configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridOptimizationConfig {
    /// Base configuration
    pub base: LayoutOptimizationConfig,
    /// Grid resolution (number of points per side)
    pub grid_resolution: usize,
    /// Enable smart start (grid-based initial positions)
    pub smart_start: bool,
    /// Smart start distance between turbines
    pub smart_start_distance: Float,
}

impl Default for GridOptimizationConfig {
    fn default() -> Self {
        Self {
            base: LayoutOptimizationConfig::default(),
            grid_resolution: 20,
            smart_start: true,
            smart_start_distance: 500.0,
        }
    }
}

/// LayoutOptimizationBoundaryGrid - Grid-based layout optimization
///
/// Optimizes turbine positions on a predefined grid within the boundary.
/// Corresponds to floris.optimization.LayoutOptimizationBoundaryGrid in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct LayoutOptimizationBoundaryGrid {
    /// FlorisModel for simulation
    fmodel: FlorisModel,
    /// Configuration
    config: GridOptimizationConfig,
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

impl LayoutOptimizationBoundaryGrid {
    /// Create new boundary grid layout optimizer
    pub fn new(fmodel: &FlorisModel, boundaries: Boundary, grid_resolution: usize) -> Result<Self> {
        let (xmin, xmax, ymin, ymax) = boundaries.bounds();
        let n_turbines = fmodel.farm.n_turbines();

        Ok(Self {
            fmodel: fmodel.clone(),
            config: GridOptimizationConfig {
                base: LayoutOptimizationConfig::default(),
                grid_resolution,
                smart_start: true,
                smart_start_distance: 500.0,
            },
            xmin,
            xmax,
            ymin,
            ymax,
            n_turbines,
            initial_value: 0.0,
        })
    }

    /// Create with custom grid configuration
    pub fn with_config(mut self, config: GridOptimizationConfig) -> Self {
        self.config = config;
        self
    }

    /// Set minimum distance between turbines
    pub fn with_min_dist(mut self, min_dist: Float) -> Self {
        self.config.base.min_dist = Some(min_dist);
        self
    }

    /// Set smart start parameters
    pub fn with_smart_start(mut self, enabled: bool, distance: Float) -> Self {
        self.config.smart_start = enabled;
        self.config.smart_start_distance = distance;
        self
    }

    /// Generate grid points within the boundary
    fn generate_grid_points(&self) -> Vec<(Float, Float)> {
        let resolution = self.config.grid_resolution;
        let dx = (self.xmax - self.xmin) / (resolution as Float + 1.0);
        let dy = (self.ymax - self.ymin) / (resolution as Float + 1.0);

        let mut points = Vec::new();
        for i in 1..=resolution {
            for j in 1..=resolution {
                let x = self.xmin + i as Float * dx;
                let y = self.ymin + j as Float * dy;
                if self.config.base.boundaries.contains(x, y) {
                    points.push((x, y));
                }
            }
        }
        points
    }

    /// Generate smart start positions using grid-based initialization
    fn generate_smart_start(&self) -> (Array1, Array1) {
        let grid_points = self.generate_grid_points();
        let _min_dist = self.config.smart_start_distance;

        if grid_points.is_empty() {
            // Fallback to random initialization
            let mut rng = rand::thread_rng();
            let mut x = Array1::zeros(self.n_turbines);
            let mut y = Array1::zeros(self.n_turbines);

            for i in 0..self.n_turbines {
                x[i] = self.xmin + (self.xmax - self.xmin) * rng.gen::<Float>();
                y[i] = self.ymin + (self.ymax - self.ymin) * rng.gen::<Float>();
            }

            return (x, y);
        }

        // Greedy algorithm to select positions
        let mut selected_indices = Vec::new();
        let mut selected_points = Vec::new();

        // Select first point at center or random
        let center_idx = grid_points.len() / 2;
        selected_indices.push(center_idx);
        selected_points.push(grid_points[center_idx]);

        // Select remaining points based on maximum minimum distance
        while selected_points.len() < self.n_turbines && !grid_points.is_empty() {
            let mut best_idx = 0;
            let mut best_min_dist = 0.0;

            for (idx, point) in grid_points.iter().enumerate() {
                if selected_indices.contains(&idx) {
                    continue;
                }

                // Calculate minimum distance to already selected points
                let min_dist_to_selected = selected_points.iter()
                    .map(|selected| {
                        let dx = point.0 - selected.0;
                        let dy = point.1 - selected.1;
                        (dx * dx + dy * dy).sqrt()
                    })
                    .fold(Float::INFINITY, |a, b| a.min(b));

                if min_dist_to_selected > best_min_dist {
                    best_min_dist = min_dist_to_selected;
                    best_idx = idx;
                }
            }

            selected_indices.push(best_idx);
            selected_points.push(grid_points[best_idx]);
        }

        // Convert to arrays
        let mut x = Array1::zeros(self.n_turbines);
        let mut y = Array1::zeros(self.n_turbines);

        for (i, point) in selected_points.iter().enumerate().take(self.n_turbines) {
            x[i] = point.0;
            y[i] = point.1;
        }

        (x, y)
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

    /// Find the best positions on the grid using exhaustive search
    fn find_optimal_grid_positions(&self, initial_x: &Array1, initial_y: &Array1) -> (Array1, Array1) {
        let grid_points = self.generate_grid_points();
        let min_dist = self.config.base.min_dist.unwrap_or(500.0);

        // Use smart start or initial layout
        let (mut best_x, mut best_y) = if self.config.smart_start {
            self.generate_smart_start()
        } else {
            (initial_x.clone(), initial_y.clone())
        };

        // If grid is coarse, try grid-based search
        if grid_points.len() >= self.n_turbines && grid_points.len() < 1000 {
            // Exhaustive search on grid points
            let mut best_value = self.calculate_objective(&best_x, &best_y);

            // Generate all combinations of grid points
            let n_combinations = if grid_points.len() > 20 {
                10000 // Limit for large grids
            } else {
                // Use binomial coefficient approximation
                let _count = 0usize;
                let mut c = vec![0.0; self.n_turbines + 1];
                c[0] = 1.0;
                for i in 1..=self.n_turbines {
                    c[i] = c[i - 1] * (grid_points.len() - i + 1) as Float / i as Float;
                }
                c[self.n_turbines].min(10000.0) as usize
            };

            // Random sampling of grid combinations
            let mut rng = rand::thread_rng();
            for _ in 0..n_combinations.min(1000) {
                // Randomly select n_turbines points
                let mut indices: Vec<usize> = (0..grid_points.len()).collect();
                indices.shuffle(&mut rng);
                let selected_indices = &indices[..self.n_turbines];

                let mut test_x = Array1::zeros(self.n_turbines);
                let mut test_y = Array1::zeros(self.n_turbines);

                // Check minimum distance
                let mut valid = true;
                for (i, &idx) in selected_indices.iter().enumerate() {
                    test_x[i] = grid_points[idx].0;
                    test_y[i] = grid_points[idx].1;
                }

                // Verify minimum distance constraint
                for i in 0..self.n_turbines {
                    for j in (i + 1)..self.n_turbines {
                        let dx = test_x[i] - test_x[j];
                        let dy = test_y[i] - test_y[j];
                        if (dx * dx + dy * dy).sqrt() < min_dist {
                            valid = false;
                            break;
                        }
                    }
                    if !valid {
                        break;
                    }
                }

                if valid {
                    let value = self.calculate_objective(&test_x, &test_y);
                    if value > best_value {
                        best_value = value;
                        best_x = test_x.clone();
                        best_y = test_y.clone();
                    }
                }
            }
        }

        (best_x, best_y)
    }
}

impl LayoutOptimizer for LayoutOptimizationBoundaryGrid {
    fn n_turbines(&self) -> usize {
        self.n_turbines
    }

    fn min_dist(&self) -> Float {
        self.config.base.min_dist.unwrap_or(500.0)
    }

    fn boundaries(&self) -> &Boundary {
        &self.config.base.boundaries
    }

    fn calculate_objective(&self, x: &Array1, y: &Array1) -> Float {
        let mut fmodel = self.fmodel.clone();

        if let Err(_e) = fmodel.set_layout(x, y) {
            return 0.0;
        }

        if let Err(_e) = fmodel.run() {
            return 0.0;
        }

        if self.config.base.use_value {
            fmodel.get_farm_avp()
        } else {
            fmodel.get_farm_aep_uniform(8760.0)
        }
    }

    fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        // Get initial layout
        let initial_x = self.farm_layout_x();
        let initial_y = self.farm_layout_y();

        // Find optimal positions using grid-based search
        let (opt_x, opt_y) = self.find_optimal_grid_positions(&initial_x, &initial_y);

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
            iterations: 1,
            improvement_pct,
        })
    }
}

/// Mixed integer optimizer combining grid and continuous optimization
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LayoutOptimizationMixedInteger {
    /// FlorisModel for simulation
    fmodel: FlorisModel,
    /// Configuration
    config: GridOptimizationConfig,
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

impl LayoutOptimizationMixedInteger {
    /// Create new mixed integer layout optimizer
    pub fn new(fmodel: &FlorisModel, boundaries: Boundary, grid_resolution: usize) -> Result<Self> {
        let (xmin, xmax, ymin, ymax) = boundaries.bounds();
        let n_turbines = fmodel.farm.n_turbines();

        Ok(Self {
            fmodel: fmodel.clone(),
            config: GridOptimizationConfig {
                base: LayoutOptimizationConfig::default(),
                grid_resolution,
                smart_start: true,
                smart_start_distance: 500.0,
            },
            xmin,
            xmax,
            ymin,
            ymax,
            n_turbines,
            initial_value: 0.0,
        })
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

impl LayoutOptimizer for LayoutOptimizationMixedInteger {
    fn n_turbines(&self) -> usize {
        self.n_turbines
    }

    fn min_dist(&self) -> Float {
        self.config.base.min_dist.unwrap_or(500.0)
    }

    fn boundaries(&self) -> &Boundary {
        &self.config.base.boundaries
    }

    fn calculate_objective(&self, x: &Array1, y: &Array1) -> Float {
        let mut fmodel = self.fmodel.clone();

        if let Err(_e) = fmodel.set_layout(x, y) {
            return 0.0;
        }

        if let Err(_e) = fmodel.run() {
            return 0.0;
        }

        if self.config.base.use_value {
            fmodel.get_farm_avp()
        } else {
            fmodel.get_farm_aep_uniform(8760.0)
        }
    }

    fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        let _initial_x = self.farm_layout_x();
        let _initial_y = self.farm_layout_y();
        let min_dist = self.min_dist();

        // First pass: grid-based selection for discrete positions
        let mut grid_optimizer = LayoutOptimizationBoundaryGrid::new(&self.fmodel, self.boundaries().clone(), self.config.grid_resolution)
            .unwrap()
            .with_min_dist(min_dist)
            .with_smart_start(self.config.smart_start, self.config.smart_start_distance);
        
        let grid_result = grid_optimizer.optimize()?;
        
        // Second pass: continuous refinement around grid positions
        let mut best_x = grid_result.x.clone();
        let mut best_y = grid_result.y.clone();
        let mut best_value = self.calculate_objective(&best_x, &best_y);

        // Small neighborhood search
        let refinement_range = (self.xmax - self.xmin) / self.config.grid_resolution as Float;
        
        // Try small perturbations around best grid positions
        let n_refinements = 50;
        let mut rng = rand::thread_rng();
        
        for _ in 0..n_refinements {
            let mut test_x = best_x.clone();
            let mut test_y = best_y.clone();
            
            // Random turbine to move
            let turbine_idx = rng.gen_range(0..self.n_turbines);
            
            // Small random displacement
            test_x[turbine_idx] += (rng.gen::<Float>() - 0.5) * 2.0 * refinement_range;
            test_y[turbine_idx] += (rng.gen::<Float>() - 0.5) * 2.0 * refinement_range;
            
            // Clamp to bounds
            test_x[turbine_idx] = test_x[turbine_idx].clamp(self.xmin, self.xmax);
            test_y[turbine_idx] = test_y[turbine_idx].clamp(self.ymin, self.ymax);
            
            // Check constraints
            if self.check_boundary_constraints(&test_x, &test_y) 
                && self.calculate_min_distance(&test_x, &test_y) >= min_dist {
                
                let value = self.calculate_objective(&test_x, &test_y);
                if value > best_value {
                    best_value = value;
                    best_x = test_x;
                    best_y = test_y;
                }
            }
        }

        Ok(LayoutOptimizationResult {
            x: best_x,
            y: best_y,
            value: best_value,
            iterations: 2,
            improvement_pct: 0.0,
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
    fn test_boundary_grid_optimizer_creation() {
        let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 3];
        
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
            max_x: 2000.0,
            min_y: 0.0,
            max_y: 2000.0,
        };
        
        let optimizer = LayoutOptimizationBoundaryGrid::new(&model, boundary, 10).unwrap();
        assert_eq!(optimizer.n_turbines(), 3);
    }

    #[test]
    fn test_mixed_integer_optimizer_creation() {
        let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 3];
        
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
            max_x: 2000.0,
            min_y: 0.0,
            max_y: 2000.0,
        };
        
        let optimizer = LayoutOptimizationMixedInteger::new(&model, boundary, 10).unwrap();
        assert_eq!(optimizer.n_turbines(), 3);
    }
}
