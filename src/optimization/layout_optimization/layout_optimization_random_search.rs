//! Layout Optimization Random Search Module
//!
//! Provides random search layout optimization using genetic algorithm-like approach.
//! This module corresponds to layout_optimization_random_search.py in Python FLORIS v4.6.

use super::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer,
};
use crate::floris_model::FlorisModel;
use crate::types::{Array1, Float};
use crate::Result;

/// LayoutOptimizationRandomSearch - Random search layout optimization
///
/// Corresponds to floris.optimization.LayoutOptimizationRandomSearch in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct LayoutOptimizationRandomSearch {
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
    /// Number of individuals
    n_individuals: usize,
    /// Initial AEP/AVP
    initial_value: Float,
    /// Random seed
    random_seed: Option<u64>,
}

impl LayoutOptimizationRandomSearch {
    /// Create new random search layout optimizer
    pub fn new(fmodel: &FlorisModel, boundaries: Boundary, n_individuals: usize) -> Result<Self> {
        let config = LayoutOptimizationConfig::default();
        let n_turbines = fmodel.farm().n_turbines();
        let (xmin, xmax, ymin, ymax) = boundaries.bounds();

        // Use 0.0 as initial value if model hasn't been run yet
        let initial_value = 0.0;

        Ok(Self {
            fmodel: fmodel.clone(),
            config,
            xmin,
            xmax,
            ymin,
            ymax,
            n_turbines,
            n_individuals,
            initial_value,
            random_seed: None,
        })
    }

    /// Set configuration
    pub fn with_config(mut self, config: LayoutOptimizationConfig) -> Self {
        self.config = config;
        self
    }

    /// Set minimum distance between turbines
    pub fn with_min_dist(mut self, min_dist: Float) -> Self {
        self.config.min_dist = Some(min_dist);
        self
    }

    /// Set random seed
    pub fn with_random_seed(mut self, seed: u64) -> Self {
        self.random_seed = Some(seed);
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

    /// Generate distance-based layout
    fn generate_distance_based_layout(&self, min_dist: Float) -> Option<(Vec<Float>, Vec<Float>)> {
        // Greedy algorithm: place each turbine as far as possible from existing ones
        let mut layout_x: Vec<Float> = Vec::new();
        let mut layout_y: Vec<Float> = Vec::new();
        let grid_step = min_dist * 0.5;

        for _ in 0..self.n_turbines {
            let mut best_x = 0.0;
            let mut best_y = 0.0;
            let mut max_min_dist = 0.0;

            // Search over a grid
            for x in (self.xmin as usize..self.xmax as usize)
                .step_by(grid_step as usize)
                .map(|v| v as Float)
            {
                for y in (self.ymin as usize..self.ymax as usize)
                    .step_by(grid_step as usize)
                    .map(|v| v as Float)
                {
                    // Check if in boundary
                    if !self.boundaries().contains(x, y) {
                        continue;
                    }

                    // Calculate minimum distance to existing turbines
                    let mut min_dist_to_existing = Float::INFINITY;
                    for (ex, ey) in layout_x.iter().zip(layout_y.iter()) {
                        let dx = x - ex;
                        let dy = y - ey;
                        let dist = (dx * dx + dy * dy).sqrt();
                        if dist < min_dist_to_existing {
                            min_dist_to_existing = dist;
                        }
                    }

                    if min_dist_to_existing > max_min_dist {
                        max_min_dist = min_dist_to_existing;
                        best_x = x;
                        best_y = y;
                    }
                }
            }

            if max_min_dist > 0.0 {
                layout_x.push(best_x);
                layout_y.push(best_y);
            } else {
                // Could not find valid position, use random
                layout_x.push(self.xmin + (self.xmax - self.xmin) * 0.5);
                layout_y.push(self.ymin + (self.ymax - self.ymin) * 0.5);
            }
        }

        Some((layout_x, layout_y))
    }
}

impl LayoutOptimizer for LayoutOptimizationRandomSearch {
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

        if let Err(e) = fmodel.set_layout(x, y) {
            eprintln!("Warning: Failed to set layout: {}", e);
            return 0.0;
        }

        if let Err(e) = fmodel.run() {
            eprintln!("Warning: Failed to run model: {}", e);
            return 0.0;
        }

        if self.config.use_value {
            fmodel.get_farm_avp()
        } else {
            fmodel.get_farm_aep_uniform(8760.0)
        }
    }

    fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        let min_dist = self.min_dist();

        // Get initial layout
        let initial_x = self.farm_layout_x();
        let initial_y = self.farm_layout_y();

        // Distance-based step probability mass function
        let d_max = ((self.xmax - self.xmin).min(self.ymax - self.ymin) / 2.0).min(min_dist * 3.0);
        let step_distances: Vec<Float> = (0..50)
            .map(|i| (i as Float + 1.0) * min_dist * 0.1)
            .filter(|&d| d < d_max)
            .chain(std::iter::once(d_max))
            .collect();

        let jump_prob = 0.05;
        let base_prob = (1.0 - jump_prob) / (step_distances.len() as Float - 1.0);
        let _step_probs: Vec<Float> = step_distances[..step_distances.len() - 1]
            .iter()
            .map(|_| base_prob)
            .chain(std::iter::once(jump_prob))
            .collect();

        // Initialize candidates with distance-based initialization
        let mut candidates_x: Vec<Vec<Float>> = Vec::new();
        let mut candidates_y: Vec<Vec<Float>> = Vec::new();
        let mut candidate_values: Vec<Float> = Vec::new();

        // First candidate is the initial layout
        candidates_x.push(initial_x.to_vec());
        candidates_y.push(initial_y.to_vec());
        candidate_values.push(self.initial_value);

        // Generate remaining candidates with distance-based initialization
        for _i in 1..self.n_individuals {
            if let Some((x, y)) = self.generate_distance_based_layout(min_dist) {
                candidates_x.push(x.clone());
                candidates_y.push(y.clone());
                let val = self.calculate_objective(&Array1::from_vec(x), &Array1::from_vec(y));
                candidate_values.push(val);
            } else {
                // Fallback to initial layout
                candidates_x.push(initial_x.to_vec());
                candidates_y.push(initial_y.to_vec());
                candidate_values.push(self.initial_value);
            }
        }

        // Random search optimization
        let mut iteration = 0;
        let max_iter = self.config.max_iterations;

        for _ in 0..max_iter {
            iteration += 1;
            let mut improved = false;

            for i in 0..self.n_individuals {
                // Select a random turbine using deterministic pseudo-random for reproducibility
                let turbine_idx = if let Some(seed) = self.random_seed {
                    ((iteration * 12345 + i * 67890 + seed as usize) % self.n_turbines) as usize
                } else {
                    (iteration * 12345 + i * 67890) % self.n_turbines
                };

                // Select random direction
                let angle = if let Some(seed) = self.random_seed {
                    ((iteration * 54321 + i * 98765 + seed as usize) % 1000) as Float / 1000.0
                        * 2.0
                        * std::f64::consts::PI
                } else {
                    ((iteration * 54321 + i * 98765) % 1000) as Float / 1000.0
                        * 2.0
                        * std::f64::consts::PI
                };

                // Select random distance
                let dist_idx = if let Some(seed) = self.random_seed {
                    ((iteration * 11111 + i * 22222 + seed as usize) % step_distances.len())
                        as usize
                } else {
                    ((iteration * 11111 + i * 22222) % step_distances.len()) as usize
                };
                let dist = step_distances[dist_idx];

                // Create new candidate layout
                let mut new_x = candidates_x[i].clone();
                let mut new_y = candidates_y[i].clone();

                new_x[turbine_idx] += angle.cos() * dist;
                new_y[turbine_idx] += angle.sin() * dist;

                // Check constraints
                let new_x_arr = Array1::from_vec(new_x.clone());
                let new_y_arr = Array1::from_vec(new_y.clone());

                if self.check_boundary_constraints(&new_x_arr, &new_y_arr)
                    && self.calculate_min_distance(&new_x_arr, &new_y_arr) >= min_dist
                {
                    let new_value = self.calculate_objective(&new_x_arr, &new_y_arr);

                    if new_value > candidate_values[i] {
                        candidates_x[i] = new_x;
                        candidates_y[i] = new_y;
                        candidate_values[i] = new_value;
                        improved = true;
                    }
                }
            }

            if !improved {
                break;
            }
        }

        // Find the best candidate
        let best_idx = candidate_values
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(idx, _)| idx)
            .unwrap_or(0);

        let opt_x = Array1::from_vec(candidates_x[best_idx].clone());
        let opt_y = Array1::from_vec(candidates_y[best_idx].clone());
        let final_value = candidate_values[best_idx];
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
    use crate::floris_config::{FarmConfig, FlorisConfig, FlowFieldConfig, SolverConfig, WakeConfig};
    use crate::floris_model::FlorisModel;
    use crate::types::Array1;

    #[test]
    fn test_random_search_optimizer_creation() {
        let farm_config = FarmConfig {
            layout_x: vec![0.0, 500.0, 1000.0],
            layout_y: vec![0.0, 0.0, 0.0],
            turbine_type: vec!["nrel_5MW".to_string(); 3],
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

        let optimizer = LayoutOptimizationRandomSearch::new(&model, boundary, 10).unwrap();
        assert_eq!(optimizer.n_turbines(), 3);
        assert_eq!(optimizer.n_individuals, 10);
    }
}
