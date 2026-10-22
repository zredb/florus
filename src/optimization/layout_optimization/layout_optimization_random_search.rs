//! Layout Optimization Random Search Module
//!
//! Provides random search layout optimization using genetic algorithm-like approach.
//! This module corresponds to layout_optimization_random_search.py in Python FLORIS v4.6.

use super::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer,
};
use crate::floris_model::FlorisModel;
use crate::numpy_random::NumpyRandomState;
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

    /// The probability mass function of the step distances (Python's
    /// `_process_dist_pmf`): a uniform distribution between 0 and 2 D, with
    /// some extra mass for a large jump.
    fn distance_pmf(&self) -> (Vec<Float>, Vec<Float>) {
        let rotor_diameter = self
            .fmodel
            .farm()
            .rotor_diameters
            .iter()
            .copied()
            .fold(0.0, Float::max);
        let jump_dist = (self.xmax - self.xmin).min(self.ymax - self.ymin) / 2.0;
        let jump_prob = 0.05;

        // np.append(np.linspace(0.0, 2.0 * D, 99), jump_dist)
        let mut distances: Vec<Float> = (0..99)
            .map(|i| 2.0 * rotor_diameter * i as Float / 98.0)
            .collect();
        distances.push(jump_dist);

        // np.append((1 - jump_prob) / len(d) * ones(len(d) - 1), jump_prob)
        let n = distances.len();
        let mut probabilities: Vec<Float> =
            (0..n - 1).map(|_| (1.0 - jump_prob) / n as Float).collect();
        probabilities.push(jump_prob);
        let total: Float = probabilities.iter().sum();
        let probabilities: Vec<Float> = probabilities.iter().map(|p| p / total).collect();

        (distances, probabilities)
    }

    /// The initial layout of an individual: the first turbine is placed at
    /// random and the remaining ones as far as possible from the others
    /// (Python's `_gen_dist_based_init`).
    fn distance_based_init(
        &self,
        rng: &mut NumpyRandomState,
        grid_step: Float,
    ) -> (Vec<Float>, Vec<Float>) {
        // The first turbine is placed randomly (and placed again if it is
        // outside of the boundary)
        let mut init_x =
            rng.randint(self.xmin as i64, self.xmax.max(self.xmin + 1.0) as i64) as Float;
        let mut init_y =
            rng.randint(self.ymin as i64, self.ymax.max(self.ymin + 1.0) as i64) as Float;
        let mut attempts = 0;
        while !self.boundaries().contains(init_x, init_y) && attempts < 1000 {
            init_x = rng.randint(self.xmin as i64, self.xmax.max(self.xmin + 1.0) as i64) as Float;
            init_y = rng.randint(self.ymin as i64, self.ymax.max(self.ymin + 1.0) as i64) as Float;
            attempts += 1;
        }

        let mut layout_x = vec![init_x];
        let mut layout_y = vec![init_y];

        // The remaining turbines are placed as far as possible from the others
        for _ in 1..self.n_turbines {
            let mut max_dist = 0.0;
            let mut save_x = 0.0;
            let mut save_y = 0.0;
            let mut x = self.xmin;
            while x < self.xmax {
                let mut y = self.ymin;
                while y < self.ymax {
                    if self.boundaries().contains(x, y) {
                        // The distance to the closest of the existing turbines
                        let mut closest = Float::INFINITY;
                        for (ex, ey) in layout_x.iter().zip(layout_y.iter()) {
                            let dx = x - ex;
                            let dy = y - ey;
                            let dist = (dx * dx + dy * dy).sqrt();
                            if dist < closest {
                                closest = dist;
                            }
                        }
                        if closest > max_dist {
                            max_dist = closest;
                            save_x = x;
                            save_y = y;
                        }
                    }
                    y += grid_step;
                }
                x += grid_step;
            }
            layout_x.push(save_x);
            layout_y.push(save_y);
        }

        (layout_x, layout_y)
    }

    /// The optimization of a single individual: the turbines are nudged in
    /// random directions until the step budget is exhausted (Python's
    /// `_single_individual_opt`).
    fn optimize_individual(
        &self,
        rng: &mut NumpyRandomState,
        layout_x: &[Float],
        layout_y: &[Float],
        pmf: &(Vec<Float>, Vec<Float>),
        n_iterations: usize,
        min_dist: Float,
    ) -> (Vec<Float>, Vec<Float>, Float) {
        let mut x = layout_x.to_vec();
        let mut y = layout_y.to_vec();
        let mut current_value =
            self.calculate_objective(&Array1::from_vec(x.clone()), &Array1::from_vec(y.clone()));

        for _ in 0..n_iterations {
            // Randomly select a turbine to nudge, a direction and a distance
            let turbine = rng.randint(0, self.n_turbines as i64) as usize;
            let direction = rng.uniform(0.0, 2.0 * std::f64::consts::PI);
            let distance = rng.choice(&pmf.0, &pmf.1);

            let test_x = x[turbine] + direction.cos() * distance;
            let test_y = y[turbine] + direction.sin() * distance;

            // The new turbine has to be inside the boundary
            if !self.boundaries().contains(test_x, test_y) {
                continue;
            }

            let original_x = x[turbine];
            let original_y = y[turbine];
            x[turbine] = test_x;
            y[turbine] = test_y;

            // The turbines have to keep their minimum distance
            let x_array = Array1::from_vec(x.clone());
            let y_array = Array1::from_vec(y.clone());
            if self.calculate_min_distance(&x_array, &y_array) < min_dist {
                x[turbine] = original_x;
                y[turbine] = original_y;
                continue;
            }

            let test_value = self.calculate_objective(&x_array, &y_array);
            if test_value > current_value {
                current_value = test_value;
            } else {
                x[turbine] = original_x;
                y[turbine] = original_y;
            }
        }

        (x, y, current_value)
    }
    /// Sort the individuals by their objective value (decreasing) and replace
    /// the worst ones with copies of the best ones (Python's
    /// `_evaluate_opt_step`).
    fn evaluate_opt_step(
        &self,
        candidates_x: &mut Vec<Vec<Float>>,
        candidates_y: &mut Vec<Vec<Float>>,
        candidate_values: &mut Vec<Float>,
    ) {
        // Sort the individuals by their objective value (decreasing)
        let mut order: Vec<usize> = (0..candidates_x.len()).collect();
        order.sort_by(|&a, &b| {
            candidate_values[b]
                .partial_cmp(&candidate_values[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let sorted_x: Vec<Vec<Float>> = order.iter().map(|&i| candidates_x[i].clone()).collect();
        let sorted_y: Vec<Vec<Float>> = order.iter().map(|&i| candidates_y[i].clone()).collect();
        let sorted_values: Vec<Float> = order.iter().map(|&i| candidate_values[i]).collect();
        *candidates_x = sorted_x;
        *candidates_y = sorted_y;
        *candidate_values = sorted_values;

        // Replace the worst performing individuals with copies of the best
        // performing ones (the relegation)
        let n = candidates_x.len();
        let relegation = self.config.relegation_number.min(n / 2);
        for k in 0..relegation {
            let dst = n - relegation + k;
            candidates_x[dst] = candidates_x[k].clone();
            candidates_y[dst] = candidates_y[k].clone();
            candidate_values[dst] = candidate_values[k];
        }
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
            fmodel.get_farm_avp().unwrap_or(0.0)
        } else {
            fmodel.get_farm_aep_uniform(8760.0).unwrap_or(0.0)
        }
    }

    fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        let min_dist = self.min_dist();

        // The initial layout and its objective
        let initial_x = self.farm_layout_x();
        let initial_y = self.farm_layout_y();
        self.initial_value = self.calculate_objective(&initial_x, &initial_y);

        // The probability mass function of the step distances
        let pmf = self.distance_pmf();
        let grid_step = self.config.grid_step_size;

        // The seed of the random number generator (Python uses the seeds
        // 23 + i for the initial layouts and 55 + generation + i for the
        // optimization of the individuals)
        let base_seed = self.random_seed.unwrap_or(0) as u32;

        // The initial layouts of the individuals (Python's
        // `_generate_initial_layouts`)
        let mut candidates_x: Vec<Vec<Float>> = Vec::new();
        let mut candidates_y: Vec<Vec<Float>> = Vec::new();
        let mut candidate_values: Vec<Float> = Vec::new();
        for i in 0..self.n_individuals {
            let mut rng = NumpyRandomState::new(base_seed.wrapping_add(23 + i as u32));
            let (x, y) = self.distance_based_init(&mut rng, grid_step);
            let value = self
                .calculate_objective(&Array1::from_vec(x.clone()), &Array1::from_vec(y.clone()));
            candidates_x.push(x);
            candidates_y.push(y);
            candidate_values.push(value);
        }

        // The initial evaluation step (Python's `_evaluate_opt_step` in the
        // constructor)
        self.evaluate_opt_step(&mut candidates_x, &mut candidates_y, &mut candidate_values);

        // The generations of the genetic algorithm (Python's
        // `_run_optimization_generation`)
        let mut iterations = 0;
        for generation in 0..self.config.n_generations {
            for i in 0..self.n_individuals {
                let mut rng = NumpyRandomState::new(
                    base_seed.wrapping_add(55 + generation as u32 + i as u32),
                );
                let (x, y, value) = self.optimize_individual(
                    &mut rng,
                    &candidates_x[i],
                    &candidates_y[i],
                    &pmf,
                    self.config.max_iterations,
                    min_dist,
                );
                iterations += self.config.max_iterations;
                candidates_x[i] = x;
                candidates_y[i] = y;
                candidate_values[i] = value;
            }

            // Sort the individuals by their objective value (decreasing) and
            // replace the worst ones with the best ones (Python's
            // `_evaluate_opt_step`)
            self.evaluate_opt_step(&mut candidates_x, &mut candidates_y, &mut candidate_values);
        }

        // The best individual (the first one after the sorting)
        let best_idx = 0;

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
            iterations,
            improvement_pct,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::floris_config::{
        FarmConfig, FlorisConfig, FlowFieldConfig, SolverConfig, WakeConfig,
    };
    use crate::floris_model::FlorisModel;
    use crate::numpy_random::NumpyRandomState;
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

        let optimizer = LayoutOptimizationRandomSearch::new(&model, boundary, 10).unwrap();
        assert_eq!(optimizer.n_turbines(), 3);
        assert_eq!(optimizer.n_individuals, 10);
    }
}
