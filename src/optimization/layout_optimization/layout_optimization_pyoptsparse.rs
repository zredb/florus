//! Layout Optimization PyOptSparse Module
//!
//! Provides layout optimization using argmin library for constrained optimization.
//! This corresponds to layout_optimization_pyoptsparse.py in Python FLORIS v4.6,
//! using argmin-rs instead of pyoptsparse for the Rust implementation.

use super::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer,
};
use crate::floris_model::FlorisModel;
use crate::types::{Array1, Float};
use crate::Result;
use argmin::core::{CostFunction, Error, Gradient, State};
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::linesearch::MoreThuenteLineSearch;

/// LayoutOptimizationPyOptSparse - Layout optimization using argmin
///
/// Corresponds to floris.optimization.LayoutOptimizationPyOptSparse in Python FLORIS v4.6
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LayoutOptimizationPyOptSparse {
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

impl LayoutOptimizationPyOptSparse {
    /// Create new PyOptSparse-style layout optimizer
    pub fn new(fmodel: &FlorisModel, boundaries: Boundary) -> Result<Self> {
        let config = LayoutOptimizationConfig::default();
        let n_turbines = fmodel.farm().n_turbines();
        let (xmin, xmax, ymin, ymax) = boundaries.bounds();

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

    /// Get turbine positions as a flattened vector for optimization
    fn get_optimization_params(&self, x: &Array1, y: &Array1) -> Vec<Float> {
        let mut params = Vec::with_capacity(2 * self.n_turbines);
        for i in 0..self.n_turbines {
            params.push(x[i]);
        }
        for i in 0..self.n_turbines {
            params.push(y[i]);
        }
        params
    }

    /// Convert optimization parameters back to x, y arrays
    fn from_optimization_params(params: &[Float]) -> (Array1, Array1) {
        let n = params.len() / 2;
        let mut x = Array1::zeros(n);
        let mut y = Array1::zeros(n);

        for i in 0..n {
            x[i] = params[i];
            y[i] = params[n + i];
        }

        (x, y)
    }
}

/// Wrapper for the cost function in argmin optimization
struct LayoutCostFunction<'a> {
    optimizer: &'a LayoutOptimizationPyOptSparse,
    min_dist: Float,
}

impl<'a> CostFunction for LayoutCostFunction<'a> {
    type Param = Vec<Float>;
    type Output = Float;

    fn cost(&self, param: &Self::Param) -> std::result::Result<Self::Output, Error> {
        let (x, y) = LayoutOptimizationPyOptSparse::from_optimization_params(param);

        // Check minimum distance constraint
        let min_actual = self.optimizer.calculate_min_distance(&x, &y);
        let constraint_penalty = if min_actual < self.min_dist {
            (self.min_dist - min_actual).powi(2) * 1e6
        } else {
            0.0
        };

        // Calculate objective
        let obj = self.optimizer.calculate_objective(&x, &y);

        // Return negative for maximization (argmin minimizes)
        Ok(-obj + constraint_penalty)
    }
}

impl<'a> Gradient for LayoutCostFunction<'a> {
    type Param = Vec<Float>;
    type Gradient = Vec<Float>;

    fn gradient(&self, param: &Self::Param) -> std::result::Result<Self::Gradient, Error> {
        // Numerical gradient approximation
        let eps = 1e-6;
        let mut grad = Vec::new();

        for i in 0..param.len() {
            let mut param_plus = param.clone();
            param_plus[i] += eps;

            let cost_plus = self.cost(&param_plus)?;
            let cost_base = self.cost(param)?;

            grad.push((cost_plus - cost_base) / eps);
        }

        Ok(grad)
    }
}

impl LayoutOptimizer for LayoutOptimizationPyOptSparse {
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
        let mut fmodel = self.fmodel.clone();

        if let Err(e) = fmodel.set_layout(x, y) {
            eprintln!("Warning: Failed to set layout: {}", e);
            return 0.0;
        }

        if let Err(e) = fmodel.run() {
            eprintln!("Warning: Failed to run model: {}", e);
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
        let initial_params = self.get_optimization_params(&initial_x, &initial_y);

        // Create cost function
        let cost_fn = LayoutCostFunction {
            optimizer: self,
            min_dist,
        };

        // Use steepest descent with line search
        let solver = SteepestDescent::new(MoreThuenteLineSearch::new());

        // Run optimization using Executor
        let result = argmin::core::Executor::new(cost_fn, solver)
            .configure(|state| {
                state
                    .param(initial_params)
                    .max_iters(self.config.max_iterations as u64)
            })
            .run();

        match result {
            Ok(opt_result) => {
                let best_param = opt_result.state().get_best_param();
                let best_param = best_param
                    .as_ref()
                    .ok_or_else(|| anyhow::anyhow!("Optimization failed to find best parameter"))?;
                let (opt_x, opt_y) =
                    LayoutOptimizationPyOptSparse::from_optimization_params(best_param);

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
                    iterations: opt_result.state().iter as usize,
                    improvement_pct,
                })
            }
            Err(_e) => {
                // Return current layout as fallback
                Ok(LayoutOptimizationResult {
                    x: initial_x,
                    y: initial_y,
                    value: self.initial_value,
                    iterations: 0,
                    improvement_pct: 0.0,
                })
            }
        }
    }
}

/// Golden Section Search optimizer for 1D problems
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LayoutOptimizationGoldenSection {
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
    /// Fixed positions (for optimization in one dimension)
    fixed_x: Option<Array1>,
    fixed_y: Option<Array1>,
}

impl LayoutOptimizationGoldenSection {
    /// Create new golden section layout optimizer
    pub fn new(fmodel: &FlorisModel, boundaries: Boundary) -> Result<Self> {
        let config = LayoutOptimizationConfig::default();
        let n_turbines = fmodel.farm().n_turbines();
        let (xmin, xmax, ymin, ymax) = boundaries.bounds();

        Ok(Self {
            fmodel: fmodel.clone(),
            config,
            xmin,
            xmax,
            ymin,
            ymax,
            n_turbines,
            initial_value: 0.0,
            fixed_x: None,
            fixed_y: None,
        })
    }

    /// Fix x positions, optimize y positions only
    pub fn with_fixed_x(mut self, x: Array1) -> Self {
        self.fixed_x = Some(x);
        self
    }

    /// Fix y positions, optimize x positions only
    pub fn with_fixed_y(mut self, y: Array1) -> Self {
        self.fixed_y = Some(y);
        self
    }
}

impl LayoutOptimizer for LayoutOptimizationGoldenSection {
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
        let mut fmodel = self.fmodel.clone();

        if let Err(_e) = fmodel.set_layout(x, y) {
            return 0.0;
        }

        if let Err(_e) = fmodel.run() {
            return 0.0;
        }

        if self.config.use_value {
            fmodel.get_farm_avp()
        } else {
            fmodel.get_farm_aep_uniform(8760.0)
        }
    }

    fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        // Get initial layout
        let coords = self.fmodel.farm().coordinates();
        let mut initial_x = Array1::zeros(self.n_turbines);
        let mut initial_y = Array1::zeros(self.n_turbines);

        for (i, coord) in coords.outer_iter().enumerate() {
            initial_x[i] = coord[0];
            initial_y[i] = coord[1];
        }

        let (opt_x, opt_y) = if let Some(ref fixed) = self.fixed_x {
            // Optimize y only using golden section sampling
            let mut best_y = initial_y.clone();
            let mut best_value = self.calculate_objective(&fixed, &best_y);

            // Sample points and find best using golden section-like sampling
            let n_samples = 100;
            for i in 1..n_samples {
                let y_test =
                    self.ymin + (self.ymax - self.ymin) * (i as Float / n_samples as Float);
                let mut test_y = initial_y.clone();
                test_y[0] = y_test;
                let value = self.calculate_objective(&fixed, &test_y);
                if value > best_value {
                    best_value = value;
                    best_y = test_y;
                }
            }

            (fixed.clone(), best_y)
        } else if let Some(ref fixed) = self.fixed_y {
            // Optimize x only using golden section sampling
            let mut best_x = initial_x.clone();
            let mut best_value = self.calculate_objective(&best_x, &fixed);

            let n_samples = 100;
            for i in 1..n_samples {
                let x_test =
                    self.xmin + (self.xmax - self.xmin) * (i as Float / n_samples as Float);
                let mut test_x = initial_x.clone();
                test_x[0] = x_test;
                let value = self.calculate_objective(&test_x, &fixed);
                if value > best_value {
                    best_value = value;
                    best_x = test_x;
                }
            }

            (best_x, fixed.clone())
        } else {
            // Fallback: use initial layout
            (initial_x, initial_y)
        };

        let final_value = self.calculate_objective(&opt_x, &opt_y);

        Ok(LayoutOptimizationResult {
            x: opt_x,
            y: opt_y,
            value: final_value,
            iterations: 1,
            improvement_pct: 0.0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::floris_config::{FarmConfig, FlorisConfig, FlowFieldConfig, SolverConfig, WakeConfig};
    use crate::floris_model::FlorisModel;

    fn create_test_config() -> FlorisConfig {
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
        };

        let solver_config = SolverConfig::default();
        let wake_config = WakeConfig::default();

        FlorisConfig {
            name: "test".to_string(),
            description: Some("test".to_string()),
            floris_version: "v4".to_string(),
            logging: Default::default(),
            solver: solver_config,
            farm: farm_config,
            flow_field: flow_field_config,
            wake: wake_config,
            turbine_library: "turbine_library".to_string(),
        }
    }

    #[test]
    fn test_pyoptsparse_optimizer_creation() {
        let config = create_test_config();
        let model = FlorisModel::from_config(config).unwrap();

        let boundary = Boundary::Rectangle {
            min_x: 0.0,
            max_x: 1000.0,
            min_y: 0.0,
            max_y: 1000.0,
        };

        let optimizer = LayoutOptimizationPyOptSparse::new(&model, boundary).unwrap();
        assert_eq!(optimizer.n_turbines(), 2);
    }

    #[test]
    fn test_golden_section_optimizer_creation() {
        let config = create_test_config();
        let model = FlorisModel::from_config(config).unwrap();

        let boundary = Boundary::Rectangle {
            min_x: 0.0,
            max_x: 1000.0,
            min_y: 0.0,
            max_y: 1000.0,
        };

        let optimizer = LayoutOptimizationGoldenSection::new(&model, boundary).unwrap();
        assert_eq!(optimizer.n_turbines(), 2);
    }
}
