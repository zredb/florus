//! Layout Optimization Scipy Module
//!
//! Provides scipy-style layout optimization with the SLSQP solver, which is a
//! port of scipy's solver (see [`super::slsqp`]).
//! This module corresponds to layout_optimization_scipy.py in Python FLORIS v4.6.

use super::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer,
};
use super::slsqp::minimize_slsqp;
use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2, Float};
use crate::Result;

/// The yaw angles of the geometric yaw optimization for the current layout
/// (Python's `_get_geoyaw_angles`).
///
/// The layout optimization uses the geometric yaw optimizer with the yaw
/// bounds of -30 and 30 degrees (Python's `layout_optimization_base.py:92-98`).
pub fn geometric_yaw_angles(fmodel: &mut FlorisModel) -> Result<Array2> {
    use crate::optimization::yaw_optimization::{
        YawOptimization, YawOptimizationConfig, YawOptimizationGeometric,
    };

    let config = YawOptimizationConfig {
        minimum_yaw_angle: -30.0,
        maximum_yaw_angle: 30.0,
        ..Default::default()
    };
    let mut optimizer = YawOptimizationGeometric::new();
    let result = optimizer.optimize(fmodel, Some(config))?;
    Ok(result.yaw_angles)
}

/// LayoutOptimizationScipy - Layout optimization using scipy's SLSQP solver
///
/// The optimization variables are the normalized coordinates of the turbines
/// (the x coordinates followed by the y coordinates), the objective is the
/// negative annual energy (or value) production normalized with the value of
/// the initial layout, and the constraints are the Kolmogorov-Smirnov
/// constraint of the minimum distance between the turbines and the distance of
/// every turbine to the boundary. The derivatives are the forward differences
/// with a step of 0.01 in the normalized coordinates.
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

    /// The physical coordinates of the normalized optimization variables (the x
    /// coordinates followed by the y coordinates).
    fn denormalize_variables(&self, x: &[Float]) -> (Vec<Float>, Vec<Float>) {
        let n = self.n_turbines;
        let px = x[..n]
            .iter()
            .map(|&v| self.denormalize(v, self.xmin, self.xmax))
            .collect();
        let py = x[n..]
            .iter()
            .map(|&v| self.denormalize(v, self.ymin, self.ymax))
            .collect();
        (px, py)
    }

    /// The Kolmogorov-Smirnov constraint of the minimum distance between the
    /// turbines (Python's `_space_constraint`).
    ///
    /// The constraint is satisfied when it is not positive.
    fn space_constraint(&self, x: &[Float], y: &[Float], min_dist: Float) -> Float {
        // The distance of every turbine to the closest other turbine
        let mut g = vec![0.0; x.len()];
        for i in 0..x.len() {
            let mut closest = Float::INFINITY;
            for j in 0..x.len() {
                if i == j {
                    continue;
                }
                let dx = x[i] - x[j];
                let dy = y[i] - y[j];
                let distance = (dx * dx + dy * dy).sqrt();
                if distance < closest {
                    closest = distance;
                }
            }
            g[i] = 1.0 - closest / min_dist;
        }

        // The Kolmogorov-Smirnov aggregation (OpenMDAO's KSComp)
        let rho = 500.0;
        let g_max = g.iter().copied().fold(Float::NEG_INFINITY, Float::max);
        let sum: Float = g.iter().map(|&value| (rho * (value - g_max)).exp()).sum();
        let ks_constraint = g_max + (1.0 / rho) * sum.ln();

        // The constraint is satisfied when it is not positive
        -ks_constraint
    }
}

/// The gradient of a function by forward differences (scipy's
/// `approx_derivative` with the '2-point' method and an absolute step).
///
/// The step is a fixed `+step`: scipy computes the gradient of the objective
/// WITHOUT the bounds of the variables (`_slsqp_py.py:361`), so the step is
/// never flipped (the perturbed point may lie outside of the unit interval).
fn approx_gradient<F>(fun: &mut F, x: &[Float], step: Float) -> Vec<Float>
where
    F: FnMut(&[Float]) -> Float,
{
    let f0 = fun(x);
    let n = x.len();
    let mut gradient = vec![0.0; n];
    let mut perturbed = x.to_vec();

    for i in 0..n {
        let original = perturbed[i];
        perturbed[i] = original + step;
        let f1 = fun(&perturbed);
        let dx = perturbed[i] - original;
        perturbed[i] = original;
        gradient[i] = (f1 - f0) / dx;
    }

    gradient
}

/// The Jacobian of a function by forward differences (scipy's
/// `approx_derivative` with the '2-point' method and an absolute step).
///
/// scipy computes the Jacobian of the constraints WITH the bounds of the
/// variables (`_slsqp_py.py:361`), so the step is adjusted to the bounds with
/// the '1-sided' scheme (scipy's `_adjust_scheme_to_bounds`).
fn approx_jacobian<F>(fun: &mut F, x: &[Float], step: Float) -> Vec<Vec<Float>>
where
    F: FnMut(&[Float]) -> Vec<Float>,
{
    let f0 = fun(x);
    let n = x.len();
    let m = f0.len();
    let mut jacobian = vec![vec![0.0; n]; m];
    let mut perturbed = x.to_vec();

    for i in 0..n {
        let h = adjusted_step(x[i], step);
        let original = perturbed[i];
        perturbed[i] = original + h;
        let f1 = fun(&perturbed);
        let dx = perturbed[i] - original;
        perturbed[i] = original;
        for k in 0..m {
            jacobian[k][i] = (f1[k] - f0[k]) / dx;
        }
    }

    jacobian
}

/// The step of the finite differences of a variable of the unit interval.
///
/// The step is flipped when the perturbed point would leave the interval
/// (scipy's `_adjust_scheme_to_bounds` with the '1-sided' scheme). When the
/// step does not fit in either direction, it is shrunk to the distance to the
/// closest bound.
fn adjusted_step(value: Float, step: Float) -> Float {
    let lower_dist = value; // the lower bound of the variables is 0
    let upper_dist = 1.0 - value;
    if value + step > 1.0 || value + step < 0.0 {
        // The step leaves the interval: it is flipped when it still fits in
        // the other direction, and shrunk to the closest bound otherwise.
        if step.abs() <= lower_dist.max(upper_dist) {
            -step
        } else if upper_dist >= lower_dist {
            upper_dist
        } else {
            -lower_dist
        }
    } else {
        step
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

        // The value of the initial layout (Python's `initial_AEP_or_AVP`,
        // which normalizes the objective of the optimization)
        let initial_value = self.calculate_objective(&initial_x, &initial_y);
        self.initial_value = initial_value;

        // The optimization variables are the normalized coordinates of the
        // turbines: the x coordinates followed by the y coordinates (Python's
        // `x0`).
        let n = self.n_turbines;
        let mut x0: Vec<Float> = Vec::with_capacity(2 * n);
        x0.extend((0..n).map(|i| self.normalize(initial_x[i], self.xmin, self.xmax)));
        x0.extend((0..n).map(|i| self.normalize(initial_y[i], self.ymin, self.ymax)));

        // The bounds of the variables (Python's `_set_opt_bounds`)
        let xl = vec![0.0; 2 * n];
        let xu = vec![1.0; 2 * n];

        // The step of the finite differences of the derivatives (Python's
        // `eps` option)
        let step = self.config.eps;

        // The objective: the negative annual energy (or value) production,
        // normalized with the value of the initial layout (Python's
        // `_obj_func`).
        let mut objective = |x: &[Float]| -> Float {
            let (px, py) = self.denormalize_variables(x);
            -self.calculate_objective(&Array1::from_vec(px), &Array1::from_vec(py)) / initial_value
        };
        let mut objective_gradient = |x: &[Float]| -> Vec<Float> {
            approx_gradient(
                &mut |point: &[Float]| {
                    let (px, py) = self.denormalize_variables(point);
                    -self.calculate_objective(&Array1::from_vec(px), &Array1::from_vec(py))
                        / initial_value
                },
                x,
                step,
            )
        };

        // The constraints: the Kolmogorov-Smirnov constraint of the minimum
        // distance between the turbines (Python's `_space_constraint`) and
        // the distance of every turbine to the boundary (Python's
        // `_distance_from_boundaries`).
        let mut constraints = |x: &[Float]| -> Vec<Float> {
            let (px, py) = self.denormalize_variables(x);
            let mut values = vec![self.space_constraint(&px, &py, min_dist)];
            values.extend((0..n).map(|i| self.boundaries().signed_distance(px[i], py[i])));
            values
        };
        let mut constraints_jacobian = |x: &[Float]| -> Vec<Vec<Float>> {
            approx_jacobian(
                &mut |point: &[Float]| {
                    let (px, py) = self.denormalize_variables(point);
                    let mut values = vec![self.space_constraint(&px, &py, min_dist)];
                    values.extend((0..n).map(|i| self.boundaries().signed_distance(px[i], py[i])));
                    values
                },
                x,
                step,
            )
        };

        // The optimization (Python's `scipy.optimize.minimize` with the
        // SLSQP method)
        let result = minimize_slsqp(
            &x0,
            &mut objective,
            &mut objective_gradient,
            &mut constraints,
            &mut constraints_jacobian,
            &xl,
            &xu,
            self.config.max_iterations,
            self.config.tolerance,
        );

        // Convert back to physical coordinates
        let (opt_x, opt_y) = self.denormalize_variables(&result.x);
        let opt_x = Array1::from_vec(opt_x);
        let opt_y = Array1::from_vec(opt_y);

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
            iterations: result.nit,
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
