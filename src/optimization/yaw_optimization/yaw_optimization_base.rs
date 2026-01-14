use crate::types::{Array1, Array2, Array4, Float};
use crate::floris_model::FlorisModel;
use crate::core::grid::GridBase;

/// Result of yaw optimization
#[derive(Debug, Clone)]
pub struct YawOptimizationResult {
    pub yaw_angles: Array2,
    pub powers: Array2,
    pub baseline_power: Float,
    pub optimized_power: Float,
    pub power_improvement: Float,
    pub improvement_percentage: Float,
}

impl Default for YawOptimizationResult {
    fn default() -> Self {
        Self {
            yaw_angles: Array2::zeros((0, 0)),
            powers: Array2::zeros((0, 0)),
            baseline_power: 0.0,
            optimized_power: 0.0,
            power_improvement: 0.0,
            improvement_percentage: 0.0,
        }
    }
}

/// Configuration for yaw optimization
#[derive(Debug, Clone)]
pub struct YawOptimizationConfig {
    pub minimum_yaw_angle: Float,
    pub maximum_yaw_angle: Float,
    pub yaw_angles_baseline: Option<Array2>,
    pub turbine_weights: Option<Vec<Float>>,
    pub exclude_downstream_turbines: bool,
    pub verify_convergence: bool,
}

impl Default for YawOptimizationConfig {
    fn default() -> Self {
        Self {
            minimum_yaw_angle: 0.0,
            maximum_yaw_angle: 25.0,
            yaw_angles_baseline: None,
            turbine_weights: None,
            exclude_downstream_turbines: true,
            verify_convergence: false,
        }
    }
}

/// Yaw angle bounds for optimization
#[derive(Debug, Clone)]
pub struct YawAngleBounds {
    pub min_yaw: Float,
    pub max_yaw: Float,
}

impl Default for YawAngleBounds {
    fn default() -> Self {
        Self {
            min_yaw: -45.0,
            max_yaw: 45.0,
        }
    }
}

impl YawAngleBounds {
    pub fn new(min_yaw: Float, max_yaw: Float) -> Self {
        Self { min_yaw, max_yaw }
    }
}

/// Trapezoid bounds for geometric yaw optimization
#[derive(Debug, Clone, Default)]
pub struct TrapezoidBounds {
    pub left_x: Float,
    pub top_left_y: Float,
    pub right_x: Float,
    pub top_right_y: Float,
    pub top_left_yaw_upper: Float,
    pub top_right_yaw_upper: Float,
    pub bottom_left_yaw_upper: Float,
    pub bottom_right_yaw_upper: Float,
    pub top_left_yaw_lower: Float,
    pub top_right_yaw_lower: Float,
    pub bottom_left_yaw_lower: Float,
    pub bottom_right_yaw_lower: Float,
}

impl TrapezoidBounds {
    pub fn new() -> Self {
        Self {
            left_x: 0.0,
            top_left_y: 1.0,
            right_x: 25.0,
            top_right_y: 1.0,
            top_left_yaw_upper: 30.0,
            top_right_yaw_upper: 0.0,
            bottom_left_yaw_upper: 30.0,
            bottom_right_yaw_upper: 0.0,
            top_left_yaw_lower: -30.0,
            top_right_yaw_lower: 0.0,
            bottom_left_yaw_lower: -30.0,
            bottom_right_yaw_lower: 0.0,
        }
    }
}

/// Yaw Optimization trait - common interface for all yaw optimization methods
pub trait YawOptimization {
    fn optimize(
        &mut self,
        fmodel: &mut FlorisModel,
        config: Option<YawOptimizationConfig>,
    ) -> anyhow::Result<YawOptimizationResult>;
}

/// Helper functions for variable handling
pub mod variable_handling {
    use super::*;

    /// Expand scalar value to full shape (n_findex, n_turbines)
    pub fn unpack_scalar(value: Float, n_findex: usize, n_turbines: usize) -> Array2 {
        Array2::from_elem((n_findex, n_turbines), value)
    }

    /// Unpack array to target shape
    pub fn unpack_array(
        value: &Array2,
        n_findex: usize,
        n_turbines: usize,
    ) -> anyhow::Result<Array2> {
        if value.shape() == [n_findex, n_turbines] {
            Ok(value.clone())
        } else {
            let result = value.broadcast((n_findex, n_turbines))
                .ok_or_else(|| anyhow::anyhow!("Cannot broadcast array from {:?} to ({}, {})", value.shape(), n_findex, n_turbines))?;
            Ok(result.to_owned())
        }
    }
}

/// Helper functions for reducing the control problem
pub mod control_problem {
    use super::*;

    /// Remove downstream turbines from optimization variables
    pub fn exclude_downstream(
        turbine_indices: &[usize],
        downstream_indices: &[usize],
    ) -> Vec<usize> {
        turbine_indices
            .iter()
            .filter(|&ti| !downstream_indices.contains(&ti))
            .copied()
            .collect()
    }

    /// Normalize optimization variables for better numerical properties
    pub fn normalize_yaw_angles(yaw_angles: &Array2) -> Array2 {
        let n_findex = yaw_angles.shape()[0];
        let n_turbines = yaw_angles.shape()[1];

        if n_turbines == 0 {
            return yaw_angles.clone();
        }

        let mut mean_yaw = Array2::zeros((1, n_turbines));
        for ti in 0..n_turbines {
            let sum: Float = (0..n_findex).map(|fi| yaw_angles[[fi, ti]]).sum();
            mean_yaw[[0, ti]] = sum / n_findex as Float;
        }

        let mut normalized = yaw_angles.clone();
        for ti in 0..n_turbines {
            let mean = mean_yaw[[0, ti]];
            for fi in 0..n_findex {
                normalized[[fi, ti]] -= mean;
            }
        }

        normalized
    }
}

/// Simple yaw optimization using coordinate descent
pub fn simple_yaw_optimization(
    fmodel: &mut FlorisModel,
    config: Option<YawOptimizationConfig>,
) -> anyhow::Result<YawOptimizationResult> {
    let config = config.unwrap_or_default();

    if !fmodel.state.initialized {
        anyhow::bail!("FlorisModel must be run before optimization");
    }

    let n_findex = fmodel.flow_field.n_findex;
    let n_turbines = fmodel.farm.n_turbines();

    let baseline_yaw_angles = config.yaw_angles_baseline.clone()
        .unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));
    fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
    fmodel.run()?;
    let baseline_power = fmodel.get_farm_power().sum();

    let mut yaw_angles = baseline_yaw_angles.clone();

    let downstream_indices = if config.exclude_downstream_turbines {
        super::derive_downstream_turbines(fmodel, 0.3, false)?
    } else {
        vec![]
    };

    let active_turbines: Vec<usize> = (0..n_turbines)
        .filter(|ti| !downstream_indices.contains(ti))
        .collect();

    let turbine_weights = config.turbine_weights.clone()
        .unwrap_or_else(|| vec![1.0; n_turbines]);

    let min_yaw = config.minimum_yaw_angle;
    let max_yaw = config.maximum_yaw_angle;
    let convergence_tolerance = 1e-6;
    let max_optimization_iterations = 10;

    for fi in 0..n_findex {
        let mut iteration = 0;
        let mut prev_power = fmodel.get_farm_power()[fi];
        let mut converged = false;

        while iteration < max_optimization_iterations && !converged {
            let mut max_change: Float = 0.0;

            for &ti in &active_turbines {
                let current_yaw = yaw_angles[[fi, ti]];
                let weight = turbine_weights[ti];

                let delta_yaw = 1.0;
                let yaw_plus = (current_yaw + delta_yaw).clamp(min_yaw, max_yaw);
                let yaw_minus = (current_yaw - delta_yaw).clamp(min_yaw, max_yaw);

                yaw_angles[[fi, ti]] = yaw_plus;
                fmodel.set_yaw_angles(yaw_angles.clone())?;
                fmodel.run()?;
                let power_plus = fmodel.get_farm_power()[fi];

                yaw_angles[[fi, ti]] = yaw_minus;
                fmodel.set_yaw_angles(yaw_angles.clone())?;
                fmodel.run()?;
                let power_minus = fmodel.get_farm_power()[fi];

                let gradient = (power_plus - power_minus) / (2.0 * delta_yaw);
                let learning_rate = 2.0 * weight;
                let step = gradient * learning_rate;

                let new_yaw = (current_yaw - step).clamp(min_yaw, max_yaw);
                let change = (new_yaw - current_yaw).abs();
                max_change = max_change.max(change);

                yaw_angles[[fi, ti]] = new_yaw;
            }

            fmodel.set_yaw_angles(yaw_angles.clone())?;
            fmodel.run()?;
            let current_power = fmodel.get_farm_power()[fi];

            if (current_power - prev_power).abs() < convergence_tolerance * prev_power.abs() {
                converged = true;
            }

            prev_power = current_power;
            iteration += 1;
        }
    }

    fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
    fmodel.run()?;
    let baseline_power_verify = fmodel.get_farm_power().sum();

    fmodel.set_yaw_angles(yaw_angles.clone())?;
    fmodel.run()?;
    let optimized_power = fmodel.get_farm_power().sum();

    let power_improvement = optimized_power - baseline_power_verify;
    let improvement_percentage = if baseline_power_verify > 0.0 {
        (power_improvement / baseline_power_verify) * 100.0
    } else {
        0.0
    };

    let turbine_powers = fmodel.get_turbine_powers();

    Ok(YawOptimizationResult {
        yaw_angles,
        powers: turbine_powers,
        baseline_power: baseline_power_verify,
        optimized_power,
        power_improvement,
        improvement_percentage,
    })
}

/// Geometric yaw optimization - standalone function
pub fn geometric_yaw(
    turbine_x: &[Float],
    turbine_y: &[Float],
    wind_direction: Float,
    rotor_diameter: Float,
    bounds: Option<TrapezoidBounds>,
) -> Vec<Float> {
    let bounds = bounds.unwrap_or_default();
    let n_turbines = turbine_x.len();
    let mut yaw_angles = Vec::with_capacity(n_turbines);

    let wind_direction_rad = wind_direction.to_radians();
    let cos_wd = wind_direction_rad.cos();
    let sin_wd = wind_direction_rad.sin();

    // Rotate coordinates
    let rotated_x: Vec<Float> = turbine_x
        .iter()
        .zip(turbine_y.iter())
        .map(|(&x, &y)| x * cos_wd + y * sin_wd)
        .collect();
    let rotated_y: Vec<Float> = turbine_x
        .iter()
        .zip(turbine_y.iter())
        .map(|(&x, &y)| -x * sin_wd + y * cos_wd)
        .collect();

    let min_x = rotated_x.iter().fold(Float::INFINITY, |m, &x| x.min(m));
    let max_x = rotated_x.iter().fold(Float::NEG_INFINITY, |m, &x| x.max(m));

    for i in 0..n_turbines {
        let x_norm = if max_x > min_x {
            (rotated_x[i] - min_x) / (max_x - min_x)
        } else {
            0.5
        };

        let d_norm = rotor_diameter;
        let x_d = rotated_x[i] / d_norm;

        let upper_yaw = if x_d <= bounds.left_x {
            bounds.top_left_yaw_upper
        } else if x_d >= bounds.right_x {
            bounds.top_right_yaw_upper
        } else {
            let t = (x_d - bounds.left_x) / (bounds.right_x - bounds.left_x);
            bounds.top_left_yaw_upper + t * (bounds.top_right_yaw_upper - bounds.top_left_yaw_upper)
        };

        let lower_yaw = if x_d <= bounds.left_x {
            bounds.top_left_yaw_lower
        } else if x_d >= bounds.right_x {
            bounds.top_right_yaw_lower
        } else {
            let t = (x_d - bounds.left_x) / (bounds.right_x - bounds.left_x);
            bounds.top_left_yaw_lower + t * (bounds.top_right_yaw_lower - bounds.top_left_yaw_lower)
        };

        let yaw = (upper_yaw + lower_yaw) / 2.0;
        yaw_angles.push(yaw.clamp(lower_yaw, upper_yaw));
    }

    yaw_angles
}
