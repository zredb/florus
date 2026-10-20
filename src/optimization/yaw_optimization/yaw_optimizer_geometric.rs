//! Geometric yaw optimizer — port of Python
//! `floris.optimization.yaw_optimization.yaw_optimizer_geometric`.
//!
//! Provides a rough estimate of optimal yaw angles based purely on the wind
//! farm geometry (Stanley et al. 2023): for each wind direction, the layout
//! is rotated into the wind frame, each turbine's distance to its nearest
//! downstream waked turbine is computed (Jensen-like spread, `spread = 0.1`),
//! and the yaw angle follows from the turbine's position relative to a
//! trapezoid of waked turbines.
//!
//! The yaw bounds come from the optimizer configuration
//! (`minimum_yaw_angle` / `maximum_yaw_angle`), as in Python.

use super::{YawOptimization, YawOptimizationConfig, YawOptimizationResult};
use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2, Float};
use crate::utilities::rotate_coordinates_rel_west;

/// Jensen-like wake spread used to decide which turbines are "waked"
/// (Python `_process_layout`, spread=0.1).
const SPREAD: Float = 0.1;

/// Geometric yaw optimizer
///
/// Computes yaw angles based on farm geometry and wind direction.
#[derive(Debug, Clone)]
pub struct YawOptimizationGeometric {
    config: YawOptimizationConfig,
}

impl YawOptimizationGeometric {
    pub fn new() -> Self {
        Self {
            config: YawOptimizationConfig::default(),
        }
    }

    /// Python `_process_layout` (yaw_optimizer_geometric.py:146-184): for each
    /// turbine, the distance (in rotor diameters) to the nearest downstream
    /// turbine inside a Jensen-like wake.
    fn process_layout(
        turbine_x: &[Float],
        turbine_y: &[Float],
        rotor_diameter: Float,
    ) -> (Vec<Float>, Vec<Float>) {
        let n = turbine_x.len();
        let mut dx = vec![f64::INFINITY; n];
        let mut dy = vec![0.0; n];

        for i in 0..n {
            let mut best = f64::INFINITY;
            let mut best_y = 0.0;
            for j in 0..n {
                // x_dists = turbine_x[j] - turbine_x[i] (downstream is positive)
                let x_dist = turbine_x[j] - turbine_x[i];
                let y_dist = turbine_y[j] - turbine_y[i];
                if x_dist <= 0.0 {
                    continue; // upstream or at the same location
                }
                // Check within Jensen model spread
                if y_dist.abs() < SPREAD * x_dist + rotor_diameter / 2.0 && x_dist < best {
                    best = x_dist;
                    best_y = y_dist;
                }
            }
            dx[i] = best / rotor_diameter;
            dy[i] = best_y / rotor_diameter;
        }

        (dx, dy)
    }

    /// Python `_get_yaw_angles` (yaw_optimizer_geometric.py:187-243).
    #[allow(clippy::too_many_arguments)]
    fn get_yaw_angles(
        x: Float,
        y: Float,
        left_x: Float,
        top_left_y: Float,
        right_x: Float,
        top_right_y: Float,
        top_left_yaw_upper: Float,
        top_right_yaw_upper: Float,
        bottom_left_yaw_upper: Float,
        bottom_right_yaw_upper: Float,
        top_left_yaw_lower: Float,
        top_right_yaw_lower: Float,
        bottom_left_yaw_lower: Float,
        bottom_right_yaw_lower: Float,
    ) -> Float {
        let dx = (x - left_x) / (right_x - left_x);
        if dx >= 1.0 {
            return 0.0;
        }
        let edge_y = top_left_y + (top_right_y - top_left_y) * dx;
        if y.abs() > edge_y {
            0.0
        } else if y >= -0.01 {
            // Tolerance to handle numerical issues
            let top_yaw = top_left_yaw_upper + (top_right_yaw_upper - top_left_yaw_upper) * dx;
            let bottom_yaw =
                bottom_left_yaw_upper + (bottom_right_yaw_upper - bottom_left_yaw_upper) * dx;
            bottom_yaw + (top_yaw - bottom_yaw) * y.abs() / edge_y
        } else {
            let top_yaw = top_left_yaw_lower + (top_right_yaw_lower - top_left_yaw_lower) * dx;
            let bottom_yaw =
                bottom_left_yaw_lower + (bottom_right_yaw_lower - bottom_left_yaw_lower) * dx;
            bottom_yaw + (top_yaw - bottom_yaw) * y.abs() / edge_y
        }
    }

    /// Python `geometric_yaw` (yaw_optimizer_geometric.py:69-144).
    fn geometric_yaw(
        turbine_x: &[Float],
        turbine_y: &[Float],
        wind_direction: Float,
        rotor_diameter: Float,
        min_yaw: Float,
        max_yaw: Float,
    ) -> Vec<Float> {
        let n = turbine_x.len();
        let mut coordinates = Array2::zeros((n, 3));
        for i in 0..n {
            coordinates[[i, 0]] = turbine_x[i];
            coordinates[[i, 1]] = turbine_y[i];
        }

        let wd = Array1::from_vec(vec![wind_direction]);
        let (rotated_x, rotated_y, _, _, _) =
            rotate_coordinates_rel_west(&wd, &coordinates).expect("rotation failed");
        let rx: Vec<Float> = rotated_x.row(0).to_vec();
        let ry: Vec<Float> = rotated_y.row(0).to_vec();

        let (dx, dy) = Self::process_layout(&rx, &ry, rotor_diameter);

        (0..n)
            .map(|i| {
                // Python passes the configuration's yaw bounds as the
                // trapezoid's upper/lower yaw corners.
                Self::get_yaw_angles(
                    dx[i], dy[i], 0.0,     // left_x
                    1.0,     // top_left_y
                    25.0,    // right_x
                    1.0,     // top_right_y
                    max_yaw, // top_left_yaw_upper
                    0.0,     // top_right_yaw_upper
                    max_yaw, // bottom_left_yaw_upper
                    0.0,     // bottom_right_yaw_upper
                    min_yaw, // top_left_yaw_lower
                    0.0,     // top_right_yaw_lower
                    min_yaw, // bottom_left_yaw_lower
                    0.0,     // bottom_right_yaw_lower
                )
            })
            .collect()
    }
}

impl Default for YawOptimizationGeometric {
    fn default() -> Self {
        Self::new()
    }
}

impl YawOptimization for YawOptimizationGeometric {
    fn optimize(
        &mut self,
        fmodel: &mut FlorisModel,
        config: Option<YawOptimizationConfig>,
    ) -> anyhow::Result<YawOptimizationResult> {
        if let Some(cfg) = config {
            self.config = cfg;
        }

        if !fmodel.state().initialized {
            anyhow::bail!("FlorisModel must be run before optimization");
        }

        let n_findex = fmodel.flow_field().n_findex;
        let n_turbines = fmodel.farm().n_turbines();
        let rotor_diameter = fmodel.farm().rotor_diameters[0];

        let layout_x = fmodel.farm().layout_x.clone();
        let layout_y = fmodel.farm().layout_y.clone();
        let wind_directions = fmodel.flow_field().wind_directions.clone();

        // Active turbines: power setpoints above the disabled threshold
        // (Python yaw_optimizer_geometric.py:50-61)
        let power_setpoints = fmodel.farm().power_setpoints.clone();
        let disabled_threshold = crate::core::turbines::operation_models::POWER_SETPOINT_DISABLED;

        let mut yaw_angles = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            let active: Vec<usize> = (0..n_turbines)
                .filter(|&ti| {
                    power_setpoints.shape()[0] == 1
                        || power_setpoints[[fi, ti]] > disabled_threshold
                })
                .collect();

            let xs: Vec<Float> = active.iter().map(|&t| layout_x[t]).collect();
            let ys: Vec<Float> = active.iter().map(|&t| layout_y[t]).collect();
            let yaws = Self::geometric_yaw(
                &xs,
                &ys,
                wind_directions[fi],
                rotor_diameter,
                self.config.minimum_yaw_angle,
                self.config.maximum_yaw_angle,
            );

            for (k, &t) in active.iter().enumerate() {
                yaw_angles[[fi, t]] = yaws[k];
            }
        }

        // Python: calc_baseline_power=False for the geometric optimizer, so
        // the reported baseline is the baseline-yaw farm power and the
        // optimized power is evaluated with the geometric yaw angles.
        let baseline_yaw_angles = self
            .config
            .yaw_angles_baseline
            .clone()
            .unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));

        fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
        fmodel.run()?;
        let farm_powers_baseline = fmodel.get_farm_power()?;
        let baseline_power: Float = farm_powers_baseline.sum();

        fmodel.set_yaw_angles(yaw_angles.clone())?;
        fmodel.run()?;
        let farm_powers_opt = fmodel.get_farm_power()?;
        let optimized_power: Float = farm_powers_opt.sum();

        let power_improvement = optimized_power - baseline_power;
        let improvement_percentage = if baseline_power > 0.0 {
            (power_improvement / baseline_power) * 100.0
        } else {
            0.0
        };

        let turbine_powers = fmodel.get_turbine_powers()?;

        Ok(YawOptimizationResult {
            yaw_angles,
            powers: turbine_powers,
            baseline_power,
            optimized_power,
            power_improvement,
            improvement_percentage,
            farm_powers_opt,
            farm_powers_baseline,
        })
    }
}
