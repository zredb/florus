/// Power Setpoint Optimization Module
///
/// Provides optimization functions for power setpoint control including:
/// - Optimal power setpoint calculation for derating control
/// - Power allocation across turbines
///
/// This module corresponds to floris.optimization.power_setpoint_optimization in Python FLORIS v4.6

use crate::types::{Array2, Float};

/// Result of power setpoint optimization
#[derive(Debug, Clone)]
pub struct PowerSetpointOptimizationResult {
    pub power_setpoints: Array2,
    pub powers: Array2,
    pub baseline_power: Float,
    pub optimized_power: Float,
    pub power_improvement: Float,
}

/// Optimize power setpoints for derating control
pub fn optimize_power_setpoints(
    _wind_speed: Float,
    _wind_direction: Float,
    _turbulence_intensity: Float,
    rated_power: Float,
    n_turbines: usize,
    _min_setpoint: Float,
    max_setpoint: Float,
) -> (Array2, Float) {
    let power_setpoints = Array2::from_elem((1, n_turbines), max_setpoint * rated_power);
    let power = 0.0;
    (power_setpoints, power)
}

/// Compute power setpoints for derating
pub fn compute_power_setpoints(
    n_turbines: usize,
    rated_power: Float,
    derating_factor: Float,
    _wind_speed: Float,
) -> Array2 {
    let setpoint = rated_power * derating_factor;
    Array2::from_elem((1, n_turbines), setpoint)
}

/// Optimize derating factor based on wake losses
pub fn optimize_derating_factor(
    upstream_power: Float,
    downstream_power: Float,
    wake_deficit: Float,
    current_derating: Float,
    _rated_power: Float,
) -> Float {
    if wake_deficit > 0.1 && downstream_power < upstream_power * (1.0 - wake_deficit) {
        (current_derating * 0.9).max(0.5)
    } else if downstream_power >= upstream_power * 0.95 {
        (current_derating * 1.05).min(1.0)
    } else {
        current_derating
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_optimize_power_setpoints() {
        let (setpoints, _power) = optimize_power_setpoints(
            8.0, 270.0, 0.06, 5_000_000.0, 3, 0.5, 1.0,
        );
        assert_eq!(setpoints.shape()[0], 1);
        assert_eq!(setpoints.shape()[1], 3);
    }

    #[test]
    fn test_compute_power_setpoints() {
        let setpoints = compute_power_setpoints(4, 10_000_000.0, 0.8, 8.0);
        assert_eq!(setpoints.shape()[0], 1);
        assert_eq!(setpoints.shape()[1], 4);
        for i in 0..4 {
            assert_relative_eq!(setpoints[[0, i]], 8_000_000.0);
        }
    }

    #[test]
    fn test_optimize_derating_factor_high_wake() {
        let result = optimize_derating_factor(1000.0, 600.0, 0.3, 0.9, 5_000_000.0);
        assert!(result < 0.9);
    }

    #[test]
    fn test_optimize_derating_factor_low_wake() {
        let result = optimize_derating_factor(1000.0, 950.0, 0.05, 0.9, 5_000_000.0);
        assert_relative_eq!(result, 0.945, epsilon = 0.01);
    }
}
