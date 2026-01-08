/// Derating Optimization Module
///
/// Provides optimization functions for turbine derating control including:
/// - Optimal derating factor calculation
/// - Wake-aware derating strategies
///
/// This module corresponds to floris.optimization.derating in Python FLORIS v4.6

use crate::types::Float;

/// Calculate optimal derating factor for wake-affected turbines
pub fn optimize_derating(
    upstream_power: Float,
    downstream_power: Float,
    wake_deficit: Float,
    derating_factor: Float,
) -> Float {
    if wake_deficit > 0.1 && downstream_power < upstream_power * (1.0 - wake_deficit) {
        (derating_factor * 0.9).max(0.5)
    } else {
        derating_factor
    }
}

/// Optimize derating for maximum farm power
pub fn optimize_derating_for_farm_power(
    upstream_rated_power: Float,
    downstream_rated_power: Float,
    wake_deficit: Float,
    cosine_loss_exponent: Float,
) -> Float {
    let mut best_derating = 1.0;
    let mut best_power = 0.0;

    for d in (50..=100).step_by(1) {
        let derating = d as Float / 100.0;
        let upstream_power = upstream_rated_power * derating;
        let effective_thrust = 0.8 * derating.powf(cosine_loss_exponent);
        let downstream_recovery = 1.0 - wake_deficit * (1.0 - effective_thrust);
        let downstream_power = downstream_rated_power * downstream_recovery;
        let total_power = upstream_power + downstream_power;

        if total_power > best_power {
            best_power = total_power;
            best_derating = derating;
        }
    }

    best_derating
}

/// Simple derating model
pub fn simple_derating(available_power: Float, power_setpoint: Float) -> Float {
    available_power.min(power_setpoint)
}

/// Calculate power reduction due to derating
pub fn derating_power_reduction(available_power: Float, power_setpoint: Float) -> Float {
    (available_power - power_setpoint).max(0.0)
}

/// Estimate optimal setpoint based on wind conditions
pub fn estimate_optimal_setpoint(
    wind_speed: Float,
    rated_power: Float,
    rated_wind_speed: Float,
    cut_in_wind_speed: Float,
    derating_factor: Float,
) -> Float {
    if wind_speed < cut_in_wind_speed {
        return 0.0;
    }

    let power_fraction = if wind_speed >= rated_wind_speed {
        1.0
    } else {
        let speed_fraction = (wind_speed - cut_in_wind_speed) / (rated_wind_speed - cut_in_wind_speed);
        speed_fraction.powi(3)
    };

    let setpoint = rated_power * power_fraction * derating_factor;
    setpoint.min(rated_power)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_optimize_derating_reduces_wake() {
        let result = optimize_derating(1000.0, 600.0, 0.3, 0.9);
        assert!(result < 0.9);

        let result2 = optimize_derating(1000.0, 950.0, 0.05, 0.9);
        assert_relative_eq!(result2, 0.9);
    }

    #[test]
    fn test_simple_derating() {
        let result = simple_derating(5000.0, 4000.0);
        assert_relative_eq!(result, 4000.0);

        let result2 = simple_derating(3000.0, 4000.0);
        assert_relative_eq!(result2, 3000.0);
    }

    #[test]
    fn test_derating_power_reduction() {
        let reduction = derating_power_reduction(5000.0, 4000.0);
        assert_relative_eq!(reduction, 1000.0);

        let reduction2 = derating_power_reduction(3000.0, 4000.0);
        assert_relative_eq!(reduction2, 0.0);
    }

    #[test]
    fn test_estimate_optimal_setpoint_below_rated() {
        let setpoint = estimate_optimal_setpoint(7.0, 5_000_000.0, 11.0, 3.0, 1.0);
        assert!(setpoint < 5_000_000.0);
        assert!(setpoint > 0.0);
    }

    #[test]
    fn test_estimate_optimal_setpoint_above_rated() {
        let setpoint = estimate_optimal_setpoint(12.0, 5_000_000.0, 11.0, 3.0, 1.0);
        assert_relative_eq!(setpoint, 5_000_000.0);
    }

    #[test]
    fn test_estimate_optimal_setpoint_below_cut_in() {
        let setpoint = estimate_optimal_setpoint(2.0, 5_000_000.0, 11.0, 3.0, 1.0);
        assert_relative_eq!(setpoint, 0.0);
    }

    #[test]
    fn test_optimize_derating_for_farm_power() {
        let derating = optimize_derating_for_farm_power(
            5_000_000.0, 5_000_000.0, 0.2, 1.0,
        );
        assert!(derating >= 0.5);
        assert!(derating <= 1.0);
    }
}
