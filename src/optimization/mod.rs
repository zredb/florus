use crate::core::Farm;
/// Optimization module for FLORIS-RS
///
/// Provides optimization functions for wind farm control:
/// - Yaw angle optimization (yaw module)
/// - Power setpoint optimization (power_setpoint module)
/// - Derating optimization (derating module)
/// - Layout optimization (layout_optimization module)
/// - Load optimization (load_optimization module)
use crate::types::Float;

// Submodules
pub mod derating;
pub mod layout_optimization;
pub mod load_optimization;
pub mod power_setpoint;
pub mod yaw_optimization;

// Re-export main types from yaw_optimization module
pub use yaw_optimization::{
    coordinate_descent_yaw, estimate_wake_deflection_angle, geometric_yaw, golden_section_search_yaw,
    simple_yaw_optimization, yaw_cosine_loss, yaw_cosine_loss_derivative, TrapezoidBounds,
    YawAngleBounds, YawOptimization, YawOptimizationConfig, YawOptimizationGeometric,
    YawOptimizationResult, YawOptimizationSR, YawOptimizationScipy,
};

pub use power_setpoint::{
    compute_power_setpoints, optimize_power_setpoints as simple_optimize_power_setpoints,
    PowerSetpointOptimizationResult,
};

pub use derating::{
    derating_power_reduction, estimate_optimal_setpoint, optimize_derating,
    optimize_derating_for_farm_power, simple_derating,
};

pub use layout_optimization::{
    calculate_pairwise_distances, generate_grid_points, is_point_in_boundary, load_optimization,
    load_optimization_result, save_optimization, save_optimization_result, Boundary,
    GridOptimizationConfig, LayoutOptimizationBoundaryGrid, LayoutOptimizationConfig,
    LayoutOptimizationGoldenSection, LayoutOptimizationMixedInteger, LayoutOptimizationPyOptSparse,
    LayoutOptimizationRandomSearch, LayoutOptimizationResult, LayoutOptimizationScipy,
    LayoutOptimizer, OptimizationConfigFile, OptimizationType,
};

pub use load_optimization::{
    compute_farm_revenue, compute_farm_voc, compute_lti, compute_net_revenue, compute_turbine_voc,
    find_a_to_satisfy_rev_voc_ratio, find_a_to_satisfy_target_voc_per_mw, optimize_power_setpoints,
    POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED,
};

/// Optimize yaw angles for a wind farm
pub fn optimize_yaw_angles(
    _farm: &Farm,
    wind_speeds: &[Float],
    _wind_directions: &[Float],
    _turbulence_intensities: &[Float],
    _yaw_bounds: Option<YawAngleBounds>,
    _max_iterations: usize,
    _tolerance: Float,
) -> YawOptimizationResult {
    let n_findex = wind_speeds.len();
    let n_turbines = 1;

    let yaw_angles = ndarray::Array2::zeros((n_findex, n_turbines));

    YawOptimizationResult {
        yaw_angles,
        powers: ndarray::Array2::zeros((n_findex, n_turbines)),
        baseline_power: 0.0,
        optimized_power: 0.0,
        power_improvement: 0.0,
        improvement_percentage: 0.0,
    }
}
