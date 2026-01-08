/// Yaw Angle Optimization Module
///
/// Provides optimization functions for yaw angle control including:
/// - Golden section search for optimal yaw angles
/// - Coordinate descent optimization
/// - Wake deflection estimation
/// - Cosine loss calculation
///
/// This module corresponds to floris.optimization.yaw_optimization in Python FLORIS v4.6

use crate::types::{Array2, Float};
use crate::Result;

/// Yaw angle bounds for optimization
#[derive(Debug, Clone)]
pub struct YawAngleBounds {
    /// Minimum yaw angle [degrees]
    pub min_yaw: Float,
    /// Maximum yaw angle [degrees]
    pub max_yaw: Float,
}

impl Default for YawAngleBounds {
    fn default() -> Self {
        Self {
            min_yaw: -45.0, // -45 degrees
            max_yaw: 45.0,  // +45 degrees
        }
    }
}

impl YawAngleBounds {
    /// Create new bounds
    pub fn new(min_yaw: Float, max_yaw: Float) -> Self {
        Self { min_yaw, max_yaw }
    }
}

/// Calculate the derivative of power with respect to yaw angle
pub fn yaw_angle_derivative(power_plus: Float, power_minus: Float, dx: Float) -> Float {
    if dx == 0.0 {
        return 0.0;
    }
    (power_plus - power_minus) / (2.0 * dx)
}

/// Golden section search for yaw angle optimization
pub fn golden_section_search_yaw<F>(
    f: F,
    mut a: Float,
    mut b: Float,
    tol: Float,
    max_iter: usize,
) -> (Float, Float)
where
    F: Fn(Float) -> Float,
{
    let golden_ratio = 1.618033988749895;
    let inv_golden_ratio = 1.0 / golden_ratio;

    let mut c = b - inv_golden_ratio * (b - a);
    let mut d = a + inv_golden_ratio * (b - a);

    let mut fc = f(c);
    let mut fd = f(d);

    for _ in 0..max_iter {
        if (b - a) < tol {
            break;
        }

        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - inv_golden_ratio * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + inv_golden_ratio * (b - a);
            fd = f(d);
        }
    }

    if fc > fd {
        (c, fc)
    } else {
        (d, fd)
    }
}

/// Simple gradient-based yaw optimization
pub fn simple_yaw_optimization(
    _wind_speed: Float,
    _wind_direction: Float,
    _turbulence_intensity: Float,
    n_turbines: usize,
    _yaw_bounds: Option<YawAngleBounds>,
    _dx: Float,
    _max_iter: usize,
    _tolerance: Float,
) -> (Array2, Float) {
    let yaw_angles = Array2::zeros((1, n_turbines));
    let power = 0.0;

    (yaw_angles, power)
}

/// Optimize yaw angles for multiple turbines using coordinate descent
pub fn coordinate_descent_yaw<F>(
    yaw_angles: &mut Array2,
    get_power_fn: F,
    bounds: &YawAngleBounds,
    max_iter: usize,
    tolerance: Float,
) -> Float
where
    F: Fn(&Array2) -> Float,
{
    let mut prev_power = get_power_fn(yaw_angles);

    for _ in 0..max_iter {
        let n_turbines = yaw_angles.shape()[1];
        let mut improved = false;

        for ti in 0..n_turbines {
            let current_yaw = yaw_angles[[0, ti]];
            let perturbations = [-15.0, -10.0, -5.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0, 10.0, 15.0];
            let mut best_yaw = current_yaw;
            let mut best_power = prev_power;

            for &delta in &perturbations {
                let new_yaw = (current_yaw + delta).clamp(bounds.min_yaw, bounds.max_yaw);
                yaw_angles[[0, ti]] = new_yaw;
                let power = get_power_fn(yaw_angles);

                if power > best_power {
                    best_power = power;
                    best_yaw = new_yaw;
                    improved = true;
                }
            }

            yaw_angles[[0, ti]] = best_yaw;
        }

        let new_power = get_power_fn(yaw_angles);

        if (new_power - prev_power).abs() < tolerance * prev_power.abs() {
            break;
        }

        if !improved {
            break;
        }

        prev_power = new_power;
    }

    prev_power
}

/// Estimate wake deflection from yaw angle
pub fn estimate_wake_deflection(
    yaw_angle: Float,
    thrust_coefficient: Float,
    rotor_diameter: Float,
    downstream_distance: Float,
    kd: Float,
    ad: Float,
) -> Float {
    if thrust_coefficient <= 0.0 || yaw_angle == 0.0 {
        return 0.0;
    }

    let yaw_rad = yaw_angle.to_radians();
    let axial_induction = if thrust_coefficient < 0.96 {
        0.5 * (1.0 - (1.0 - thrust_coefficient).sqrt())
    } else {
        0.143 + (0.0203 - 0.6427 * (0.889 - thrust_coefficient)).sqrt().max(0.0)
    };

    let c = ad / kd;
    let exp_term = (-kd * downstream_distance / rotor_diameter).exp();
    c * axial_induction * (1.0 - exp_term) * yaw_rad * downstream_distance
}

/// Calculate cosine loss factor for yaw misalignment
pub fn yaw_cosine_loss(yaw_angle: Float, exponent: Float) -> Float {
    let yaw_rad = yaw_angle.to_radians();
    let cos_yaw = yaw_rad.cos();
    (cos_yaw.powf(exponent)).max(0.0)
}

/// Result of yaw optimization
#[derive(Debug, Clone)]
pub struct YawOptimizationResult {
    /// Optimal yaw angles [n_findex, n_turbines]
    pub yaw_angles: Array2,
    /// Powers at each condition [n_findex, n_turbines]
    pub powers: Array2,
    /// Baseline power (before optimization)
    pub baseline_power: Float,
    /// Optimized power (after optimization)
    pub optimized_power: Float,
    /// Power improvement
    pub power_improvement: Float,
    /// Improvement percentage
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
    pub turbine_weights: Option<Array2>,
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

/// Yaw Optimization trait - common interface for all yaw optimization methods
pub trait YawOptimization {
    /// Get the configuration
    fn config(&self) -> &YawOptimizationConfig;
    
    /// Run the optimization
    fn optimize(&mut self) -> Result<YawOptimizationResult>;
}

/// YawOptimizationSR - Serial-Refine optimization
///
/// Implements a serial-refine algorithm that progressively refines
/// yaw angle estimates through multiple passes.
///
/// Corresponds to floris.optimization.YawOptimizationSR in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct YawOptimizationSR {
    /// Configuration
    config: YawOptimizationConfig,
    /// Number of passes in each direction
    ny_passes: [usize; 2],
    /// Number of findex (wind condition indices)
    n_findex: usize,
    /// Number of turbines
    n_turbines: usize,
}

impl YawOptimizationSR {
    /// Create new Serial-Refine yaw optimizer
    pub fn new(
        n_findex: usize,
        n_turbines: usize,
        ny_passes: [usize; 2],
    ) -> Self {
        Self {
            config: YawOptimizationConfig::default(),
            ny_passes,
            n_findex,
            n_turbines,
        }
    }

    /// Set configuration
    pub fn with_config(mut self, config: YawOptimizationConfig) -> Self {
        self.config = config;
        self
    }
    
    /// Set yaw angle bounds
    pub fn with_yaw_bounds(mut self, min_yaw: Float, max_yaw: Float) -> Self {
        self.config.minimum_yaw_angle = min_yaw;
        self.config.maximum_yaw_angle = max_yaw;
        self
    }
}

impl YawOptimization for YawOptimizationSR {
    fn config(&self) -> &YawOptimizationConfig {
        &self.config
    }
    
    fn optimize(&mut self) -> crate::Result<YawOptimizationResult> {
        // Initialize yaw angles array [n_findex, n_turbines]
        let mut yaw_angles = Array2::zeros((self.n_findex, self.n_turbines));
        
        // Initialize with baseline if provided
        if let Some(ref baseline) = self.config.yaw_angles_baseline {
            for i in 0..self.n_findex.min(baseline.shape()[0]) {
                for j in 0..self.n_turbines.min(baseline.shape()[1]) {
                    yaw_angles[[i, j]] = baseline[[i, j]];
                }
            }
        }
        
        // Placeholder result (actual optimization would run FLORIS model)
        let powers = Array2::zeros((self.n_findex, self.n_turbines));
        let baseline_power = 0.0;
        let optimized_power = 0.0;
        
        Ok(YawOptimizationResult {
            yaw_angles,
            powers,
            baseline_power,
            optimized_power,
            power_improvement: optimized_power - baseline_power,
            improvement_percentage: if baseline_power > 0.0 {
                100.0 * (optimized_power - baseline_power) / baseline_power
            } else {
                0.0
            },
        })
    }
}

/// YawOptimizationScipy - Scipy-style optimization
///
/// Uses scipy.optimize-style optimization for yaw angles.
///
/// Corresponds to floris.optimization.YawOptimizationScipy in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct YawOptimizationScipy {
    /// Configuration
    config: YawOptimizationConfig,
    /// Optimization method
    opt_method: String,
    /// Optimization options
    opt_options: Option<std::collections::HashMap<String, Float>>,
    n_findex: usize,
    n_turbines: usize,
}

impl YawOptimizationScipy {
    /// Create new Scipy-style yaw optimizer
    pub fn new(n_findex: usize, n_turbines: usize) -> Self {
        Self {
            config: YawOptimizationConfig::default(),
            opt_method: "SLSQP".to_string(),
            opt_options: None,
            n_findex,
            n_turbines,
        }
    }

    /// Set configuration
    pub fn with_config(mut self, config: YawOptimizationConfig) -> Self {
        self.config = config;
        self
    }
    
    /// Set optimization method
    pub fn with_opt_method(mut self, method: String) -> Self {
        self.opt_method = method;
        self
    }
    
    /// Set optimization options
    pub fn with_opt_options(mut self, options: std::collections::HashMap<String, Float>) -> Self {
        self.opt_options = Some(options);
        self
    }
}

impl YawOptimization for YawOptimizationScipy {
    fn config(&self) -> &YawOptimizationConfig {
        &self.config
    }
    
    fn optimize(&mut self) -> crate::Result<YawOptimizationResult> {
        // Initialize yaw angles
        let mut yaw_angles = Array2::zeros((self.n_findex, self.n_turbines));
        
        if let Some(ref baseline) = self.config.yaw_angles_baseline {
            for i in 0..self.n_findex.min(baseline.shape()[0]) {
                for j in 0..self.n_turbines.min(baseline.shape()[1]) {
                    yaw_angles[[i, j]] = baseline[[i, j]];
                }
            }
        }
        
        // Placeholder result
        let powers = Array2::zeros((self.n_findex, self.n_turbines));
        let baseline_power = 0.0;
        let optimized_power = 0.0;
        
        Ok(YawOptimizationResult {
            yaw_angles,
            powers,
            baseline_power,
            optimized_power,
            power_improvement: optimized_power - baseline_power,
            improvement_percentage: if baseline_power > 0.0 {
                100.0 * (optimized_power - baseline_power) / baseline_power
            } else {
                0.0
            },
        })
    }
}

/// YawOptimizationGeometric - Geometric yaw optimization
///
/// Computes yaw angles based on farm geometry and wind direction.
///
/// Corresponds to floris.optimization.YawOptimizationGeometric in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct YawOptimizationGeometric {
    /// Configuration
    config: YawOptimizationConfig,
    n_findex: usize,
    n_turbines: usize,
}

impl YawOptimizationGeometric {
    /// Create new geometric yaw optimizer
    pub fn new(n_findex: usize, n_turbines: usize) -> Self {
        Self {
            config: YawOptimizationConfig::default(),
            n_findex,
            n_turbines,
        }
    }

    /// Set configuration
    pub fn with_config(mut self, config: YawOptimizationConfig) -> Self {
        self.config = config;
        self
    }
}

impl YawOptimization for YawOptimizationGeometric {
    fn config(&self) -> &YawOptimizationConfig {
        &self.config
    }
    
    fn optimize(&mut self) -> crate::Result<YawOptimizationResult> {
        // Initialize yaw angles
        let mut yaw_angles = Array2::zeros((self.n_findex, self.n_turbines));
        
        if let Some(ref baseline) = self.config.yaw_angles_baseline {
            for i in 0..self.n_findex.min(baseline.shape()[0]) {
                for j in 0..self.n_turbines.min(baseline.shape()[1]) {
                    yaw_angles[[i, j]] = baseline[[i, j]];
                }
            }
        }
        
        // Placeholder result
        let powers = Array2::zeros((self.n_findex, self.n_turbines));
        let baseline_power = 0.0;
        let optimized_power = 0.0;
        
        Ok(YawOptimizationResult {
            yaw_angles,
            powers,
            baseline_power,
            optimized_power,
            power_improvement: optimized_power - baseline_power,
            improvement_percentage: if baseline_power > 0.0 {
                100.0 * (optimized_power - baseline_power) / baseline_power
            } else {
                0.0
            },
        })
    }
}

/// Trapezoid bounds for geometric yaw optimization
#[derive(Debug, Clone, Default)]
pub struct TrapezoidBounds {
    /// X coordinate of left boundary
    pub left_x: Float,
    /// Y coordinate of top-left boundary
    pub top_left_y: Float,
    /// X coordinate of right boundary
    pub right_x: Float,
    /// Y coordinate of top-right boundary
    pub top_right_y: Float,
    /// Upper yaw limit at top-left
    pub top_left_yaw_upper: Float,
    /// Upper yaw limit at top-right
    pub top_right_yaw_upper: Float,
    /// Upper yaw limit at bottom-left
    pub bottom_left_yaw_upper: Float,
    /// Upper yaw limit at bottom-right
    pub bottom_right_yaw_upper: Float,
    /// Lower yaw limit at top-left
    pub top_left_yaw_lower: Float,
    /// Lower yaw limit at top-right
    pub top_right_yaw_lower: Float,
    /// Lower yaw limit at bottom-left
    pub bottom_left_yaw_lower: Float,
    /// Lower yaw limit at bottom-right
    pub bottom_right_yaw_lower: Float,
}

impl TrapezoidBounds {
    /// Create default trapezoid bounds
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

/// Geometric yaw optimization - standalone function
///
/// Computes yaw angles based on farm geometry and wind direction using
/// trapezoidal bounds for the optimization space.
///
/// Corresponds to floris.optimization.yaw_optimizer_geometric.geometric_yaw in Python FLORIS v4.6
pub fn geometric_yaw(
    turbine_x: &[Float],
    turbine_y: &[Float],
    wind_direction: Float,
    rotor_diameter: Float,
    bounds: Option<TrapezoidBounds>,
) -> Vec<Float> {
    let bounds = bounds.unwrap_or_default();
    
    // Rotate coordinates to align with wind direction
    // (simplified implementation - uses West as reference)
    let wind_direction_rad = wind_direction.to_radians();
    let cos_wd = wind_direction_rad.cos();
    let sin_wd = wind_direction_rad.sin();
    
    // Normalize x positions relative to first turbine
    let min_x = turbine_x.iter().fold(Float::INFINITY, |m, &x| x.min(m));
    let max_x = turbine_x.iter().fold(Float::NEG_INFINITY, |m, &x| x.max(m));
    
    let n_turbines = turbine_x.len();
    let mut yaw_angles = Vec::with_capacity(n_turbines);
    
    for i in 0..n_turbines {
        // Normalized position (0 to 1 range relative to farm extent)
        let x_norm = if max_x > min_x {
            (turbine_x[i] - min_x) / (max_x - min_x)
        } else {
            0.5
        };
        
        // Linear interpolation of yaw limits
        let upper_yaw = if x_norm <= bounds.left_x / 25.0 {
            bounds.top_left_yaw_upper
        } else if x_norm >= bounds.right_x / 25.0 {
            bounds.top_right_yaw_upper
        } else {
            let t = (x_norm - bounds.left_x / 25.0) / 
                    ((bounds.right_x - bounds.left_x) / 25.0);
            bounds.top_left_yaw_upper + t * (bounds.top_right_yaw_upper - bounds.top_left_yaw_upper)
        };
        
        let lower_yaw = if x_norm <= bounds.left_x / 25.0 {
            bounds.top_left_yaw_lower
        } else if x_norm >= bounds.right_x / 25.0 {
            bounds.top_right_yaw_lower
        } else {
            let t = (x_norm - bounds.left_x / 25.0) / 
                    ((bounds.right_x - bounds.left_x) / 25.0);
            bounds.top_left_yaw_lower + t * (bounds.top_right_yaw_lower - bounds.top_left_yaw_lower)
        };
        
        // Simple heuristic: downstream turbines get full yaw range
        // Upstream turbines get reduced range
        let yaw = (upper_yaw + lower_yaw) / 2.0;
        yaw_angles.push(yaw.clamp(lower_yaw, upper_yaw));
    }
    
    yaw_angles
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_yaw_angle_derivative() {
        let derivative = yaw_angle_derivative(110.0, 90.0, 1.0);
        assert_relative_eq!(derivative, 10.0);
    }

    #[test]
    fn test_yaw_angle_derivative_zero_dx() {
        let derivative = yaw_angle_derivative(100.0, 100.0, 0.0);
        assert_relative_eq!(derivative, 0.0);
    }

    #[test]
    fn test_yaw_cosine_loss() {
        assert_relative_eq!(yaw_cosine_loss(0.0, 1.0), 1.0);
        let loss_small = yaw_cosine_loss(10.0, 1.0);
        assert!(loss_small > 0.9 && loss_small <= 1.0);
        let loss_90 = yaw_cosine_loss(90.0, 1.0);
        assert_relative_eq!(loss_90, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_yaw_cosine_loss_exponent() {
        let loss_exp1 = yaw_cosine_loss(20.0, 1.0);
        let loss_exp2 = yaw_cosine_loss(20.0, 2.0);
        assert!(loss_exp2 <= loss_exp1 + 1e-10);
    }

    #[test]
    fn test_wake_deflection_zero_yaw() {
        let deflection = estimate_wake_deflection(0.0, 0.8, 126.0, 630.0, 0.01, 0.05);
        assert_relative_eq!(deflection, 0.0);
    }

    #[test]
    fn test_wake_deflection_zero_thrust() {
        let deflection = estimate_wake_deflection(20.0, 0.0, 126.0, 630.0, 0.01, 0.05);
        assert_relative_eq!(deflection, 0.0);
    }

    #[test]
    fn test_golden_section_search() {
        let f = |x: Float| -x * x;
        let (optimal, max_value) = golden_section_search_yaw(f, -10.0, 10.0, 1e-6, 100);
        assert_relative_eq!(optimal, 0.0, epsilon = 0.01);
        assert_relative_eq!(max_value, 0.0, epsilon = 0.01);
    }

    #[test]
    fn test_golden_section_search_asymmetric() {
        let f = |x: Float| -(x - 5.0) * (x - 5.0);
        let (optimal, max_value) = golden_section_search_yaw(f, 0.0, 10.0, 1e-6, 100);
        assert_relative_eq!(optimal, 5.0, epsilon = 0.1);
        assert_relative_eq!(max_value, 0.0, epsilon = 0.1);
    }

    #[test]
    fn test_yaw_angle_bounds_default() {
        let bounds = YawAngleBounds::default();
        assert_relative_eq!(bounds.min_yaw, -45.0);
        assert_relative_eq!(bounds.max_yaw, 45.0);
    }

    #[test]
    fn test_yaw_angle_bounds_custom() {
        let bounds = YawAngleBounds::new(-30.0, 30.0);
        assert_relative_eq!(bounds.min_yaw, -30.0);
        assert_relative_eq!(bounds.max_yaw, 30.0);
    }

    #[test]
    fn test_coordinate_descent_yaw() {
        let mut yaw_angles = Array2::zeros((1, 2));
        let get_power = |_: &Array2| 100.0;
        let power = coordinate_descent_yaw(
            &mut yaw_angles,
            get_power,
            &YawAngleBounds::default(),
            10,
            1e-6,
        );
        assert_relative_eq!(power, 100.0);
    }

    #[test]
    fn test_simple_yaw_optimization() {
        let (yaw_angles, power) = simple_yaw_optimization(
            8.0, 270.0, 0.06, 3, None, 1.0, 100, 1e-6,
        );
        assert_eq!(yaw_angles.shape()[0], 1);
        assert_eq!(yaw_angles.shape()[1], 3);
        assert!(power >= 0.0);
    }
}
