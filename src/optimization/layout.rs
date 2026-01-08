/// Layout Optimization Module
///
/// Provides layout optimization functions for wind farm layout optimization.
///
/// This module corresponds to floris.optimization.layout_optimization in Python FLORIS v4.6

use crate::types::{Array1, Array2, Float};
use crate::Result;

/// Boundary type for layout optimization
#[derive(Debug, Clone)]
pub enum Boundary {
    /// Polygon boundary defined by x, y coordinates
    Polygon(Vec<(Float, Float)>),
    /// Circle boundary defined by center and radius
    Circle { center_x: Float, center_y: Float, radius: Float },
    /// Rectangle boundary
    Rectangle { min_x: Float, max_x: Float, min_y: Float, max_y: Float },
}

/// Configuration for layout optimization
#[derive(Debug, Clone)]
pub struct LayoutOptimizationConfig {
    /// Boundary constraints
    pub boundaries: Boundary,
    /// Minimum distance between turbines [meters]
    pub min_dist: Option<Float>,
    /// Optimization solver
    pub solver: String,
    /// Enable geometric yaw during layout optimization
    pub enable_geometric_yaw: bool,
    /// Optimize value (AVP) instead of power (AEP)
    pub use_value: bool,
}

impl Default for LayoutOptimizationConfig {
    fn default() -> Self {
        Self {
            boundaries: Boundary::Rectangle {
                min_x: 0.0,
                max_x: 5000.0,
                min_y: 0.0,
                max_y: 5000.0,
            },
            min_dist: Some(500.0), // 5 rotor diameters for 100m turbines
            solver: "SLSQP".to_string(),
            enable_geometric_yaw: false,
            use_value: false,
        }
    }
}

/// Result of layout optimization
#[derive(Debug, Clone)]
pub struct LayoutOptimizationResult {
    /// Optimized turbine x coordinates
    pub x: Array1,
    /// Optimized turbine y coordinates
    pub y: Array1,
    /// AEP or AVP at optimized layout
    pub value: Float,
    /// Number of iterations
    pub iterations: usize,
}

impl Default for LayoutOptimizationResult {
    fn default() -> Self {
        Self {
            x: Array1::from_vec(vec![]),
            y: Array1::from_vec(vec![]),
            value: 0.0,
            iterations: 0,
        }
    }
}

/// LayoutOptimizationScipy - Layout optimization using scipy-style optimization
///
/// Corresponds to floris.optimization.LayoutOptimizationScipy in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct LayoutOptimizationScipy {
    /// Configuration
    config: LayoutOptimizationConfig,
    /// Number of turbines
    n_turbines: usize,
}

impl LayoutOptimizationScipy {
    /// Create new layout optimizer
    pub fn new(n_turbines: usize, boundaries: Boundary) -> Self {
        Self {
            config: LayoutOptimizationConfig {
                boundaries,
                ..Default::default()
            },
            n_turbines,
        }
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

    /// Set solver
    pub fn with_solver(mut self, solver: String) -> Self {
        self.config.solver = solver;
        self
    }

    /// Enable geometric yaw
    pub fn with_geometric_yaw(mut self, enabled: bool) -> Self {
        self.config.enable_geometric_yaw = enabled;
        self
    }

    /// Optimize layout
    pub fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        // Placeholder implementation
        let x = Array1::from_vec(vec![0.0; self.n_turbines]);
        let y = Array1::from_vec(vec![0.0; self.n_turbines]);

        Ok(LayoutOptimizationResult {
            x,
            y,
            value: 0.0,
            iterations: 0,
        })
    }
}

/// LayoutOptimizationRandomSearch - Random search layout optimization
///
/// Corresponds to floris.optimization.LayoutOptimizationRandomSearch in Python FLORIS v4.6
#[derive(Debug, Clone)]
pub struct LayoutOptimizationRandomSearch {
    /// Configuration
    config: LayoutOptimizationConfig,
    /// Number of turbines
    n_turbines: usize,
    /// Number of random samples
    n_samples: usize,
}

impl LayoutOptimizationRandomSearch {
    /// Create new random search layout optimizer
    pub fn new(n_turbines: usize, boundaries: Boundary, n_samples: usize) -> Self {
        Self {
            config: LayoutOptimizationConfig {
                boundaries,
                ..Default::default()
            },
            n_turbines,
            n_samples,
        }
    }

    /// Optimize layout using random search
    pub fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
        // Placeholder implementation
        let x = Array1::from_vec(vec![0.0; self.n_turbines]);
        let y = Array1::from_vec(vec![0.0; self.n_turbines]);

        Ok(LayoutOptimizationResult {
            x,
            y,
            value: 0.0,
            iterations: self.n_samples,
        })
    }
}

/// Check if a point is inside the boundary
pub fn is_point_in_boundary(x: Float, y: Float, boundary: &Boundary) -> bool {
    match boundary {
        Boundary::Polygon(points) => {
            // Ray casting algorithm
            let mut inside = false;
            let mut j = points.len() - 1;
            for (i, &(xi, yi)) in points.iter().enumerate() {
                let (xj, yj) = points[j];
                if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
                    inside = !inside;
                }
                j = i;
            }
            inside
        }
        Boundary::Circle { center_x, center_y, radius } => {
            let dx = x - center_x;
            let dy = y - center_y;
            dx * dx + dy * dy <= radius * radius
        }
        Boundary::Rectangle { min_x, max_x, min_y, max_y } => {
            x >= *min_x && x <= *max_x && y >= *min_y && y <= *max_y
        }
    }
}

/// Generate grid points within boundary
pub fn generate_grid_points(
    boundary: &Boundary,
    n_points: usize,
) -> Vec<(Float, Float)> {
    let mut points = Vec::with_capacity(n_points);
    
    match boundary {
        Boundary::Rectangle { min_x, max_x, min_y, max_y } => {
            let n_side = (n_points as Float).sqrt().ceil() as usize;
            let dx = (max_x - min_x) / (n_side as Float + 1.0);
            let dy = (max_y - min_y) / (n_side as Float + 1.0);
            
            for i in 1..=n_side {
                for j in 1..=n_side {
                    let x = min_x + i as Float * dx;
                    let y = min_y + j as Float * dy;
                    if points.len() < n_points {
                        points.push((x, y));
                    }
                }
            }
        }
        _ => {
            // For non-rectangular boundaries, generate random points
            let mut rng = fastrand::Rng::new();
            for _ in 0..n_points {
                match boundary {
                    Boundary::Circle { center_x, center_y, radius } => {
                        let angle = rng.f32() * 2.0 * std::f32::consts::PI as Float;
                        let r = radius * rng.f32().sqrt();
                        let x = center_x + r * angle.cos();
                        let y = center_y + r * angle.sin();
                        points.push((x as Float, y as Float));
                    }
                    _ => {
                        // Fallback to (0, 0)
                        points.push((0.0, 0.0));
                    }
                }
            }
        }
    }
    
    points
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_is_point_in_rectangle() {
        let boundary = Boundary::Rectangle {
            min_x: 0.0,
            max_x: 100.0,
            min_y: 0.0,
            max_y: 100.0,
        };
        
        assert!(is_point_in_boundary(50.0, 50.0, &boundary));
        assert!(is_point_in_boundary(0.0, 0.0, &boundary));
        assert!(!is_point_in_boundary(150.0, 50.0, &boundary));
        assert!(!is_point_in_boundary(50.0, 150.0, &boundary));
    }

    #[test]
    fn test_is_point_in_circle() {
        let boundary = Boundary::Circle {
            center_x: 0.0,
            center_y: 0.0,
            radius: 10.0,
        };
        
        assert!(is_point_in_boundary(0.0, 0.0, &boundary));
        assert!(is_point_in_boundary(5.0, 5.0, &boundary));
        assert!(!is_point_in_boundary(10.0, 0.0, &boundary));
        assert!(!is_point_in_boundary(15.0, 0.0, &boundary));
    }

    #[test]
    fn test_layout_optimization_config_default() {
        let config = LayoutOptimizationConfig::default();
        assert!(config.min_dist.is_some());
        assert_eq!(config.solver, "SLSQP");
        assert!(!config.enable_geometric_yaw);
        assert!(!config.use_value);
    }

    #[test]
    fn test_layout_optimization_result_default() {
        let result = LayoutOptimizationResult::default();
        assert!(result.x.is_empty());
        assert!(result.y.is_empty());
        assert_eq!(result.value, 0.0);
        assert_eq!(result.iterations, 0);
    }
}
