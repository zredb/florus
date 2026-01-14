//! Layout Optimization Base Module
//!
//! Provides shared types, boundary definitions, and the base trait for layout optimization.
//! This module corresponds to layout_optimization_base.py in Python FLORIS v4.6.

use crate::types::{Array1, Float};
use crate::Result;
use ndarray::Array;
use rand::Rng;
use serde::{Deserialize, Serialize};

/// Boundary type for layout optimization
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Boundary {
    /// Polygon boundary defined by x, y coordinates
    Polygon(Vec<(Float, Float)>),
    /// Rectangle boundary
    Rectangle { min_x: Float, max_x: Float, min_y: Float, max_y: Float },
}

impl Default for Boundary {
    fn default() -> Self {
        Boundary::Rectangle {
            min_x: 0.0,
            max_x: 5000.0,
            min_y: 0.0,
            max_y: 5000.0,
        }
    }
}

impl Boundary {
    /// Get the bounds of the boundary
    pub fn bounds(&self) -> (Float, Float, Float, Float) {
        match self {
            Boundary::Polygon(points) => {
                let min_x = points.iter().map(|&(x, _)| x).fold(Float::INFINITY, Float::min);
                let max_x = points.iter().map(|&(x, _)| x).fold(Float::NEG_INFINITY, Float::max);
                let min_y = points.iter().map(|&(_, y)| y).fold(Float::INFINITY, Float::min);
                let max_y = points.iter().map(|&(_, y)| y).fold(Float::NEG_INFINITY, Float::max);
                (min_x, max_x, min_y, max_y)
            }
            Boundary::Rectangle { min_x, max_x, min_y, max_y } => {
                (*min_x, *max_x, *min_y, *max_y)
            }
        }
    }

    /// Check if a point is inside the boundary
    pub fn contains(&self, x: Float, y: Float) -> bool {
        match self {
            Boundary::Polygon(points) => {
                // Ray casting algorithm
                let mut inside = false;
                let mut j = points.len() - 1;
                for (i, &(xi, yi)) in points.iter().enumerate() {
                    let (xj, yj) = points[j];
                    if ((yi > y) != (yj > y)) && 
                       (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
                        inside = !inside;
                    }
                    j = i;
                }
                inside
            }
            Boundary::Rectangle { min_x, max_x, min_y, max_y } => {
                x >= *min_x && x <= *max_x && y >= *min_y && y <= *max_y
            }
        }
    }

    /// Get the polygon boundary as a list of points (converts Rectangle to Polygon)
    pub fn to_polygon(&self) -> Vec<(Float, Float)> {
        match self {
            Boundary::Polygon(points) => points.clone(),
            Boundary::Rectangle { min_x, max_x, min_y, max_y } => {
                vec![
                    (*min_x, *min_y),
                    (*max_x, *min_y),
                    (*max_x, *max_y),
                    (*min_x, *max_y),
                ]
            }
        }
    }
}

/// Configuration for layout optimization
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    /// Maximum iterations
    pub max_iterations: usize,
    /// Optimization tolerance
    pub tolerance: Float,
}

impl Default for LayoutOptimizationConfig {
    fn default() -> Self {
        Self {
            boundaries: Boundary::default(),
            min_dist: Some(500.0),
            solver: "SLSQP".to_string(),
            enable_geometric_yaw: false,
            use_value: false,
            max_iterations: 100,
            tolerance: 1e-9,
        }
    }
}

/// Result of layout optimization
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayoutOptimizationResult {
    /// Optimized turbine x coordinates
    pub x: Array1,
    /// Optimized turbine y coordinates
    pub y: Array1,
    /// AEP or AVP at optimized layout
    pub value: Float,
    /// Number of iterations
    pub iterations: usize,
    /// Improvement from initial layout [%]
    pub improvement_pct: Float,
}

impl Default for LayoutOptimizationResult {
    fn default() -> Self {
        Self {
            x: Array1::from_vec(vec![]),
            y: Array1::from_vec(vec![]),
            value: 0.0,
            iterations: 0,
            improvement_pct: 0.0,
        }
    }
}

/// Base trait for layout optimization
pub trait LayoutOptimizer {
    /// Get the number of turbines
    fn n_turbines(&self) -> usize;
    
    /// Get the minimum distance between turbines
    fn min_dist(&self) -> Float;
    
    /// Get the boundaries
    fn boundaries(&self) -> &Boundary;
    
    /// Calculate the objective function value (AEP or AVP)
    fn calculate_objective(&self, x: &Array1, y: &Array1) -> Float;
    
    /// Normalize coordinates to [0, 1] range
    fn normalize(&self, val: Float, x1: Float, x2: Float) -> Float {
        (val - x1) / (x2 - x1)
    }
    
    /// Denormalize coordinates from [0, 1] range
    fn denormalize(&self, val: Float, x1: Float, x2: Float) -> Float {
        val * (x2 - x1) + x1
    }
    
    /// Calculate minimum distance between any two turbines
    fn calculate_min_distance(&self, x: &Array1, y: &Array1) -> Float {
        let mut min_dist = Float::INFINITY;
        for i in 0..x.len() {
            for j in (i + 1)..x.len() {
                let dx = x[i] - x[j];
                let dy = y[i] - y[j];
                let dist = (dx * dx + dy * dy).sqrt();
                if dist < min_dist {
                    min_dist = dist;
                }
            }
        }
        min_dist
    }
    
    /// Check if all turbines are within boundaries
    fn check_boundary_constraints(&self, x: &Array1, y: &Array1) -> bool {
        for i in 0..x.len() {
            if !self.boundaries().contains(x[i], y[i]) {
                return false;
            }
        }
        true
    }
    
    /// Calculate space constraint violation (returns 0 if satisfied)
    fn space_constraint(&self, x: &Array1, y: &Array1, rho: Float) -> Float {
        let min_dist = self.min_dist();
        let mut distances = Vec::new();
        
        for i in 0..x.len() {
            for j in (i + 1)..x.len() {
                let dx = x[i] - x[j];
                let dy = y[i] - y[j];
                let dist = (dx * dx + dy * dy).sqrt();
                distances.push(dist);
            }
        }
        
        if distances.is_empty() {
            return 0.0;
        }
        
        // Calculate constraint using KS aggregation
        let g: Vec<Float> = distances.iter().map(|&d| 1.0 - d / min_dist).collect();
        let g_max = g.iter().fold(Float::NEG_INFINITY, |a, &b| a.max(b));
        
        let sum_exp: Float = g.iter()
            .map(|&gi| ((gi - g_max) * rho).exp())
            .sum();
        
        // Constraint is satisfied when KS_constraint <= 0
        let ks_constraint = g_max + (1.0 / rho) * sum_exp.ln();
        
        -ks_constraint // Return negative for minimization
    }
    
    /// Run the optimization
    fn optimize(&mut self) -> Result<LayoutOptimizationResult>;
}

/// Check if a point is inside the boundary
pub fn is_point_in_boundary(x: Float, y: Float, boundary: &Boundary) -> bool {
    boundary.contains(x, y)
}

/// Generate grid points within boundary
pub fn generate_grid_points(
    boundary: &Boundary,
    n_points: usize,
) -> Vec<(Float, Float)> {
    let mut points = Vec::with_capacity(n_points);

    let (min_x, max_x, min_y, max_y) = boundary.bounds();
    let n_side = (n_points as Float).sqrt().ceil() as usize;
    let dx = (max_x - min_x) / (n_side as Float + 1.0);
    let dy = (max_y - min_y) / (n_side as Float + 1.0);

    for i in 1..=n_side {
        for j in 1..=n_side {
            let x = min_x + i as Float * dx;
            let y = min_y + j as Float * dy;
            if points.len() < n_points && boundary.contains(x, y) {
                points.push((x, y));
            }
        }
    }

    // Fill remaining with random points if needed
    let mut rng = rand::thread_rng();
    while points.len() < n_points {
        let x = min_x + (max_x - min_x) * rng.gen::<Float>();
        let y = min_y + (max_y - min_y) * rng.gen::<Float>();
        if boundary.contains(x, y) {
            points.push((x, y));
        }
    }

    points
}

/// Calculate pairwise distances between all turbines
pub fn calculate_pairwise_distances(x: &Array1, y: &Array1) -> Array<Float, ndarray::Dim<[usize; 2]>> {
    let n = x.len();
    let mut distances = Array::zeros((n, n));

    for i in 0..n {
        for j in 0..n {
            let dx = x[i] - x[j];
            let dy = y[i] - y[j];
            distances[[i, j]] = (dx * dx + dy * dy).sqrt();
        }
    }

    distances
}

/// Optimization type enumeration for loading
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "optimizer_type")]
pub enum OptimizationType {
    /// Scipy coordinate descent optimizer
    Scipy {
        /// Optimization configuration
        config: LayoutOptimizationConfig,
    },
    /// Random search optimizer
    RandomSearch {
        /// Number of individuals
        n_individuals: usize,
        /// Number of generations
        n_generations: usize,
        /// Configuration
        config: LayoutOptimizationConfig,
    },
    /// PyOptSparse gradient-based optimizer
    PyOptSparse {
        /// Configuration
        config: LayoutOptimizationConfig,
    },
    /// Grid-based optimizer
    Grid {
        /// Grid resolution
        grid_resolution: usize,
        /// Configuration
        config: LayoutOptimizationConfig,
    },
    /// Mixed integer optimizer (grid + continuous)
    MixedInteger {
        /// Grid resolution
        grid_resolution: usize,
        /// Configuration
        config: LayoutOptimizationConfig,
    },
}

/// Top-level optimization configuration file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationConfigFile {
    /// Optimization type and parameters
    #[serde(flatten)]
    pub optimization: OptimizationType,
}

/// Load optimization configuration from a YAML file
///
/// This function reads a YAML file containing optimization configuration
/// and returns the appropriate `LayoutOptimizationConfig`.
///
/// # Arguments
///
/// * `path` - Path to the YAML configuration file
///
/// # Returns
///
/// Returns `Ok(LayoutOptimizationConfig)` if successful, or an error otherwise.
///
/// # Example
///
/// ```ignore
/// use florus::optimization::load_optimization;
///
/// let config = load_optimization("optimization_config.yaml")?;
/// ```
pub fn load_optimization(path: impl AsRef<std::path::Path>) -> Result<LayoutOptimizationConfig> {
    use crate::utilities::load_yaml;

    let config_value = load_yaml(path)?;

    // Try to deserialize as OptimizationConfigFile first
    if let Ok(config_file) = serde_yaml::from_value::<OptimizationConfigFile>(config_value.clone()) {
        match config_file.optimization {
            OptimizationType::Scipy { config, .. } => return Ok(config),
            OptimizationType::RandomSearch { config, .. } => return Ok(config),
            OptimizationType::PyOptSparse { config, .. } => return Ok(config),
            OptimizationType::Grid { config, .. } => return Ok(config),
            OptimizationType::MixedInteger { config, .. } => return Ok(config),
        }
    }

    // Fallback: try to deserialize directly as LayoutOptimizationConfig
    let config = serde_yaml::from_value::<LayoutOptimizationConfig>(config_value)?;
    Ok(config)
}

/// Save optimization configuration to a YAML file
///
/// # Arguments
///
/// * `config` - The optimization configuration to save
/// * `path` - Path to save the YAML file
pub fn save_optimization(config: &LayoutOptimizationConfig, path: impl AsRef<std::path::Path>) -> Result<()> {
    let yaml = serde_yaml::to_string(config)?;
    std::fs::write(path, yaml)?;
    Ok(())
}

/// Load optimization result from a YAML file
///
/// # Arguments
///
/// * `path` - Path to the YAML result file
///
/// # Returns
///
/// Returns `Ok(LayoutOptimizationResult)` if successful, or an error otherwise.
pub fn load_optimization_result(path: impl AsRef<std::path::Path>) -> Result<LayoutOptimizationResult> {
    use crate::utilities::load_yaml;

    let config_value = load_yaml(path)?;
    let result = serde_yaml::from_value::<LayoutOptimizationResult>(config_value)?;
    Ok(result)
}

/// Save optimization result to a YAML file
///
/// # Arguments
///
/// * `result` - The optimization result to save
/// * `path` - Path to save the YAML file
pub fn save_optimization_result(result: &LayoutOptimizationResult, path: impl AsRef<std::path::Path>) -> Result<()> {
    let yaml = serde_yaml::to_string(result)?;
    std::fs::write(path, yaml)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_boundary_bounds() {
        let boundary = Boundary::Rectangle {
            min_x: 0.0,
            max_x: 100.0,
            min_y: 0.0,
            max_y: 100.0,
        };
        
        let (xmin, xmax, ymin, ymax) = boundary.bounds();
        assert_eq!(xmin, 0.0);
        assert_eq!(xmax, 100.0);
        assert_eq!(ymin, 0.0);
        assert_eq!(ymax, 100.0);
    }

    #[test]
    fn test_boundary_contains() {
        let boundary = Boundary::Rectangle {
            min_x: 0.0,
            max_x: 100.0,
            min_y: 0.0,
            max_y: 100.0,
        };
        
        assert!(boundary.contains(50.0, 50.0));
        assert!(boundary.contains(0.0, 0.0));
        assert!(!boundary.contains(150.0, 50.0));
        assert!(!boundary.contains(50.0, 150.0));
    }

    #[test]
    fn test_polygon_boundary() {
        let boundary = Boundary::Polygon(vec![
            (0.0, 0.0),
            (100.0, 0.0),
            (100.0, 100.0),
            (0.0, 100.0),
        ]);
        
        assert!(boundary.contains(50.0, 50.0));
        assert!(!boundary.contains(150.0, 50.0));
    }

    #[test]
    fn test_min_distance_calculation() {
        let x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
        let y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        
        let min_dist = 500.0; // Expected minimum distance
        
        struct TestOptimizer;
        impl LayoutOptimizer for TestOptimizer {
            fn n_turbines(&self) -> usize { 3 }
            fn min_dist(&self) -> Float { 500.0 }
            fn boundaries(&self) -> &Boundary {
                static BOUNDARY: Boundary = Boundary::Rectangle {
                    min_x: -1000.0, max_x: 2000.0,
                    min_y: -1000.0, max_y: 2000.0,
                };
                &BOUNDARY
            }
            fn calculate_objective(&self, _: &Array1, _: &Array1) -> Float { 0.0 }
            fn optimize(&mut self) -> Result<LayoutOptimizationResult> {
                Ok(LayoutOptimizationResult::default())
            }
        }
        
        let optimizer = TestOptimizer;
        assert_relative_eq!(optimizer.calculate_min_distance(&x, &y), min_dist);
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
