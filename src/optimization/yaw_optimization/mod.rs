//! Yaw optimization module for FLORUS
//!
//! Provides optimization algorithms for wind turbine yaw control including:
//! - Base optimizer with common functionality (yaw_optimization_base)
//! - Geometric yaw optimizer for fast layout estimates (yaw_optimizer_geometric)
//! - SciPy-based gradient optimizer for precise control (yaw_optimizer_scipy)
//! - Serial Refine optimizer for large farm optimization (yaw_optimizer_sr)
//! - Helper functions for wake analysis (yaw_optimization_tools)
//!
//! This module corresponds to floris/optimization/yaw_optimization/ in Python FLORIS v4.6
//!
//! # Module Structure
//!
//! ## Base Module
//! - **yaw_optimization_base**: Contains the `YawOptimization` trait/base class
//!   with common functionality for constraint handling, variable normalization,
//!   downstream turbine exclusion, and power calculation.
//!
//! ## Optimizer Implementations
//!
//! ### 1. Geometric Optimizer (yaw_optimizer_geometric)
//! Fast, geometry-based yaw optimization suitable for:
//! - Layout optimization quick estimates
//! - Large wind farms
//! - Cases where approximate solutions are acceptable
//!
//! Algorithm: Based on trapezoid wake model and turbine positioning.
//!
//! ### 2. SciPy Optimizer (yaw_optimizer_scipy)
//! Gradient-based optimization using coordinate descent.
//! Suitable for:
//! - Precise yaw control
//! - Small to medium wind farms
//!
//! ### 3. Serial Refine Optimizer (yaw_optimizer_sr)
//! Custom serial refinement algorithm that optimizes turbines progressively.
//! Suitable for:
//! - Large wind farms
//! - AEP optimization
//!
//! ## Helper Tools
//!
//! **yaw_optimization_tools**: Provides utility functions:
//! - `derive_downstream_turbines()`: Identify turbines unaffected by others' wakes
//! - `yaw_cosine_loss()`: Calculate cosine loss from yaw angles
//! - `estimate_wake_deflection_angle()`: Estimate wake steering effects
//!
//! # Usage Examples
//!
//! ```rust,ignore
//! use florus::optimization::yaw_optimization::YawOptimization;
//!
//! // Using the trait directly
//! let result = yaw_optimization::simple_yaw_optimization(fmodel, None)?;
//! ```
//!
//! ```rust,ignore
//! use florus::optimization::yaw_optimization::YawOptimizationGeometric;
//!
//! // Geometric optimizer for fast estimates
//! let result = YawOptimizationGeometric::optimize(fmodel, None)?;
//! ```
//!
//! # Coordinate System
//!
//! All yaw angles are in degrees relative to the wind direction.
//! Positive values rotate the turbine counterclockwise (left),
//! negative values rotate clockwise (right) when facing the wind.

pub mod yaw_optimization_base;
pub mod yaw_optimization_tools;
pub mod yaw_optimizer_geometric;
pub mod yaw_optimizer_scipy;
pub mod yaw_optimizer_sr;

// Re-export main types and functions from submodules
pub use yaw_optimization_base::{
    control_problem, geometric_yaw, simple_yaw_optimization, variable_handling, TrapezoidBounds,
    YawAngleBounds, YawOptimization, YawOptimizationConfig, YawOptimizationResult,
};

pub use yaw_optimization_tools::{
    calculate_turbine_weights, coordinate_descent_yaw, derive_downstream_turbines,
    estimate_wake_deflection_angle, golden_section_search_yaw, is_turbine_in_wake,
    optimize_yaw_single_findex, yaw_cosine_loss, yaw_cosine_loss_derivative,
};

pub use yaw_optimizer_geometric::YawOptimizationGeometric;
pub use yaw_optimizer_scipy::{SciPyOptimizationConfig, YawOptimizationScipy};
pub use yaw_optimizer_sr::{SerialRefineConfig, YawOptimizationSR};
