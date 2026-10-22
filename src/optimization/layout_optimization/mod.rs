//! Layout optimization module for wind farm turbine positioning
//!
//! Provides optimization algorithms for wind farm layout including:
//! - Random search optimization (genetic algorithm-style)
//! - Scipy-based coordinate descent optimization
//! - PyOptSparse-style gradient-based optimization (argmin)
//! - Grid-based exhaustive search optimization
//! - Mixed integer optimization (grid + continuous refinement)
//!
//! The scipy-based optimization uses the SLSQP solver of
//! [`slsqp`], which is a port of scipy's solver.
//!
//! This module corresponds to floris/optimization/layout_optimization/ in Python FLORIS v4.6

// Module declarations
pub mod layout_optimization_base;
pub mod layout_optimization_boundary_grid;
pub mod layout_optimization_gridded;
pub mod layout_optimization_pyoptsparse;
pub mod layout_optimization_random_search;
pub mod layout_optimization_scipy;
pub mod slsqp;

// Re-exports from submodules
pub use layout_optimization_base::{
    calculate_pairwise_distances, generate_grid_points, is_point_in_boundary, load_optimization,
    load_optimization_result, save_optimization, save_optimization_result, Boundary,
    LayoutOptimizationConfig, LayoutOptimizationResult, LayoutOptimizer, OptimizationConfigFile,
    OptimizationType,
};

pub use layout_optimization_gridded::LayoutOptimizationGridded;
pub use layout_optimization_random_search::LayoutOptimizationRandomSearch;

pub use layout_optimization_scipy::LayoutOptimizationScipy;

pub use layout_optimization_pyoptsparse::{
    LayoutOptimizationGoldenSection, LayoutOptimizationPyOptSparse,
};

pub use layout_optimization_boundary_grid::{
    GridOptimizationConfig, LayoutOptimizationBoundaryGrid, LayoutOptimizationMixedInteger,
};
