/// Layout optimization module for wind farm turbine positioning
///
/// Provides optimization algorithms for wind farm layout including:
/// - Random search optimization (genetic algorithm-style)
/// - Scipy-based coordinate descent optimization
/// - PyOptSparse-style gradient-based optimization (argmin)
/// - Grid-based exhaustive search optimization
/// - Mixed integer optimization (grid + continuous refinement)
///
/// This module corresponds to floris/optimization/layout_optimization/ in Python FLORIS v4.6

// Module declarations
pub mod layout_optimization_base;
pub mod layout_optimization_random_search;
pub mod layout_optimization_scipy;
pub mod layout_optimization_pyoptsparse;
pub mod layout_optimization_boundary_grid;

// Re-exports from submodules
pub use layout_optimization_base::{
    Boundary,
    LayoutOptimizationConfig,
    LayoutOptimizationResult,
    LayoutOptimizer,
    OptimizationConfigFile,
    OptimizationType,
    is_point_in_boundary,
    generate_grid_points,
    calculate_pairwise_distances,
    load_optimization,
    save_optimization,
    load_optimization_result,
    save_optimization_result,
};

pub use layout_optimization_random_search::{
    LayoutOptimizationRandomSearch,
};

pub use layout_optimization_scipy::{
    LayoutOptimizationScipy,
};

pub use layout_optimization_pyoptsparse::{
    LayoutOptimizationPyOptSparse,
    LayoutOptimizationGoldenSection,
};

pub use layout_optimization_boundary_grid::{
    LayoutOptimizationBoundaryGrid,
    LayoutOptimizationMixedInteger,
    GridOptimizationConfig,
};
