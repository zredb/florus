//! Example: Layout optimization with genetic random search
//!
//! This example shows a layout optimization using the genetic random search
//! algorithm. It provides options for the users to try different distance
//! probability mass functions for the random search perturbations.
//!
//! Scenario aligned 1:1 with
//! examples_layout_optimization/003_genetic_random_search.py.
//!
//! Note: Python runs the search until a wall-clock deadline (10 s per
//! individual and 60 s in total), so its result is not reproducible; the
//! Rust example runs a fixed number of generations and iterations instead.
//! The random wind resource is drawn with the numpy compatible generator and
//! is therefore identical to Python's.

use florus::floris_model::FlorisModel;
use florus::numpy_random::NumpyRandomState;
use florus::optimization::layout_optimization::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig,
};
use florus::optimization::layout_optimization::{LayoutOptimizationRandomSearch, LayoutOptimizer};
use florus::types::{Array1, Array2};
use florus::wind_data::WindRose;

fn main() -> anyhow::Result<()> {
    // Set up FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Setup 72 wind directions with a wind speed and a frequency
    // distribution (numpy: np.random.seed(1) and np.random.randn(72))
    let wind_directions = Array1::from_iter((0..72).map(|i| i as f64 * 5.0));
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let mut rng = NumpyRandomState::new(1);
    let mut frequencies = rng.standard_normal_vec(72);
    frequencies.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let frequencies: Vec<f64> = frequencies.iter().map(|v| v.abs()).collect();
    let total: f64 = frequencies.iter().sum();
    let freq_table = Array2::from_shape_fn((72, 1), |(i, _)| frequencies[i] / total);

    fmodel.set_wind_data(&WindRose::new(
        wind_directions,
        wind_speeds,
        Array2::from_elem((72, 1), 0.06),
        Some(freq_table),
        None,
        false,
        None,
        None,
    )?)?;

    // Set the boundaries
    // The boundaries for the turbines, specified as vertices
    let boundaries = Boundary::Polygon(vec![
        (0.0, 0.0),
        (0.0, 1000.0),
        (1000.0, 1000.0),
        (1000.0, 0.0),
        (0.0, 0.0),
    ]);

    // Set turbine locations to 4 turbines in a rectangle
    let d = 126.0; // rotor diameter for the NREL 5MW
    let layout_x = Array1::from_vec(vec![0.0, 0.0, 6.0 * d, 6.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 4.0 * d, 0.0, 4.0 * d]);
    fmodel.set_layout(&layout_x, &layout_y)?;
    fmodel.run()?;

    // Perform the optimization
    let config = LayoutOptimizationConfig {
        boundaries: boundaries.clone(),
        min_dist: Some(5.0 * d),
        // Python runs the search until a wall-clock deadline (10 s per
        // individual); the Rust example runs a fixed number of iterations
        max_iterations: 3000,
        n_generations: 3,
        ..Default::default()
    };
    let mut layout_opt = LayoutOptimizationRandomSearch::new(&fmodel, boundaries, 4)?
        .with_config(config)
        .with_random_seed(1);

    // The description of the problem (Python's `describe`)
    let initial_aep = fmodel.get_expected_farm_power(None, None)? * 8760.0 / 1e9;
    println!("Random Layout Optimization");
    println!(
        "Number of turbines to optimize = {}",
        fmodel.farm().n_turbines()
    );
    println!(
        "Minimum distance between turbines = {:.1} [D], {:.1} [m]",
        5.0,
        5.0 * d
    );
    println!("Number of individuals = 4");
    println!("Initial AEP = {initial_aep:.1} [GWh]");

    // Run the optimization
    let result = layout_opt.optimize()?;

    // The final result (Python's `_finalize_optimization`)
    println!("Optimal layout: x={:?}, y={:?}", result.x, result.y);
    println!(
        "Final AEP = {:.1} [GWh] ({:+.2}%)",
        result.value / 1e9,
        result.improvement_pct
    );

    Ok(())
}
