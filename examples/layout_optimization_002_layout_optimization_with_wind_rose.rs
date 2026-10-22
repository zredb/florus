//! Example: Layout optimization with heterogeneous inflow
//!
//! This example shows a layout optimization using the geometric yaw option.
//! It combines elements of layout optimization and heterogeneous inflow for
//! demonstrative purposes.
//!
//! Heterogeneity in the inflow provides the necessary driver for coupled yaw
//! and layout optimization to be worthwhile. First, a layout optimization is
//! run without coupled yaw optimization; then a coupled optimization is run to
//! show the benefits of coupled optimization when flows are heterogeneous.
//!
//! Scenario aligned 1:1 with
//! examples_layout_optimization/002_optimize_layout_with_heterogeneity.py.
//!
//! Note: Python optimizes the layout with scipy's SLSQP solver, while
//! FLORUS uses a pattern search, so the optimal layouts are not identical
//! (the AEP gains agree closely).

use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::HeterogeneousInflowConfig;
use florus::optimization::layout_optimization::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig,
};
use florus::optimization::layout_optimization::layout_optimization_scipy::geometric_yaw_angles;
use florus::optimization::layout_optimization::{LayoutOptimizationScipy, LayoutOptimizer};
use florus::types::{Array1, Array2};
use florus::wind_data::WindRose;

fn main() -> anyhow::Result<()> {
    // Initialize FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Setup 2 wind directions (due east and due west)
    // and 1 wind speed with uniform probability
    let wind_directions = Array1::from_vec(vec![90.0, 270.0]);
    let n_wds = wind_directions.len();
    let wind_speeds = Array1::from_vec(vec![8.0]);

    // Shape frequency distribution to match number of wind directions and wind speeds
    let freq_table = Array2::from_elem((n_wds, 1), 1.0 / n_wds as f64);

    // The boundaries for the turbines, specified as vertices
    let d = 126.0; // rotor diameter for the NREL 5MW
    let size_d = 12.0;
    let boundaries = Boundary::Polygon(vec![
        (0.0, 0.0),
        (size_d * d, 0.0),
        (size_d * d, 0.1),
        (0.0, 0.1),
        (0.0, 0.0),
    ]);

    // Set turbine locations to 4 turbines at corners of the rectangle
    // (optimal without flow heterogeneity)
    let layout_x = Array1::from_vec(vec![0.1, 0.3 * size_d * d, 0.6 * size_d * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);

    // Generate exaggerated heterogeneous inflow (same for all wind directions)
    let speed_multipliers = Array2::from_shape_fn((n_wds, 4), |(_, j)| [0.5, 1.0, 0.5, 1.0][j]);
    let x_locs = Array1::from_vec(vec![0.0, size_d * d, 0.0, size_d * d]);
    let y_locs = Array1::from_vec(vec![-d, -d, d, d]);

    // The configuration of the heterogeneous inflow (Python's
    // `heterogeneous_inflow_config_by_wd`)
    let heterogeneous_config = HeterogeneousInflowConfig {
        x: x_locs.clone(),
        y: y_locs.clone(),
        z: None,
        wind_speeds: None,
        wind_directions: Some(wind_directions.clone()),
        speed_multipliers: speed_multipliers.clone(),
    };

    // Establish a WindRose object
    let wind_rose = WindRose::new(
        wind_directions.clone(),
        wind_speeds,
        Array2::from_elem((n_wds, 1), 0.06),
        Some(freq_table),
        None,
        false,
        None,
        None,
    )?;

    fmodel.set_wind_data(&wind_rose)?;
    fmodel.set_heterogeneous_inflow_config(heterogeneous_config)?;
    fmodel.set_layout(&layout_x, &layout_y)?;

    // Setup and solve the layout optimization problem without heterogeneity
    let config = LayoutOptimizationConfig {
        boundaries: boundaries.clone(),
        min_dist: Some(2.0 * d),
        ..Default::default()
    };
    let mut layout_opt =
        LayoutOptimizationScipy::new(&fmodel, boundaries.clone())?.with_config(config);

    // Run the optimization
    let solution = layout_opt.optimize()?;

    // Get the resulting improvement in AEP
    println!("... calcuating improvement in AEP");
    fmodel.run()?;
    let base_aep = fmodel.get_expected_farm_power(None, None)? * 8760.0 / 1e6;
    fmodel.set_layout(&solution.x, &solution.y)?;
    fmodel.run()?;
    let opt_aep = fmodel.get_expected_farm_power(None, None)? * 8760.0 / 1e6;

    let percent_gain = 100.0 * (opt_aep - base_aep) / base_aep;

    // Print the results
    println!("Optimal layout: x={:?}, y={:?}", solution.x, solution.y);
    println!(
        "Optimal layout improves AEP by {percent_gain:.1}% from {base_aep:.1} MWh to {opt_aep:.1} MWh"
    );

    // Rerun the layout optimization with geometric yaw enabled
    println!("\nReoptimizing with geometric yaw enabled.");
    fmodel.set_layout(&layout_x, &layout_y)?;
    let config = LayoutOptimizationConfig {
        boundaries: boundaries.clone(),
        min_dist: Some(2.0 * d),
        enable_geometric_yaw: true,
        ..Default::default()
    };
    let mut layout_opt =
        LayoutOptimizationScipy::new(&fmodel, boundaries.clone())?.with_config(config);

    // Run the optimization
    let solution = layout_opt.optimize()?;

    // Get the resulting improvement in AEP
    println!("... calcuating improvement in AEP");

    // The base AEP is computed with zero yaw angles (Python's example)
    fmodel.set_yaw_angles(Array2::zeros((n_wds, layout_x.len())))?;
    fmodel.run()?;
    let base_aep = fmodel.get_expected_farm_power(None, None)? * 8760.0 / 1e6;

    // The optimized layout with its geometric yaw angles
    fmodel.set_layout(&solution.x, &solution.y)?;
    fmodel.run()?;
    let yaw_angles = geometric_yaw_angles(&mut fmodel)?;
    fmodel.set_yaw_angles(yaw_angles.clone())?;
    fmodel.run()?;
    let opt_aep = fmodel.get_expected_farm_power(None, None)? * 8760.0 / 1e6;

    let percent_gain = 100.0 * (opt_aep - base_aep) / base_aep;

    // Print the results
    println!("Optimal layout: x={:?}, y={:?}", solution.x, solution.y);
    println!(
        "Optimal layout improves AEP by {percent_gain:.1}% from {base_aep:.1} MWh to {opt_aep:.1} MWh"
    );

    // The geometric yaw angles of the optimized layout
    println!(
        "Turbine geometric yaw angles for wind direction {:.2} and wind speed {:.2} m/s: {:?}",
        wind_directions[1],
        8.0,
        yaw_angles.row(1).to_vec()
    );

    Ok(())
}
