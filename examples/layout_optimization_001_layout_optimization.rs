//! Example: Optimize Layout
//!
//! This example shows a simple layout optimization, optimizing for both
//! annual energy production (AEP) and annual value production (AVP).
//!
//! First, a 4 turbine array is optimized such that the layout of the
//! turbines produces the highest AEP based on the given wind resource. The
//! turbines are constrained to a square boundary and a random wind resource
//! is supplied.
//!
//! Next, with the same boundary, the same 4 turbine array is optimized to
//! maximize AVP instead of AEP, using the value table defined in the
//! WindRose object, where value represents the value of the energy produced
//! for a given wind condition.
//!
//! Scenario aligned 1:1 with examples_layout_optimization/001_optimize_layout.py.
//!
//! Note: Python optimizes the layout with scipy's SLSQP solver, while
//! FLORUS uses a coordinate descent. The optimal layouts are therefore not
//! identical, but both find the same optimum (the turbines are pushed to the
//! corners of the boundary), so the AEP and AVP agree closely.

use florus::floris_model::FlorisModel;
use florus::numpy_random::NumpyRandomState;
use florus::optimization::layout_optimization::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig,
};
use florus::optimization::layout_optimization::{LayoutOptimizationScipy, LayoutOptimizer};
use florus::types::{Array1, Array2};
use florus::wind_data::WindRose;

fn main() -> anyhow::Result<()> {
    // Initialize the FLORIS interface
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Setup 72 wind directions with 1 wind speed and a frequency distribution
    let wind_directions = Array1::from_iter((0..72).map(|i| i as f64 * 5.0));
    let wind_speeds = Array1::from_vec(vec![8.0]);

    // Shape a random frequency distribution to match the number of wind
    // directions and wind speeds (numpy: np.random.seed(1) and
    // np.random.randn(72))
    let mut rng = NumpyRandomState::new(1);
    let mut frequencies = rng.standard_normal_vec(72);
    frequencies.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let frequencies: Vec<f64> = frequencies.iter().map(|v| v.abs()).collect();
    let total: f64 = frequencies.iter().sum();
    let freq_table = Array2::from_shape_fn((72, 1), |(i, _)| frequencies[i] / total);

    // Define the value table such that the value of the energy produced is
    // significantly higher when the wind direction is close to the north or
    // south, and zero when the wind is from the east or west. Here, value is
    // given a mean value of 25 USD/MWh.
    let raw_values: Vec<f64> = wind_directions
        .iter()
        .map(|&wd| (0.5 + 0.5 * (2.0 * wd.to_radians()).cos()).powi(10))
        .collect();
    let mean_value: f64 = raw_values.iter().sum::<f64>() / raw_values.len() as f64;
    let value_table = Array2::from_shape_fn((72, 1), |(i, _)| 25.0 * raw_values[i] / mean_value);

    // Establish a WindRose object
    let wind_rose = WindRose::new(
        wind_directions,
        wind_speeds,
        Array2::from_elem((72, 1), 0.06),
        Some(freq_table),
        Some(value_table),
        false,
        None,
        None,
    )?;

    fmodel.set_wind_data(&wind_rose)?;

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

    // Setup the optimization problem to maximize AEP instead of value
    let config = LayoutOptimizationConfig {
        boundaries: boundaries.clone(),
        ..Default::default()
    };
    let mut layout_opt =
        LayoutOptimizationScipy::new(&fmodel, boundaries.clone())?.with_config(config);

    // Run the optimization
    let solution = layout_opt.optimize()?;

    // Get the resulting improvement in AEP
    println!("... calculating improvement in AEP");
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

    // reset to the original layout
    fmodel.set_layout(&layout_x, &layout_y)?;

    // Now set up the optimization problem to maximize annual value production
    // (AVP) using the value table provided in the WindRose object.
    let config = LayoutOptimizationConfig {
        boundaries: boundaries.clone(),
        use_value: true,
        ..Default::default()
    };
    let mut layout_opt =
        LayoutOptimizationScipy::new(&fmodel, boundaries.clone())?.with_config(config);

    // Run the optimization
    let solution = layout_opt.optimize()?;

    // Get the resulting improvement in AVP
    println!("... calculating improvement in annual value production (AVP)");
    fmodel.run()?;
    let base_avp = fmodel.get_expected_farm_value(None, None, None)? * 8760.0 / 1e6;
    fmodel.set_layout(&solution.x, &solution.y)?;
    fmodel.run()?;
    let opt_avp = fmodel.get_expected_farm_value(None, None, None)? * 8760.0 / 1e6;

    let percent_gain = 100.0 * (opt_avp - base_avp) / base_avp;

    // Print the results
    println!("Optimal layout: x={:?}, y={:?}", solution.x, solution.y);
    println!(
        "Optimal layout improves AVP by {percent_gain:.1}% from {base_avp:.1} dollars to {opt_avp:.1} dollars"
    );

    Ok(())
}
