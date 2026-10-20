//! Example: Heterogeneous Inflow for single case
//!
//! This example illustrates how to set up a heterogeneous inflow condition:
//! it changes the wind farm layout, sets a single wind condition, defines a
//! heterogeneous inflow (1.25x ambient for the upper turbines) and prints the
//! resulting turbine powers.
//!
//! Scenario aligned 1:1 with
//! examples_heterogeneous/001_heterogeneous_inflow_single.py.

use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::HeterogeneousInflowConfig;
use florus::types::{Array1, Array2};
use ndarray::Array;

fn main() -> anyhow::Result<()> {
    // Initialize FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Change the layout to a 4 turbine layout in a box
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 0.0, 500.0, 500.0]),
        &Array1::from_vec(vec![0.0, 500.0, 0.0, 500.0]),
    )?;

    // Set FLORIS to run for a single condition
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;

    // Define the speed-ups of the heterogeneous inflow, and their locations.
    // Heterogeneity is only applied within the bounds of the given points; the
    // inflow is 1.25x the ambient wind speed for the upper turbines at y = 500 m.
    let speed_ups = Array2::from_shape_vec((1, 4), vec![1.0, 1.25, 1.0, 1.25])?;
    let x_locs = Array1::from_vec(vec![-500.0, -500.0, 1000.0, 1000.0]);
    let y_locs = Array1::from_vec(vec![-500.0, 1000.0, -500.0, 1000.0]);

    // Set the heterogeneous inflow configuration
    fmodel.set_heterogeneous_inflow_config(HeterogeneousInflowConfig {
        x: x_locs,
        y: y_locs,
        z: None,
        wind_speeds: None,
        wind_directions: None,
        speed_multipliers: speed_ups,
    })?;

    // Run the FLORIS simulation
    fmodel.run()?;

    // Get the power output of the turbines
    let turbine_powers = fmodel.get_turbine_powers()? / 1000.0;

    println!("Turbine 0 power = {:.1} kW", turbine_powers[[0, 0]]);
    println!("Turbine 1 power = {:.1} kW", turbine_powers[[0, 1]]);
    println!("Turbine 2 power = {:.1} kW", turbine_powers[[0, 2]]);
    println!("Turbine 3 power = {:.1} kW", turbine_powers[[0, 3]]);

    // Keep the ndarray import used (the speeds are set through the API above)
    let _ = Array::from_elem(1, 0.0);

    Ok(())
}
