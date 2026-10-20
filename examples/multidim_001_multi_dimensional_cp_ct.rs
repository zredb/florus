//! Example: Multi-dimensional power/thrust coefficient data
//!
//! This example creates a FLORIS instance with a turbine definition that has
//! a multi-dimensional power/thrust coefficient table (the wave period `Tp`
//! and the wave height `Hs` select the table) and shows the turbine powers
//! for a single wind condition, several wind speeds and several wind
//! directions.
//!
//! Scenario aligned 1:1 with
//! examples_multidim/001_multi_dimensional_cp_ct.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    // Initialize FLORIS with the given input file
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct.yaml")?;

    // Convert to a simple two turbine layout
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;

    // Single wind speed and wind direction
    println!("\n========================= Single Wind Direction and Wind Speed =========================");

    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;

    // Set the yaw angles to 0
    fmodel.set_yaw_angles(Array2::zeros((1, 2)))?;

    // Calculate
    fmodel.run()?;

    // Get the turbine powers
    let turbine_powers = fmodel.get_turbine_powers()? / 1000.0;
    println!("The turbine power matrix should be of dimensions 1 findex X 2 Turbines");
    println!("{turbine_powers:?}");
    println!(
        "Shape:  ({}, {})",
        turbine_powers.shape()[0],
        turbine_powers.shape()[1]
    );

    // Single wind speed and multiple wind directions
    println!("\n========================= Single Wind Direction and Multiple Wind Speeds ===============");

    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array1::from_vec(vec![270.0, 270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06]),
    )?;
    fmodel.set_yaw_angles(Array2::zeros((3, 2)))?;
    fmodel.run()?;

    let turbine_powers = fmodel.get_turbine_powers()? / 1000.0;
    println!("The turbine power matrix should be of dimensions 3 findex X 2 Turbines");
    println!("{turbine_powers:?}");
    println!(
        "Shape:  ({}, {})",
        turbine_powers.shape()[0],
        turbine_powers.shape()[1]
    );

    // Multiple wind speeds and multiple wind directions
    println!("\n========================= Multiple Wind Directions and Multiple Wind Speeds ============");

    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0, 9.0, 10.0, 8.0, 9.0, 10.0, 8.0, 9.0, 10.0]),
        Array1::from_vec(vec![
            260.0, 260.0, 260.0, 270.0, 270.0, 270.0, 280.0, 280.0, 280.0,
        ]),
        Array1::from_vec(vec![0.06; 9]),
    )?;
    fmodel.set_yaw_angles(Array2::zeros((9, 2)))?;
    fmodel.run()?;

    let turbine_powers = fmodel.get_turbine_powers()? / 1000.0;
    println!("The turbine power matrix should be of dimensions 9 WD/WS X 2 Turbines");
    println!("{turbine_powers:?}");
    println!(
        "Shape:  ({}, {})",
        turbine_powers.shape()[0],
        turbine_powers.shape()[1]
    );

    Ok(())
}
