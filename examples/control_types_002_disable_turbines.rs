//! Example: Disabling turbines
//!
//! This example demonstrates the mixed operation model with disabled
//! turbines: the first two turbines of a three-turbine farm are shut down for
//! the second wind condition, which changes the wake of the remaining
//! turbine.
//!
//! Scenario aligned 1:1 with
//! examples_control_types/002_disable_turbines.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Change to the mixed model turbine
    fmodel.set_operation_model("mixed")?;

    // A wind farm of 3 aligned wind turbines
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // Two identical wind conditions
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0, 8.0]),
        Array1::from_vec(vec![270.0, 270.0]),
        Array1::from_vec(vec![0.06, 0.06]),
    )?;

    // Shut down the first 2 turbines for the second findex
    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((2, 3), false);
    disable_turbines[[1, 0]] = true;
    disable_turbines[[1, 1]] = true;

    fmodel.set_operation(None, None, None, None, None, Some(disable_turbines))?;
    fmodel.run()?;

    let turbine_powers = fmodel.get_turbine_powers()?;
    let effective_wind_speeds = fmodel.turbine_average_velocities();

    println!("Turbine powers [kW]:");
    println!(
        "  All on:          {:.2} {:.2} {:.2}",
        turbine_powers[[0, 0]] / 1000.0,
        turbine_powers[[0, 1]] / 1000.0,
        turbine_powers[[0, 2]] / 1000.0
    );
    println!(
        "  T0 & T1 disabled: {:.2} {:.2} {:.2}",
        turbine_powers[[1, 0]] / 1000.0,
        turbine_powers[[1, 1]] / 1000.0,
        turbine_powers[[1, 2]] / 1000.0
    );
    println!("Effective wind speeds [m/s]:");
    println!(
        "  All on:          {:.4} {:.4} {:.4}",
        effective_wind_speeds[[0, 0]],
        effective_wind_speeds[[0, 1]],
        effective_wind_speeds[[0, 2]]
    );
    println!(
        "  T0 & T1 disabled: {:.4} {:.4} {:.4}",
        effective_wind_speeds[[1, 0]],
        effective_wind_speeds[[1, 1]],
        effective_wind_speeds[[1, 2]]
    );

    Ok(())
}
