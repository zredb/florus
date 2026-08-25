//! Example: Disabling Turbines
//!
//! Demonstrates shutting down specific turbines during simulation using
//! the "mixed" operation model.
//!
//! Corresponds to: examples_control_types/002_disable_turbines.py

use florus::{FlorisModel, Result};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Disabling Turbines Example ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 3-turbine aligned layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0, 1000.0]),
        &ndarray::arr1(&[0.0, 0.0, 0.0]),
    )?;

    // Change to the mixed operation model
    fmodel.set_operation_model("mixed")?;

    // 2 identical wind conditions
    let wind_speeds = Array1::from_vec(vec![8.0, 8.0]);
    let wind_directions = Array1::from_vec(vec![270.0, 270.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06, 0.06]);
    fmodel.set_wind_conditions(wind_speeds, wind_directions, turbulence_intensities)?;

    // Disable T0 & T1 for the second findex
    let mut disable_turbines = ndarray::Array2::from_elem((2, 3), false);
    disable_turbines[[1, 0]] = true;
    disable_turbines[[1, 1]] = true;

    fmodel.set_operation(None, None, None, None, None, Some(disable_turbines))?;
    fmodel.run()?;

    let turbine_powers = fmodel.get_turbine_powers() / 1000.0;
    let effective_wind_speeds = fmodel.turbine_average_velocities();

    println!("Turbine Powers (kW):");
    println!("{:<12} | {:<12} | {:<12} | {:<12}", "Condition", "T0", "T1", "T2");
    println!("{}", "-".repeat(52));
    for fi in 0..2 {
        let label = if fi == 0 { "All on" } else { "T0&T1 off" };
        println!(
            "{:<12} | {:<12.2} | {:<12.2} | {:<12.2}",
            label, turbine_powers[[fi, 0]], turbine_powers[[fi, 1]], turbine_powers[[fi, 2]]
        );
    }

    println!("\nEffective Wind Speeds (m/s):");
    println!("{:<12} | {:<12} | {:<12} | {:<12}", "Condition", "T0", "T1", "T2");
    println!("{}", "-".repeat(52));
    for fi in 0..2 {
        let label = if fi == 0 { "All on" } else { "T0&T1 off" };
        println!(
            "{:<12} | {:<12.2} | {:<12.2} | {:<12.2}",
            label, effective_wind_speeds[[fi, 0]], effective_wind_speeds[[fi, 1]], effective_wind_speeds[[fi, 2]]
        );
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
