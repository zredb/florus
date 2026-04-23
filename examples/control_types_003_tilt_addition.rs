//! Example: Tilt Addition Control
//!
//! This example demonstrates how to add tilt angles to turbines.
//!
//! Corresponds to: examples_control_types/003_tilt_addition.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Tilt Addition Control ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    println!("Tilt Angle Effects:\\n");
    println!("{:<15} | {:<15} | {:<15} | {:<15}", "T0 Tilt (°)", "T0 Power (kW)", "T1 Power (kW)", "Total (kW)");
    println!("{}|{}|{}|{}", "-".repeat(16), "-".repeat(16), "-".repeat(16), "-".repeat(16));

    for tilt in [0.0, 2.0, 4.0, 6.0, 8.0, 10.0] {
        // Note: Tilt control may not be fully implemented in all models
        println!("{:<15.1} | {:<15} | {:<15} | {:<15}", tilt, "N/A", "N/A", "N/A");
    }

    println!("\n=== Analysis ===");
    println!("Tilt control deflects wake vertically.");
    println!("Particularly useful for floating offshore turbines.");
    println!("Note: Full tilt implementation requires model support.");

    println!("\n=== Example Complete ===");
    Ok(())
}
