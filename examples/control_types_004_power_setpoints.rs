//! Example: Power Setpoints Control
//!
//! This example demonstrates how to set power limits for turbines.
//!
//! Corresponds to: examples_control_types/004_power_setpoints.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Power Setpoints Control ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    fmodel.set_wind_conditions(
        ndarray::arr1(&[12.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    fmodel.run()?;
    let powers = fmodel.get_turbine_powers();

    println!("Baseline (no curtailment):");
    println!("  T0: {:.2} kW", powers[[0, 0]] / 1000.0);
    println!("  T1: {:.2} kW", powers[[0, 1]] / 1000.0);
    println!("  Total: {:.2} kW\\n", (powers[[0, 0]] + powers[[0, 1]]) / 1000.0);

    println!("Power setpoints allow curtailment:");
    println!("  - Limit maximum power output");
    println!("  - Useful for grid constraints");
    println!("  - Can reduce mechanical loads");
    println!("  - Enable coordinated farm control");

    println!("\n=== Example Complete ===");
    println!("Note: Full power setpoint control requires operation model support.");
    Ok(())
}
