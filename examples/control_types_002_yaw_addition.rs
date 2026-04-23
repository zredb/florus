//! Example: Yaw Addition Control
//!
//! This example demonstrates how to add yaw angles to turbines for wake steering.
//!
//! Corresponds to: examples_control_types/002_yaw_addition.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Yaw Addition Control ===\n");

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

    println!("Yaw Angle Sweep:\\n");
    println!("{:<15} | {:<15} | {:<15} | {:<15}", "T0 Yaw (°)", "T0 Power (kW)", "T1 Power (kW)", "Total (kW)");
    println!("{}|{}|{}|{}", "-".repeat(16), "-".repeat(16), "-".repeat(16), "-".repeat(16));

    for yaw in [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0] {
        fmodel.set_yaw_angles(ndarray::arr2(&[[yaw, 0.0]]))?;
        fmodel.run()?;
        let powers = fmodel.get_turbine_powers();
        let t0 = powers[[0, 0]] / 1000.0;
        let t1 = powers[[0, 1]] / 1000.0;
        let total = t0 + t1;
        println!("{:<15.1} | {:<15.2} | {:<15.2} | {:<15.2}", yaw, t0, t1, total);
    }

    println!("\n=== Analysis ===");
    println!("Yaw control redirects the wake away from downstream turbines.");
    println!("Optimal yaw angle balances upstream loss vs downstream gain.");

    println!("\n=== Example Complete ===");
    Ok(())
}
