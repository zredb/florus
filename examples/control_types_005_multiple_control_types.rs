//! Example: Multiple Control Types
//!
//! This example demonstrates combining multiple control types.
//!
//! Corresponds to: examples_control_types/005_multiple_control_types.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Multiple Control Types ===\n");

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

    println!("Combined Control Strategy:\\n");
    println!("Control Type     | Purpose");
    println!("{}|{}", "-".repeat(17), "-".repeat(40));
    println!("{:<16} | {}", "Yaw angles", "Wake steering (horizontal)");
    println!("{:<16} | {}", "Tilt angles", "Wake deflection (vertical)");
    println!("{:<16} | {}", "Power limits", "Curtailment and load control");
    println!("{:<16} | {}", "AWC modes", "Active wake mixing");

    println!("\nBenefits of combined control:");
    println!("  - Multi-dimensional wake manipulation");
    println!("  - Optimized farm-level performance");
    println!("  - Reduced structural loads");
    println!("  - Improved grid integration");

    // Demonstrate yaw control as example
    fmodel.set_yaw_angles(ndarray::arr2(&[[15.0, 0.0]]))?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers();
    
    println!("\nExample: T0 yawed 15°");
    println!("  T0: {:.2} kW", powers[[0, 0]] / 1000.0);
    println!("  T1: {:.2} kW", powers[[0, 1]] / 1000.0);
    println!("  Total: {:.2} kW", (powers[[0, 0]] + powers[[0, 1]]) / 1000.0);

    println!("\n=== Example Complete ===");
    Ok(())
}
