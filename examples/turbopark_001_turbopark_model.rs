//! Example: TurboPark Wake Model
//!
//! This example demonstrates the TurboPark wake model.
//!
//! Corresponds to: examples_turbopark/001_turbopark_model.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== TurboPark Wake Model ===\n");
    println!("This example demonstrates the TurboPark wake model.\n");

    // Try to load TurboPark configuration
    match FlorisModel::from_file("examples/inputs/turbopark.yaml") {
        Ok(fmodel) => {
            println!("TurboPark model loaded successfully!\n");
            
            fmodel.set_layout(
                &ndarray::arr1(&[0.0, 630.0]),
                &ndarray::arr1(&[0.0, 0.0]),
            )?;

            fmodel.set_wind_conditions(
                ndarray::arr1(&[8.0]),
                ndarray::arr1(&[270.0]),
                ndarray::arr1(&[0.06]),
            )?;

            fmodel.run()?;
            let powers = fmodel.get_turbine_powers();
            
            println!("Results:");
            println!("  T0: {:.2} kW", powers[[0, 0]] / 1000.0);
            println!("  T1: {:.2} kW", powers[[0, 1]] / 1000.0);
            println!("  Total: {:.2} kW\n", (powers[[0, 0]] + powers[[0, 1]]) / 1000.0);
        }
        Err(e) => {
            println!("Note: TurboPark configuration not available.\n");
            println!("Error: {}\n", e);
        }
    }

    println!("TurboPark Model Features:\n");
    println!("1. Gaussian-based wake model");
    println!("   - Smooth velocity deficit profiles");
    println!("   - Analytical solutions");
    println!("   - Fast computation\n");

    println!("2. Wake superposition");
    println!("   - Sum-of-squares method");
    println!("   - Handles multiple wakes");
    println!("   - Physically consistent\n");

    println!("3. Applications");
    println!("   - Large wind farms");
    println!("   - Layout optimization");
    println!("   - AEP calculations\n");

    println!("=== Example Complete ===");
    Ok(())
}
