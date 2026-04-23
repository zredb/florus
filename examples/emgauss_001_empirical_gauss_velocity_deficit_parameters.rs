//! Example: Empirical Gauss Velocity Deficit Parameters
//!
//! This example demonstrates the empirical Gauss wake model parameters.
//!
//! Corresponds to: examples_emgauss/001_empirical_gauss_velocity_deficit_parameters.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Empirical Gauss Velocity Deficit Parameters ===\n");

    let fmodel = FlorisModel::from_file("examples/inputs/emgauss.yaml")?;

    println!("Empirical Gauss Wake Model Parameters:\\n");
    println!("The empirical Gauss model uses measured data to define:");
    println!("  - Velocity deficit profiles");
    println!("  - Wake expansion rates");
    println!("  - Recovery characteristics");

    println!("\nKey Parameters:");
    println!("  - sigma_z0: Initial wake width");
    println!("  - sigma_y0: Initial wake width (lateral)");
    println!("  - kz: Vertical expansion rate");
    println!("  - ky: Lateral expansion rate");
    println!("  - alpha: Deficit decay rate");
    println!("  - beta: Deficit shape parameter");

    println!("\nAdvantages:");
    println!("  - Based on field measurements");
    println!("  - Captures real wake behavior");
    println!("  - More accurate for specific sites");

    println!("\n=== Example Complete ===");
    println!("Note: Full parameter sweep requires emgauss model configuration.");
    Ok(())
}
