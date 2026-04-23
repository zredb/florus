//! Example: Uncertain Floris Model
//!
//! This example demonstrates uncertainty quantification in FLORIS.
//!
//! Corresponds to: examples_uncertain/001_uncertain_floris_model.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Uncertain Floris Model ===\n");

    println!("Uncertainty Quantification:\\n");
    println!("Real wind farms face uncertainties in:");
    println!("  - Wind direction (±5° typical)");
    println!("  - Wind speed (measurement error)");
    println!("  - Turbulence intensity (spatial variation)");
    println!("  - Turbine performance (degradation)");

    println!("\nUncertainModel Approach:");
    println!("  - Sample multiple scenarios");
    println!("  - Weight by probability");
    println!("  - Compute expected values");
    println!("  - Quantify risk/variance");

    println!("\nApplications:");
    println!("  - Robust yaw optimization");
    println!("  - AEP estimation with confidence");
    println!("  - Risk-aware control strategies");

    println!("\n=== Example Complete ===");
    println!("Note: Full uncertain modeling requires UncertainFlorisModel class.");
    Ok(())
}
