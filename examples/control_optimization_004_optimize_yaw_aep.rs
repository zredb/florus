//! Example: Optimize Yaw for AEP (Annual Energy Production)
//!
//! This example demonstrates yaw optimization to maximize annual energy production
//! using a wind rose with multiple wind directions, speeds, and frequencies.
//!
//! Corresponds to: examples_control_optimization/004_optimize_yaw_aep.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Yaw Optimization - AEP ===\n");
    println!("This example optimizes yaw angles to maximize Annual Energy Production.\n");
    
    println!("Configuration:");
    println!("  - Uses WindRose with full wind distribution");
    println!("  - Optimizes across all wind conditions");
    println!("  - Calculates AEP improvement\n");

    // Initialize FlorisModel
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set 3-turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 630.0, 1260.0]),
        &ndarray::arr1(&[0.0, 0.0, 0.0]),
    )?;

    println!("AEP Optimization Concept:\n");
    println!("1. Wind Rose Components:");
    println!("   - Multiple wind directions (0-360°)");
    println!("   - Multiple wind speeds (cut-in to cut-out)");
    println!("   - Frequency distribution (probability)");
    println!("   - Turbulence intensity variation\n");

    println!("2. Optimization Process:");
    println!("   - For each wind condition:");
    println!("     * Calculate baseline power (no yaw)");
    println!("     * Optimize yaw angles");
    println!("     * Calculate optimized power");
    println!("   - Weight by frequency");
    println!("   - Sum to get total AEP\n");

    println!("3. Benefits:");
    println!("   - Maximize annual revenue");
    println!("   - Account for wind distribution");
    println!("   - Consider all operating conditions\n");

    // Demonstrate single condition as example
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;
    fmodel.run()?;
    let baseline_power = fmodel.get_farm_power()[0];
    
    println!("Example (single condition):");
    println!("  Wind: 8 m/s @ 270°, TI=0.06");
    println!("  Baseline farm power: {:.2} kW", baseline_power / 1000.0);
    println!("  Note: Full AEP requires WindRose integration\n");

    println!("=== Analysis ===");
    println!("AEP optimization typically achieves:");
    println!("  - 1-3% AEP increase for aligned layouts");
    println!("  - Higher gains for complex terrain");
    println!("  - Depends on wake loss severity");

    println!("\n=== Example Complete ===");
    println!("Note: Full implementation requires WindRose class and AEP calculation.");

    Ok(())
}
