//! Example: Load Optimization
//!
//! This example demonstrates wind farm load optimization.
//!
//! Corresponds to: examples_load_optimization/001_load_optimization.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Load Optimization ===\n");

    println!("Wind Farm Load Optimization:\\n");
    
    println!("Objective:");
    println!("  - Minimize structural loads on turbines");
    println!("  - Extend turbine lifetime");
    println!("  - Reduce maintenance costs");
    println!("  - Balance with power production\\n");

    println!("Load Types:");
    println!("  1. Fatigue loads:");
    println!("     - Cyclic loading over time");
    println!("     - Damage accumulation");
    println!("     - Key for lifetime estimation\\n");

    println!("  2. Extreme loads:");
    println!("     - Maximum instantaneous loads");
    println!("     - Design limit states");
    println!("     - Safety margins\\n");

    println!("  3. Component-specific:");
    println!("     - Blade root bending moments");
    println!("     - Tower base moments");
    println!("     - Shaft torque\\n");

    println!("Optimization Strategies:");
    println!("  - Derate upstream turbines");
    println!("  - Adjust yaw angles for load reduction");
    println!("  - Coordinate control across farm");
    println!("  - Trade-off: power vs. loads\\n");

    println!("Benefits:");
    println!("  - Extended turbine life (10-20%)");
    println!("  - Reduced O&M costs");
    println!("  - Lower insurance premiums");
    println!("  - Improved reliability\\n");

    println!("Challenges:");
    println!("  - Complex load modeling");
    println!("  - Multi-objective optimization");
    println!("  - Uncertainty in load predictions");
    println!("  - Validation difficulty\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires load modeling and multi-objective optimization.");
    Ok(())
}
