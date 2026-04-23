//! Example: Layout Optimization
//!
//! This example demonstrates wind farm layout optimization.
//!
//! Corresponds to: examples_layout_optimization/001_layout_optimization.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Layout Optimization ===\n");

    println!("Wind Farm Layout Optimization:\\n");
    
    println!("Objective:");
    println!("  - Maximize Annual Energy Production (AEP)");
    println!("  - Minimize wake losses");
    println!("  - Optimize turbine placement\\n");

    println!("Optimization Variables:");
    println!("  - Turbine x-coordinates");
    println!("  - Turbine y-coordinates");
    println!("  - Number of turbines (optional)\\n");

    println!("Constraints:");
    println!("  1. Boundary constraints:");
    println!("     - Farm boundary polygon");
    println!("     - Exclusion zones");
    println!("     - Setback distances\\n");

    println!("  2. Spacing constraints:");
    println!("     - Minimum turbine spacing");
    println!("     - Avoid excessive clustering");
    println!("     - Maintenance access\\n");

    println!("  3. Regulatory constraints:");
    println!("     - Environmental restrictions");
    println!("     - Property boundaries");
    println!("     - Aviation requirements\\n");

    println!("Optimization Methods:");
    println!("  - Gradient-based (fast, local optimum)");
    println!("  - Genetic algorithms (slow, global search)");
    println!("  - Random search (simple baseline)");
    println!("  - Particle swarm optimization\\n");

    println!("Benefits:");
    println!("  - 5-15% AEP improvement typical");
    println!("  - Reduced wake losses");
    println!("  - Better land utilization");
    println!("  - Higher revenue\\n");

    println!("=== Example Complete ===");
    println!("Note: Full implementation requires optimization framework.");
    Ok(())
}
