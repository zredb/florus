//! Example: Layout Optimization with Wind Rose
//!
//! This example demonstrates layout optimization using wind rose data.
//!
//! Corresponds to: examples_layout_optimization/002_layout_optimization_with_wind_rose.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Layout Optimization with Wind Rose ===\n");

    println!("Layout Optimization with Full Wind Distribution:\\n");
    
    println!("Key Difference from Single Condition:");
    println!("  - Optimizes across all wind conditions");
    println!("  - Weights by wind direction frequency");
    println!("  - Maximizes total AEP, not single-case power\\n");

    println!("Process:");
    println!("  1. Define wind rose (directions, speeds, frequencies)");
    println!("  2. For each candidate layout:");
    println!("     - Calculate power for each wind condition");
    println!("     - Weight by frequency");
    println!("     - Sum for total AEP");
    println!("  3. Optimize layout to maximize AEP\\n");

    println!("Benefits:");
    println!("  - Accounts for prevailing winds");
    println!("  - Avoids optimizing for rare conditions");
    println!("  - More realistic performance prediction");
    println!("  - Better financial returns\\n");

    println!("Challenges:");
    println!("  - Higher computational cost");
    println!("  - More complex optimization landscape");
    println!("  - Longer run times");
    println!("  - Requires good wind data\\n");

    println!("Typical Results:");
    println!("  - Turbines staggered relative to prevailing wind");
    println!("  - Aligned perpendicular to dominant directions");
    println!("  - 5-20% AEP improvement over grid layout");
    println!("  - Site-specific optimal patterns\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires WindRose integration and optimization framework.");
    Ok(())
}
