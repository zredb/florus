//! Example: Layout Optimization - Gridded
//!
//! This example demonstrates gridded layout optimization approach.
//!
//! Corresponds to: examples_layout_optimization/003_layout_optimization_gridded.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Layout Optimization - Gridded ===\n");

    println!("Gridded Layout Optimization:\\n");
    
    println!("Concept:");
    println!("  - Define a grid of possible turbine locations");
    println!("  - Select optimal subset of grid points");
    println!("  - Binary optimization problem\\n");

    println!("Advantages:");
    println!("  - Discrete search space");
    println!("  - Easier to implement constraints");
    println!("  - Good for irregular boundaries");
    println!("  - Can use combinatorial optimization\\n");

    println!("Methods:");
    println!("  - Greedy selection");
    println!("  - Genetic algorithms");
    println!("  - Mixed-integer programming");
    println!("  - Sequential addition/removal\\n");

    println!("Grid Resolution Trade-offs:");
    println!("  - Fine grid: More options, slower");
    println!("  - Coarse grid: Faster, less flexible");
    println!("  - Typical: 100-500m spacing\\n");

    println!("Applications:");
    println!("  - Complex terrain");
    println!("  - Irregular farm boundaries");
    println!("  - Exclusion zones");
    println!("  - Phased development\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires grid-based optimization framework.");
    Ok(())
}
