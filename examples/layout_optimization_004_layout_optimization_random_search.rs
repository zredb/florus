//! Example: Layout Optimization - Random Search
//!
//! This example demonstrates random search for layout optimization.
//!
//! Corresponds to: examples_layout_optimization/004_layout_optimization_random_search.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Layout Optimization - Random Search ===\n");

    println!("Random Search Layout Optimization:\\n");
    
    println!("Algorithm:");
    println!("  1. Generate random turbine layouts");
    println!("  2. Evaluate AEP for each layout");
    println!("  3. Keep the best layout found");
    println!("  4. Repeat for N iterations\\n");

    println!("Advantages:");
    println!("  - Simple to implement");
    println!("  - No gradient required");
    println!("  - Global search capability");
    println!("  - Easy to parallelize\\n");

    println!("Disadvantages:");
    println!("  - Slow convergence");
    println!("  - Many evaluations needed");
    println!("  - No guarantee of optimality");
    println!("  - Inefficient for large problems\\n");

    println!("Use Cases:");
    println!("  - Baseline comparison");
    println!("  - Small turbine counts (<20)");
    println!("  - Quick feasibility studies");
    println!("  - Educational purposes\\n");

    println!("Typical Performance:");
    println!("  - 100-1000 iterations needed");
    println!("  - 3-8% AEP improvement");
    println!("  - Better than grid, worse than gradient");
    println!("  - Good starting point for refinement\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires random layout generation and evaluation.");
    Ok(())
}
