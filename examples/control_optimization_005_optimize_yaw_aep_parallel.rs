//! Example: Parallel Yaw Optimization for AEP
//!
//! This example demonstrates parallel computation for yaw optimization.
//!
//! Corresponds to: examples_control_optimization/005_optimize_yaw_aep_parallel.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Parallel Yaw Optimization - AEP ===\n");

    println!("Parallel Optimization Concept:\\n");
    println!("AEP optimization is computationally expensive:");
    println!("  - Multiple wind directions (e.g., 36)");
    println!("  - Multiple wind speeds (e.g., 10)");
    println!("  - Total: 360+ conditions to optimize\\n");

    println!("Parallelization Strategies:");
    println!("  1. Wind condition parallelism:");
    println!("     - Each core optimizes different conditions");
    println!("     - Embarrassingly parallel");
    println!("     - Linear speedup\\n");

    println!("  2. Optimization algorithm parallelism:");
    println!("     - Parallel gradient evaluation");
    println!("     - Multi-start optimization");
    println!("     - Population-based methods\\n");

    println!("Benefits:");
    println!("  - Reduce wall-clock time");
    println!("  - Enable fine-grained wind roses");
    println!("  - Make AEP optimization practical\\n");

    println!("Implementation Notes:");
    println!("  - Use Rayon for shared-memory parallelism");
    println!("  - Consider distributed computing for large cases");
    println!("  - Balance load across cores\\n");

    println!("=== Example Complete ===");
    println!("Note: Full implementation requires parallel optimization framework.");
    Ok(())
}
