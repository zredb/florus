//! Example: Compare Yaw Optimizers
//!
//! This example compares different yaw optimization algorithms.
//!
//! Corresponds to: examples_control_optimization/006_compare_yaw_optimizers.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Compare Yaw Optimizers ===\n");

    println!("Optimization Algorithm Comparison:\\n");
    
    println!("1. Serial Refine (SR):");
    println!("   - Sequential turbine optimization");
    println!("   - Fast convergence");
    println!("   - Good for real-time control\\n");

    println!("2. Genetic Algorithm (GA):");
    println!("   - Population-based search");
    println!("   - Global optimization");
    println!("   - Slower but more robust\\n");

    println!("3. Gradient-based methods:");
    println!("   - Use analytical gradients");
    println!("   - Very fast convergence");
    println!("   - May get stuck in local optima\\n");

    println!("4. Random Search:");
    println!("   - Simple baseline");
    println!("   - No gradient needed");
    println!("   - Slow convergence\\n");

    println!("Comparison Metrics:");
    println!("  - Solution quality (power gain)");
    println!("  - Computation time");
    println!("  - Robustness to initial conditions");
    println!("  - Scalability with turbine count\\n");

    println!("Recommendations:");
    println!("  - Real-time control: Serial Refine");
    println!("  - Offline AEP: Genetic Algorithm");
    println!("  - Research: Compare multiple methods\\n");

    println!("=== Example Complete ===");
    println!("Note: Full comparison requires multiple optimizer implementations.");
    Ok(())
}
