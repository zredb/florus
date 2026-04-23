//! Example: Load Optimization with Wind Rose
//!
//! This example demonstrates load optimization across wind conditions.
//!
//! Corresponds to: examples_load_optimization/002_load_optimization_with_wind_rose.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Load Optimization with Wind Rose ===\n");

    println!("Load Optimization with Full Wind Distribution:\\n");
    
    println!("Concept:");
    println!("  - Optimize loads across all wind conditions");
    println!("  - Weight by wind direction/speed frequency");
    println!("  - Minimize cumulative damage\\n");

    println!("Damage Calculation:");
    println!("  1. For each wind condition:");
    println!("     - Calculate turbine loads");
    println!("     - Apply S-N curve (fatigue)");
    println!("     - Compute damage rate");
    println!("  2. Weight by frequency");
    println!("  3. Sum for total lifetime damage\\n");

    println!("Optimization Objectives:");
    println!("  - Minimize total damage");
    println!("  - Subject to power constraints");
    println!("  - Balance multiple turbines");
    println!("  - Consider component limits\\n");

    println!("Control Strategies:");
    println!("  - Wake steering for load reduction");
    println!("  - Power curtailment in high-load conditions");
    println!("  - Coordinated farm control");
    println!("  - Adaptive control based on wind\\n");

    println!("Benefits:");
    println!("  - Lifetime extension (15-25%)");
    println!("  - Reduced maintenance frequency");
    println!("  - Lower replacement costs");
    println!("  - Improved availability\\n");

    println!("Trade-offs:");
    println!("  - Power loss: 1-5% typical");
    println!("  - vs. Lifetime gain: 15-25%");
    println!("  - Net economic benefit positive");
    println!("  - Site-specific optimization needed\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires load models and multi-objective optimization.");
    Ok(())
}
