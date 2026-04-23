//! Example: Heterogeneous Multi-Turbine
//!
//! This example demonstrates heterogeneous inflow with multiple turbine types.
//!
//! Corresponds to: examples_heterogeneous/003_het_multi_turbine.rs

use florus::Result;

fn main() -> Result<()> {
    println!("=== Heterogeneous Multi-Turbine ===\n");

    println!("Heterogeneous Inflow with Multiple Turbine Types:\\n");
    
    println!("Scenario:");
    println!("  - Wind farm with different turbine models");
    println!("  - Spatially varying wind conditions");
    println!("  - Complex interactions\\n");

    println!("Challenges:");
    println!("  1. Different rotor diameters:");
    println!("     - Sample wind at different heights");
    println!("     - Varying hub heights");
    println!("     - Different swept areas\\n");

    println!("  2. Different power curves:");
    println!("     - Varying CP/CT characteristics");
    println!("     - Different cut-in/cut-out speeds");
    println!("     - Unique performance profiles\\n");

    println!("  3. Wake interactions:");
    println!("     - Upstream turbines affect downstream");
    println!("     - Different wake characteristics");
    println!("     - Complex superposition\\n");

    println!("Modeling Approach:");
    println!("  - Define heterogeneous wind field");
    println!("  - Assign turbine types to locations");
    println!("  - Sample wind at each turbine hub height");
    println!("  - Calculate individual turbine powers");
    println!("  - Sum for total farm power\\n");

    println!("Applications:");
    println!("  - Farm repowering projects");
    println!("  - Mixed technology deployments");
    println!("  - Phased construction");
    println!("  - Optimization studies\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires heterogeneous inflow and multi-turbine support.");
    Ok(())
}
