//! Example: WRG with Heterogeneous Inflow
//!
//! This example demonstrates using WRG data for heterogeneous inflow.
//!
//! Corresponds to: examples_wind_resource_grid/003_wrg_with_heterogeneous.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== WRG with Heterogeneous Inflow ===\n");

    println!("WRG-Driven Heterogeneous Inflow:\\n");
    
    println!("Integration Concept:");
    println!("  - Use WRG grid as heterogeneous wind field");
    println!("  - Interpolate WRG data to turbine locations");
    println!("  - Apply spatially varying wind conditions\\n");

    println!("Process:");
    println!("  1. Load WRG file");
    println!("  2. Extract wind resource grid");
    println!("  3. For each turbine:");
    println!("     - Find nearest WRG grid points");
    println!("     - Interpolate wind speed/direction");
    println!("     - Apply TI and shear");
    println!("  4. Run FLORUS simulation\\n");

    println!("Benefits:");
    println!("  - Realistic wind variation");
    println!("  - Industry-standard data source");
    println!("  - Comprehensive wind characterization");
    println!("  - Improved AEP accuracy\\n");

    println!("Challenges:");
    println!("  - Coordinate system conversion (lat/lon to x/y)");
    println!("  - Height interpolation");
    println!("  - Data resolution vs. farm size");
    println!("  - Computational overhead\\n");

    println!("Applications:");
    println!("  - Final energy yield assessment");
    println!("  - Bankable P50/P90 calculations");
    println!("  - Complex terrain sites");
    println!("  - Large offshore farms\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires WRG parser and heterogeneous inflow integration.");
    Ok(())
}
