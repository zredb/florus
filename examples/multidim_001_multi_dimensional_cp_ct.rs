//! Example: Multi-Dimensional CP/CT Tables
//!
//! This example demonstrates multi-dimensional power/thrust coefficient tables.
//!
//! Corresponds to: examples_multidim/001_multi_dimensional_cp_ct.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Multi-Dimensional CP/CT Tables ===\n");

    println!("Multi-Dimensional Look-up Tables:\\n");
    println!("Standard turbines use 1D tables (wind speed only).");
    println!("Multi-dimensional tables add dependencies on:");
    println!("  - Turbulence intensity (TI)");
    println!("  - Wave height (Hs) for floating turbines");
    println!("  - Wave period (Tp) for floating turbines");

    println!("\nBenefits:");
    println!("  - More accurate power predictions");
    println!("  - Captures environmental effects");
    println!("  - Essential for floating offshore wind");

    println!("\nTable Dimensions:");
    println!("  - 2D: Wind speed × TI");
    println!("  - 3D: Wind speed × Hs × Tp");
    println!("  - 4D+: Multiple environmental parameters");

    println!("\n=== Example Complete ===");
    println!("Note: Full multi-dim support requires turbine files with multi-dim tables.");
    Ok(())
}
