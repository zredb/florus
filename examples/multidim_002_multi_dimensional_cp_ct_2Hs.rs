//! Example: Multi-Dimensional CP/CT with 2Hs
//!
//! This example demonstrates dual significant wave height dimensions.
//!
//! Corresponds to: examples_multidim/002_multi_dimensional_cp_ct_2Hs.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Multi-Dimensional CP/CT with 2Hs ===\n");

    println!("Dual Wave Height Dimensions:\\n");
    println!("Some floating turbine models use two Hs values:");
    println!("  - Hs_current: Current wave conditions");
    println!("  - Hs_design: Design wave conditions");

    println!("\nPurpose:");
    println!("  - Capture hysteresis effects");
    println!("  - Model structural response");
    println!("  - Account for platform motion history");

    println!("\nApplications:");
    println!("  - Semi-submersible platforms");
    println!("  - Spar buoy designs");
    println!("  - Complex floating systems");

    println!("\n=== Example Complete ===");
    Ok(())
}
