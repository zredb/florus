//! Example: Floating Turbine Models
//!
//! This example demonstrates floating wind turbine modeling in FLORUS.
//!
//! Corresponds to: examples_floating/001_floating_turbine_models.rs

use florus::Result;

fn main() -> Result<()> {
    println!("=== Floating Turbine Models ===\n");

    println!("Floating Wind Turbine Modeling:\\n");
    
    println!("Key Differences from Fixed-Bottom:\\n");
    
    println!("1. Platform Motion:");
    println!("   - Pitch, roll, heave motions");
    println!("   - Affects rotor orientation");
    println!("   - Changes effective wind speed\\n");

    println!("2. Tilt Control:");
    println!("   - Platform pitch changes tilt angle");
    println!("   - Vertical wake deflection");
    println!("   - Can be used for wake steering\\n");

    println!("3. Multi-Dimensional CP/CT:");
    println!("   - Depends on wave conditions");
    println!("   - Wave height (Hs) effects");
    println!("   - Wave period (Tp) effects\\n");

    println!("4. Dynamic Response:");
    println!("   - Time-varying platform motion");
    println!("   - Coupled aero-hydro-servo-elastic");
    println!("   - Requires specialized models\\n");

    println!("Floating Platform Types:");
    println!("  - Semi-submersible");
    println!("  - Spar buoy");
    println!("  - Tension-leg platform (TLP)");
    println!("  - Barge\\n");

    println!("Applications:");
    println!("  - Deep water sites (>60m depth)");
    println!("  - Offshore wind farms");
    println!("  - Higher capacity factors");
    println!("  - Reduced visual impact\\n");

    println!("=== Example Complete ===");
    println!("Note: Full floating support requires platform motion models.");
    Ok(())
}
