//! Example: Tilt-Driven Vertical Wake Deflection
//!
//! This example demonstrates using tilt for vertical wake steering.
//!
//! Corresponds to: examples_floating/003_tilt_driven_vertical_wake_deflection.rs

use florus::Result;

fn main() -> Result<()> {
    println!("=== Tilt-Driven Vertical Wake Deflection ===\n");

    println!("Vertical Wake Steering Concept:\\n");
    
    println!("Traditional Yaw Control:");
    println!("  - Horizontal wake deflection");
    println!("  - Redirects wake left or right");
    println!("  - Well-established technique\\n");

    println!("Tilt Control:");
    println!("  - Vertical wake deflection");
    println!("  - Redirects wake up or down");
    println!("  - Emerging control strategy\\n");

    println!("Mechanism:");
    println!("  1. Change rotor tilt angle");
    println!("  2. Alters wake trajectory vertically");
    println!("  3. Wake moves above/below downstream rotors");
    println!("  4. Reduces wake losses\\n");

    println!("Floating Turbine Advantage:");
    println!("  - Platform pitch naturally changes tilt");
    println!("  - Can be actively controlled");
    println!("  - No additional actuators needed");
    println!("  - Dual benefit: stability + control\\n");

    println!("Benefits:");
    println!("  - Complements yaw control");
    println!("  - 3D wake manipulation");
    println!("  - Higher optimization potential");
    println!("  - Reduced structural loads\\n");

    println!("Challenges:");
    println!("  - Slower response than yaw");
    println!("  - Platform dynamics coupling");
    println!("  - Control complexity");
    println!("  - Model accuracy needed\\n");

    println!("Applications:");
    println!("  - Aligned turbine rows");
    println!("  - Floating offshore farms");
    println!("  - Complex terrain");
    println!("  - Multi-dimensional optimization\\n");

    println!("=== Example Complete ===");
    println!("Note: Full implementation requires tilt control API.");
    Ok(())
}
