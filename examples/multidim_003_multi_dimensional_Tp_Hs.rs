//! Example: Multi-Dimensional Tp-Hs Tables
//!
//! This example demonstrates wave period and height dependencies.
//!
//! Corresponds to: examples_multidim/003_multi_dimensional_Tp_Hs.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Multi-Dimensional Tp-Hs Tables ===\n");

    println!("Wave Period and Height Dependencies:\\n");
    println!("Floating turbines are affected by:");
    println!("  - Hs (Significant wave height)");
    println!("  - Tp (Wave peak period)");

    println!("\nEffects on Performance:");
    println!("  - Platform pitch/roll motion");
    println!("  - Rotor effective wind speed");
    println!("  - Power and thrust variations");

    println!("\nTable Structure:");
    println!("  CP[wind_speed, Hs, Tp]");
    println!("  CT[wind_speed, Hs, Tp]");

    println!("\n=== Example Complete ===");
    Ok(())
}
