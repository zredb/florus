//! Example: Control Types Overview
//!
//! This example demonstrates the different control types available in FLORUS:
//! yaw angles, tilt angles, and power setpoints.
//!
//! Corresponds to: examples_control_types/001_control_types_overview.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Control Types Overview ===\n");

    // Initialize FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set a simple 2-turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    // Set wind conditions
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    println!("Available Control Types:\n");
    println!("1. Yaw Angles (degrees)");
    println!("   - Controls turbine orientation relative to wind");
    println!("   - Used for wake steering");
    println!("   - Range: typically -30° to +30°\n");

    println!("2. Tilt Angles (degrees)");
    println!("   - Controls rotor tilt from horizontal");
    println!("   - Used for vertical wake deflection");
    println!("   - Important for floating turbines\n");

    println!("3. Power Setpoints (Watts)");
    println!("   - Maximum power output limit");
    println!("   - Used for curtailment");
    println!("   - Can reduce power below available\n");

    println!("4. Active Wake Mixing Control");
    println!("   - Induces additional mixing in wake");
    println!("   - Accelerates wake recovery");
    println!("   - Requires specific turbine capabilities\n");

    // Demonstrate setting yaw angles
    println!("--- Demonstrating Yaw Control ---\n");
    fmodel.set_yaw_angles(ndarray::arr2(&[[10.0, 0.0]]))?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers();
    println!("With T0 yawed 10°:");
    println!("  T0: {:.2} kW", powers[[0, 0]] / 1000.0);
    println!("  T1: {:.2} kW", powers[[0, 1]] / 1000.0);
    println!("  Total: {:.2} kW\n", (powers[[0, 0]] + powers[[0, 1]]) / 1000.0);

    println!("=== Analysis ===");
    println!("Control types enable:");
    println!("  - Wake steering (yaw)");
    println!("  - Vertical wake control (tilt)");
    println!("  - Power curtailment (setpoints)");
    println!("  - Advanced wake mixing (AWC)");

    println!("\n=== Example Complete ===");

    Ok(())
}
