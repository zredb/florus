//! Example: Multiple Turbine Types
//!
//! This example uses an input file where multiple turbine types are defined.
//! The first two turbines are the NREL 5MW, and the third turbine is the IEA 10MW.
//!
//! Corresponds to: examples_turbine/002_multiple_turbine_types.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Multiple Turbine Types ===\n");

    // Initialize FLORIS with the given input file
    let fmodel = FlorisModel::from_file("examples/inputs/gch_multiple_turbine_types.yaml")?;

    println!("Farm Configuration:");
    println!("  Number of turbines: {}", fmodel.core().farm.turbines.len());
    println!("  Layout X: {:?}", fmodel.core().farm.layout_x);
    println!("  Layout Y: {:?}", fmodel.core().farm.layout_y);
    println!("  Rotor diameters: {:?}", fmodel.core().farm.rotor_diameters);

    // Run simulation
    println!("\nRunning simulation...\n");
    let mut fmodel = fmodel;
    fmodel.run()?;

    // Get turbine powers
    let powers = fmodel.get_turbine_powers();
    
    println!("Turbine Powers:");
    for ti in 0..powers.ncols() {
        let power_kw = powers[[0, ti]] / 1000.0;
        println!("  T{}: {:.2} kW", ti, power_kw);
    }
    
    let farm_power: f64 = powers.row(0).sum();
    println!("\nFarm Total Power: {:.2} kW", farm_power / 1000.0);

    println!("\n=== Analysis ===");
    println!("This example demonstrates:");
    println!("  - Multiple turbine types in a single farm");
    println!("  - Different turbine characteristics (rotor diameter, hub height)");
    println!("  - Wake interactions between different turbine types");

    println!("\n=== Example Complete ===");
    println!("\nNote: Full visualization would show horizontal, y-plane, and cross-plane cuts.");

    Ok(())
}
