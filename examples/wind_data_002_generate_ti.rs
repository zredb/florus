//! Example: Generate Turbulence Intensity (TI) Table
//!
//! This example demonstrates how to work with turbulence intensity data
//! and shows typical TI values for different wind conditions.
//!
//! Corresponds to: examples_wind_data/002_generate_ti.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Generate Turbulence Intensity Table ===\n");

    // Initialize FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set a simple 2-turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    // Demonstrate different TI levels
    let ti_levels = vec![0.04, 0.06, 0.08, 0.10, 0.12];
    let ws = 8.0;
    let wd = 270.0;

    println!("Turbulence Intensity Effects on Power Output");
    println!("Wind Speed: {:.1} m/s, Wind Direction: {:.0}°\n", ws, wd);
    println!("{:<15} | {:<15} | {:<15} | {:<15}", "TI", "T0 Power (kW)", "T1 Power (kW)", "Farm Total (kW)");
    println!("{}|{}|{}|{}", "-".repeat(16), "-".repeat(16), "-".repeat(16), "-".repeat(16));

    for ti in &ti_levels {
        fmodel.set_wind_conditions(
            ndarray::arr1(&[ws]),
            ndarray::arr1(&[wd]),
            ndarray::arr1(&[*ti]),
        )?;
        fmodel.run()?;

        let powers = fmodel.get_turbine_powers();
        let t0_power = powers[[0, 0]] / 1000.0;
        let t1_power = powers[[0, 1]] / 1000.0;
        let farm_total = (t0_power + t1_power);

        println!("{:<15.2} | {:<15.2} | {:<15.2} | {:<15.2}", ti, t0_power, t1_power, farm_total);
    }

    println!("\n=== Analysis ===");
    println!("Turbulence intensity affects:");
    println!("  - Wake recovery rate (higher TI = faster recovery)");
    println!("  - Downstream turbine power (higher TI = more power)");
    println!("  - Overall farm efficiency");
    println!("\nTypical TI values:");
    println!("  - Offshore: 0.04 - 0.06 (low turbulence)");
    println!("  - Onshore flat terrain: 0.06 - 0.08");
    println!("  - Complex terrain: 0.08 - 0.12+ (high turbulence)");

    println!("\n=== Example Complete ===");
    println!("\nNote: Full TI table generation would create lookup tables");
    println!("mapping wind speed/direction to expected turbulence intensity.");

    Ok(())
}
