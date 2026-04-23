//! Example: Generate Value Table
//!
//! This example demonstrates how to generate and work with wind condition value tables,
//! showing the relationship between wind speed, direction, and power output.
//!
//! Corresponds to: examples_wind_data/003_generate_value.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Generate Wind Condition Value Table ===\n");

    // Initialize FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set a simple 2-turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    // Create a grid of wind conditions
    let wind_speeds = vec![6.0, 8.0, 10.0, 12.0, 14.0];
    let wind_directions = vec![240.0, 270.0, 300.0];
    let ti = 0.06;

    println!("Power Output Matrix (kW)");
    println!("Wind Direction vs Wind Speed\n");

    // Header
    print!("{:<12} |", "WD \\ WS");
    for ws in &wind_speeds {
        print!(" {:>10.1}", ws);
    }
    println!();
    print!("{}", "-".repeat(13));
    for _ in &wind_speeds {
        print!("|{}", "-".repeat(11));
    }
    println!();

    // Data rows
    for wd in &wind_directions {
        print!("{:>10.0}° |", wd);
        
        for ws in &wind_speeds {
            fmodel.set_wind_conditions(
                ndarray::arr1(&[*ws]),
                ndarray::arr1(&[*wd]),
                ndarray::arr1(&[ti]),
            )?;
            fmodel.run()?;

            let powers = fmodel.get_turbine_powers();
            let farm_power: f64 = powers.row(0).sum() / 1000.0;
            print!(" {:>10.1}", farm_power);
        }
        println!();
    }

    println!("\n=== Analysis ===");
    println!("This table shows:");
    println!("  - Power varies with both wind speed and direction");
    println!("  - Wake effects reduce power at certain directions");
    println!("  - Higher wind speeds generally produce more power");
    println!("  - Direction affects wake alignment between turbines");

    println!("\n=== Example Complete ===");
    println!("\nNote: Full value table generation would create comprehensive");
    println!("lookup tables for AEP calculations and optimization.");

    Ok(())
}
