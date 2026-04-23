//! Example: Operation Model Overview
//!
//! This example demonstrates the operation model settings in FLORUS,
//! showing how different operational parameters affect turbine behavior.
//!
//! Corresponds to: examples_operation_models/001_operation_model.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Operation Model Overview ===\n");

    // Initialize FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set a simple 2-turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    // Set wind conditions
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0, 10.0, 12.0]),
        ndarray::arr1(&[270.0, 270.0, 270.0]),
        ndarray::arr1(&[0.06, 0.06, 0.06]),
    )?;

    println!("Farm Configuration:");
    println!("  Turbines: 2 (NREL 5MW)");
    println!("  Spacing: 500m (3.97D)");
    println!("  Wind Direction: 270°");
    println!("  Wind Speeds: 8.0, 10.0, 12.0 m/s\n");

    // Test with different yaw angles
    println!("Yaw Angle Effects:\n");
    println!("{:<15} | {:<15} | {:<15} | {:<15}", "Condition", "T0 Power (kW)", "T1 Power (kW)", "Farm Total (kW)");
    println!("{}|{}|{}|{}", "-".repeat(16), "-".repeat(16), "-".repeat(16), "-".repeat(16));

    let yaw_configs = vec![
        ("No Yaw", [0.0, 0.0]),
        ("T0 Yaw 10°", [10.0, 0.0]),
        ("T0 Yaw 20°", [20.0, 0.0]),
        ("T0 Yaw 30°", [30.0, 0.0]),
    ];

    for (config_name, yaw_angles) in &yaw_configs {
        fmodel.set_yaw_angles(ndarray::arr2(&[
            *yaw_angles,
            *yaw_angles,
            *yaw_angles,
        ]))?;
        fmodel.run()?;

        let powers = fmodel.get_turbine_powers();
        
        for i in 0..3 {
            let t0_power = powers[[i, 0]] / 1000.0;
            let t1_power = powers[[i, 1]] / 1000.0;
            let farm_total = t0_power + t1_power;
            
            if i == 0 {
                print!("{:<15} |", config_name);
            } else {
                print!("{:<15} |", "");
            }
            println!(" {:>13.2} | {:>13.2} | {:>13.2}", t0_power, t1_power, farm_total);
        }
        println!();
    }

    println!("=== Analysis ===");
    println!("Operation models control:");
    println!("  - Yaw angles (wake steering)");
    println!("  - Tilt angles (vertical wake deflection)");
    println!("  - Power setpoints (curtailment)");
    println!("  - Active wake mixing control");
    println!("\nKey observations:");
    println!("  - Upstream turbine yaw reduces its power");
    println!("  - Downstream turbine gains from wake deflection");
    println!("  - Optimal yaw balances total farm power");

    println!("\n=== Example Complete ===");

    Ok(())
}
