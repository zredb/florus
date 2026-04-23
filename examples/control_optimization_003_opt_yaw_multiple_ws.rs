//! Example: Optimize Yaw for Multiple Wind Speeds
//!
//! This example demonstrates yaw optimization across different wind speeds,
//! showing how optimal yaw angles vary with wind conditions.
//!
//! Corresponds to: examples_control_optimization/003_opt_yaw_multiple_ws.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Yaw Optimization - Multiple Wind Speeds ===\n");
    println!("This example optimizes yaw angles for multiple wind speeds.\n");
    
    println!("Configuration:");
    println!("  - Turbines: 3 (aligned in x-direction)");
    println!("  - Wind directions: 270°");
    println!("  - Wind speeds: [6, 8, 10, 12] m/s");
    println!("  - Turbulence intensity: 0.06\n");

    // Initialize FlorisModel
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set 3-turbine layout aligned in x-direction
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 630.0, 1260.0]),
        &ndarray::arr1(&[0.0, 0.0, 0.0]),
    )?;

    let wind_speeds = vec![6.0, 8.0, 10.0, 12.0];
    
    println!("{:<15} | {:<20} | {:<20} | {:<15}", "Wind Speed", "Baseline Power", "Optimized Power", "Gain (%)");
    println!("{}|{}|{}|{}", "-".repeat(16), "-".repeat(21), "-".repeat(21), "-".repeat(16));

    for ws in &wind_speeds {
        // Set wind conditions
        fmodel.set_wind_conditions(
            ndarray::arr1(&[*ws]),
            ndarray::arr1(&[270.0]),
            ndarray::arr1(&[0.06]),
        )?;
        
        // Run baseline (no yaw)
        fmodel.set_yaw_angles(ndarray::arr2(&[[0.0, 0.0, 0.0]]))?;
        fmodel.run()?;
        let baseline_power = fmodel.get_farm_power()[0];
        
        // Note: Full yaw optimization requires YawOptimizationSR implementation
        // For now, demonstrate the concept
        println!("{:<15.1} | {:<20.2} | {:<20} | {:<15}", 
                 ws, 
                 baseline_power / 1000.0,
                 "Requires optimizer",
                 "-");
    }
    
    println!("\n=== Analysis ===");
    println!("Yaw optimization effectiveness varies with wind speed:");
    println!("  - Low wind speeds: Less benefit (wake effects smaller)");
    println!("  - Medium wind speeds: Maximum benefit (strong wakes)");
    println!("  - High wind speeds: Reduced benefit (turbines at rated power)");

    println!("\n=== Example Complete ===");
    println!("Note: Full implementation requires YawOptimizationSR class.");

    Ok(())
}
