//! Example: Optimize yaw for a single wind speed with uncertainty consideration.
//!
//! This example demonstrates yaw optimization across multiple wind directions,
//! comparing baseline and optimized performance.
//!
//! Note: This is a simplified version. Full uncertainty modeling requires
//! UncertainFlorisModel which may not be fully implemented yet.

use florus::{FlorisModel, TimeSeries, Result};
use florus::optimization::yaw_optimization::{YawOptimizationSR, YawOptimizationConfig, YawOptimization};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Yaw Optimization - Single Wind Speed with Uncertainty ===\n");

    // Load the FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Define wind directions to sweep (250° to 290°, step 1°)
    let wind_directions: Vec<f64> = (250..290).map(|d| d as f64).collect();
    let n_dirs = wind_directions.len();
    
    let wind_speeds = vec![8.0; n_dirs];
    let turbulence_intensities = vec![0.06; n_dirs];

    // Set layout: 3 turbines aligned in x-direction, 5D spacing
    let d = 126.0; // Rotor diameter for NREL 5MW
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // Set wind conditions
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions),
        Array1::from_vec(wind_speeds),
        Array1::from_vec(turbulence_intensities),
    )?;
    fmodel.set_wind_data(&time_series)?;

    println!("Wind farm configuration:");
    println!("  - Turbines: 3 (aligned in x-direction)");
    println!("  - Spacing: 5D ({:.0} m) between turbines", 5.0 * d);
    println!("  - Wind directions: {} conditions (250° to 289°)", n_dirs);
    println!("  - Wind speed: 8.0 m/s (constant)");
    println!("  - Turbulence intensity: 0.06 (constant)\n");

    // Initialize optimizer
    println!("Running yaw optimization (Serial-Refine method)...");
    let config = YawOptimizationConfig {
        minimum_yaw_angle: -30.0,
        maximum_yaw_angle: 30.0,
        ..Default::default()
    };
    
    let mut yaw_opt = YawOptimizationSR::new(&fmodel)?;
    let result = yaw_opt.optimize(&mut fmodel, Some(config))?;

    println!("\nOptimization completed!");
    println!("  - Baseline farm power: {:.0} W", result.baseline_farm_power);
    println!("  - Optimized farm power: {:.0} W", result.optimized_farm_power);
    println!("  - Power improvement: {:.0} W", result.power_improvement);
    println!("  - Improvement percentage: {:.2}%", result.improvement_percentage);

    // Show optimal yaw angles for first few wind directions
    println!("\nOptimal yaw angles (first 5 wind directions):");
    for i in 0..5.min(result.optimal_yaw_angles.nrows()) {
        print!("  WD={:.0}°: ", result.wind_directions[i]);
        for ti in 0..result.optimal_yaw_angles.ncols() {
            print!("T{}={:.1}° ", ti, result.optimal_yaw_angles[[i, ti]]);
        }
        println!();
    }

    println!("\n=== Yaw Optimization Complete ===");
    println!("\nNote: For full uncertainty modeling, use UncertainFlorisModel");
    println!("with wind direction standard deviation (wd_std).");

    Ok(())
}
