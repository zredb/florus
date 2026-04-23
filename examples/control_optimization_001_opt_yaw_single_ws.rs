//! Example: Optimize Yaw for Single Wind Speed (from examples_control_optimization folder)
//!
//! Use the serial-refine method to optimize the yaw angles for a 3-turbine wind farm.
//! This corresponds to `examples_control_optimization/001_opt_yaw_single_ws.py` in Python FLORIS.

use florus::{FlorisModel, Result};
use florus::wind_data::TimeSeries;
use florus::optimization::yaw_optimization::{
    YawOptimizationSR, YawOptimizationConfig, YawOptimization
};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Yaw Optimization - Single Wind Speed Example ===\n");
    println!("Optimizing yaw angles using Serial-Refine method for a 3-turbine farm.\n");

    // Load the default example floris object
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Define an inflow that keeps wind speed and TI constant while sweeping wind directions
    let wind_directions: Vec<f64> = (0..120).map(|i| i as f64 * 3.0).collect(); // 0 to 357 deg
    let n = wind_directions.len();
    let wind_speeds = vec![8.0; n];
    let turbulence_intensities = vec![0.06; n];

    let time_series = TimeSeries::new(
        Array1::from(wind_directions),
        Array1::from(wind_speeds),
        Array1::from(turbulence_intensities),
    )?;

    // Reinitialize as a 3-turbine using the above inflow
    let d = 126.0; // Rotor diameter for the NREL 5 MW
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 5.0 * d, 10.0 * d]),
        &ndarray::arr1(&[0.0, 0.0, 0.0]),
    )?;
    fmodel.set_wind_data(&time_series)?;

    println!("Wind farm configuration:");
    println!("  - Turbines: 3 (aligned in x-direction)");
    println!("  - Spacing: 5D (630 m) between turbines");
    println!("  - Wind directions: {} conditions (0° to 357°, step 3°)", n);
    println!("  - Wind speed: 8.0 m/s (constant)");
    println!("  - Turbulence intensity: 0.06 (constant)\n");

    // Run the model before optimization
    println!("Running baseline simulation...");
    fmodel.run()?;
    println!("Baseline simulation completed.\n");

    // Initialize optimizer object and run optimization using the Serial-Refine method
    println!("Running yaw optimization (Serial-Refine method)...");
    let mut yaw_opt = YawOptimizationSR::new();
    
    let config = YawOptimizationConfig {
        minimum_yaw_angle: -30.0,
        maximum_yaw_angle: 30.0,
        ..Default::default()
    };

    let result = yaw_opt.optimize(&mut fmodel, Some(config))?;

    println!("\nOptimization completed!");
    println!("  - Baseline farm power: {:.0} W", result.baseline_power);
    println!("  - Optimized farm power: {:.0} W", result.optimized_power);
    println!("  - Power improvement: {:.0} W", result.power_improvement);
    println!("  - Improvement percentage: {:.2}%", result.improvement_percentage);

    // Print yaw angles for each turbine
    println!("\nOptimal yaw angles (first wind direction):");
    if result.yaw_angles.dim().0 > 0 {
        for (t, &yaw) in result.yaw_angles.row(0).iter().enumerate() {
            println!("  Turbine {}: {:.1}°", t, yaw);
        }
    }

    println!("\n=== Yaw Optimization Complete ===");
    println!("\nKey observations:");
    println!("  - Upstream turbine (T0) typically yaws to deflect wake");
    println!("  - Downstream turbines benefit from wake deflection");
    println!("  - Power uplift varies with wind direction");
    println!("  - Maximum uplift occurs when turbines are closely aligned with wind");

    Ok(())
}
