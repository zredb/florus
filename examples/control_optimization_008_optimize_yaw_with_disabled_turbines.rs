//! Example: Optimize yaw with disabled turbines.
//!
//! This example demonstrates yaw optimization when some turbines are disabled,
//! showing how the optimizer adapts to different farm configurations.

use florus::{FlorisModel, TimeSeries, Result};
use florus::optimization::yaw_optimization::{YawOptimizationSR, YawOptimizationConfig, YawOptimization};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Yaw Optimization - With Disabled Turbines ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set layout: 3 turbines aligned
    let d = 126.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // Set wind conditions
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    fmodel.set_wind_data(&time_series)?;

    println!("Configuration:");
    println!("  - Turbines: 3 (aligned in x-direction)");
    println!("  - Wind direction: 270°");
    println!("  - Wind speed: 8.0 m/s\n");

    // Test 1: All turbines enabled
    println!("Test 1: All turbines enabled");
    fmodel.disable_turbines(&[])?; // Enable all
    
    let config = YawOptimizationConfig {
        minimum_yaw_angle: -30.0,
        maximum_yaw_angle: 30.0,
        ..Default::default()
    };
    
    let mut yaw_opt = YawOptimizationSR::new(&fmodel)?;
    let result = yaw_opt.optimize(&mut fmodel, Some(config))?;
    
    println!("  Baseline power: {:.0} W", result.baseline_farm_power);
    println!("  Optimized power: {:.0} W", result.optimized_farm_power);
    println!("  Improvement: {:.2}%\n", result.improvement_percentage);

    // Test 2: Middle turbine disabled
    println!("Test 2: Middle turbine (T1) disabled");
    fmodel.disable_turbines(&[1])?;
    
    let mut yaw_opt = YawOptimizationSR::new(&fmodel)?;
    let result = yaw_opt.optimize(&mut fmodel, Some(config))?;
    
    println!("  Baseline power: {:.0} W", result.baseline_farm_power);
    println!("  Optimized power: {:.0} W", result.optimized_farm_power);
    println!("  Improvement: {:.2}%\n", result.improvement_percentage);

    println!("=== Optimization Complete ===");
    println!("\nKey observations:");
    println!("  - Disabling upstream turbine changes wake patterns");
    println!("  - Optimizer adapts yaw angles accordingly");
    println!("  - Power improvement varies with configuration");

    Ok(())
}
