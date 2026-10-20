//! Example: Compare yaw optimizers
//!
//! This example performs yaw optimizations with three different optimizers
//! (Serial Refine, Geometric and SciPy) over a full range of wind
//! directions at a single wind speed, and compares their results.
//!
//! Scenario aligned 1:1 with
//! examples_control_optimization/006_compare_yaw_optimizers.py.

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_geometric::YawOptimizationGeometric;
use florus::optimization::yaw_optimization::yaw_optimizer_scipy::YawOptimizationScipy;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::YawOptimization;
use florus::types::{Array1, Float};
use ndarray::Array;
use std::time::Instant;

fn main() -> anyhow::Result<()> {
    // Load the default example floris object
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Reinitialize as a 3-turbine farm with a range of WDs and one WS
    let d = 126.0; // Rotor diameter for the NREL 5 MW
    let wd_array = Array::range(0.0, 360.0, 3.0);
    let ws_array = Array::from_elem(wd_array.len(), 8.0);
    let turbulence_intensities = Array::from_elem(wd_array.len(), 0.06);
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;
    fmodel.set_wind_conditions(ws_array, wd_array, turbulence_intensities)?;
    fmodel.run()?;

    println!("Performing optimizations with SciPy...");
    let start = Instant::now();
    let mut yaw_opt_scipy = YawOptimizationScipy::new();
    let df_opt_scipy = yaw_opt_scipy.optimize(&mut fmodel, None)?;
    let time_scipy = start.elapsed().as_secs_f64();

    println!("Performing optimizations with Serial Refine...");
    let start = Instant::now();
    let mut yaw_opt_sr = YawOptimizationSR::new();
    let df_opt_sr = yaw_opt_sr.optimize(&mut fmodel, None)?;
    let time_sr = start.elapsed().as_secs_f64();

    println!("Performing optimizations with Geometric Yaw...");
    let start = Instant::now();
    let mut yaw_opt_geo = YawOptimizationGeometric::new();
    let df_opt_geo = yaw_opt_geo.optimize(&mut fmodel, None)?;
    let time_geo = start.elapsed().as_secs_f64();

    // Print time spent
    println!("\n Time spent, Geometric Yaw: {time_geo:.2} s.");
    println!(" Time spent, Serial Refine: {time_sr:.2} s.");
    println!(" Time spent, SciPy (SLSQP): {time_scipy:.2} s.\n");

    // The geometric optimizer does not compute power within the optimization
    // (Python floris_model.get_farm_power() after setting the geometric yaw)
    let geo_farm_power = df_opt_geo.farm_powers_opt.clone();

    // Finally, compare the overall power gains
    let total_power_uplift_geo: Float = (0..geo_farm_power.len())
        .map(|i| geo_farm_power[i] - df_opt_sr.farm_powers_baseline[i])
        .sum();
    let total_power_uplift_sr: Float = (0..df_opt_sr.farm_powers_opt.len())
        .map(|i| df_opt_sr.farm_powers_opt[i] - df_opt_sr.farm_powers_baseline[i])
        .sum();
    let total_power_uplift_scipy: Float = (0..df_opt_scipy.farm_powers_opt.len())
        .map(|i| df_opt_scipy.farm_powers_opt[i] - df_opt_scipy.farm_powers_baseline[i])
        .sum();

    println!("Total power uplift, Geometric Yaw: {total_power_uplift_geo:.3e} W");
    println!("Total power uplift, Serial Refine: {total_power_uplift_sr:.3e} W");
    println!("Total power uplift, SciPy (SLSQP): {total_power_uplift_scipy:.3e} W");

    Ok(())
}
