//! Example: Peak-Shaving Operation Model
//!
//! Demonstrates peak-shaving thrust reduction near rated wind speed.
//!
//! Corresponds to: examples_control_types/005_peak_shaving.py

use florus::{FlorisModel, Result};
use florus::wind_data::TimeSeries;
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Peak-Shaving Operation Model Example ===\n");

    let n_ws = 100;
    let wind_speeds: Vec<f64> = (0..n_ws).map(|i| i as f64 * 30.0 / (n_ws - 1) as f64).collect();

    // Simple model (baseline)
    let mut fmodel_simple = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_simple.set_layout(
        &ndarray::arr1(&[0.0, 1000.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;
    let ts = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_ws]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.10; n_ws]),
    )?;
    fmodel_simple.set_wind_data(&ts)?;
    fmodel_simple.set_operation_model("simple")?;
    fmodel_simple.run()?;
    let powers_simple = fmodel_simple.get_turbine_powers() / 1000.0;
    let ct_simple = fmodel_simple.get_turbine_thrust_coefficients();

    // Peak-shaving model
    let mut fmodel_ps = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_ps.set_layout(
        &ndarray::arr1(&[0.0, 1000.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;
    let ts2 = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_ws]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.10; n_ws]),
    )?;
    fmodel_ps.set_wind_data(&ts2)?;
    fmodel_ps.set_operation_model("peak-shaving")?;
    fmodel_ps.run()?;
    let powers_ps = fmodel_ps.get_turbine_powers() / 1000.0;
    let ct_ps = fmodel_ps.get_turbine_thrust_coefficients();

    println!("Upstream turbine (T0) comparison:");
    println!("{:<12} | {:<15} | {:<15} | {:<15} | {:<15}", "WS (m/s)", "Simple P (kW)", "PeakShave P", "Simple Ct", "PeakShave Ct");
    println!("{}", "-".repeat(78));
    for i in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90] {
        println!(
            "{:<12.1} | {:<15.2} | {:<15.2} | {:<15.4} | {:<15.4}",
            wind_speeds[i], powers_simple[[i, 0]], powers_ps[[i, 0]], ct_simple[[i, 0]], ct_ps[[i, 0]]
        );
    }

    // Show where peak-shaving reduces thrust
    let mut diff_count = 0;
    for i in 0..n_ws {
        if (ct_simple[[i, 0]] - ct_ps[[i, 0]]).abs() > 0.001 {
            diff_count += 1;
        }
    }
    println!("\nPeak-shaving active at {}/{} wind speed points", diff_count, n_ws);

    println!("\n=== Example Complete ===");
    Ok(())
}
