//! Example: Simple-Derating Control
//!
//! Demonstrates the simple-derating operation model in FLORUS.
//! The simple-derating model allows specifying a power setpoint for each turbine,
//! limiting power output to at most the setpoint.
//!
//! Corresponds to: examples_control_types/001_derating_control.py

use florus::{FlorisModel, Result};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Simple-Derating Control Example ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 2-turbine layout first (must be before set_operation_model)
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 1000.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    // Change to the simple-derating operation model (AFTER set_layout)
    fmodel.set_operation_model("simple-derating")?;

    // Set constant wind conditions over N=50 time steps
    let n = 50;
    let wind_speeds = Array1::from_vec(vec![10.0; n]);
    let wind_directions = Array1::from_vec(vec![270.0; n]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; n]);
    fmodel.set_wind_conditions(wind_speeds, wind_directions, turbulence_intensities)?;

    fmodel.run()?;
    let turbine_powers_orig = fmodel.get_turbine_powers();

    println!("Baseline (no derating):");
    println!("  Upstream (T0):   {:.2} kW", turbine_powers_orig[[0, 0]] / 1000.0);
    println!("  Downstream (T1): {:.2} kW", turbine_powers_orig[[0, 1]] / 1000.0);

    // Create power setpoints sweeping from 1 kW to 6000 kW for both turbines
    let setpoints: Vec<f64> = (0..n)
        .map(|i| 1000.0 + (i as f64) * (5999000.0 / (n - 1) as f64))
        .collect();

    // Both turbines get the same setpoint sweep
    let mut power_setpoints = ndarray::Array2::zeros((n, 2));
    for i in 0..n {
        power_setpoints[[i, 0]] = setpoints[i];
        power_setpoints[[i, 1]] = setpoints[i];
    }

    fmodel.set_operation(
        None,
        Some(power_setpoints.clone()),
        None, None, None, None,
    )?;
    fmodel.run()?;
    let turbine_powers_derated = fmodel.get_turbine_powers();

    println!("\nWith derating (power setpoints applied to both turbines):");
    println!("{:<18} | {:<15} | {:<15} | {:<15}", "Setpoint (kW)", "T0 Power (kW)", "T1 Power (kW)", "Farm (kW)");
    println!("{}", "-".repeat(68));
    for i in [0, n / 4, n / 2, 3 * n / 4, n - 1] {
        println!(
            "{:<18.1} | {:<15.2} | {:<15.2} | {:<15.2}",
            setpoints[i] / 1000.0,
            turbine_powers_derated[[i, 0]] / 1000.0,
            turbine_powers_derated[[i, 1]] / 1000.0,
            (turbine_powers_derated[[i, 0]] + turbine_powers_derated[[i, 1]]) / 1000.0,
        );
    }

    // Now compute available power at downstream turbine
    // Set upstream turbine setpoints to sweep, downstream set to default (unlimited)
    let mut power_setpoints_avail = ndarray::Array2::zeros((n, 2));
    for i in 0..n {
        power_setpoints_avail[[i, 0]] = setpoints[i];
        power_setpoints_avail[[i, 1]] = 1e12; // unlimited
    }

    fmodel.set_operation(
        None,
        Some(power_setpoints_avail),
        None, None, None, None,
    )?;
    fmodel.run()?;
    let turbine_powers_avail = fmodel.get_turbine_powers();

    println!("\nAvailable power at downstream turbine (unlimited setpoint):");
    println!("{:<18} | {:<15} | {:<15}", "Setpoint (kW)", "T0 Power (kW)", "T1 Available (kW)");
    println!("{}", "-".repeat(52));
    for i in [0, n / 4, n / 2, 3 * n / 4, n - 1] {
        println!(
            "{:<18.1} | {:<15.2} | {:<15.2}",
            setpoints[i] / 1000.0,
            turbine_powers_avail[[i, 0]] / 1000.0,
            turbine_powers_avail[[i, 1]] / 1000.0,
        );
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
