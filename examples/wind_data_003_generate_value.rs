//! Example: Value Generation Methods
//!
//! Demonstrates generating electricity value tables using custom function
//! and piecewise linear methods.
//!
//! Corresponds to: examples_wind_data/003_generate_value.py

use florus::wind_data::WindRose;
use ndarray::Array1;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Value Generation Methods ===\n");

    let wind_directions = Array1::from_vec(vec![0.0, 120.0, 240.0]);
    let wind_speeds = Array1::from_vec(vec![4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]);
    let ti_table = ndarray::arr2(&[
        [0.12, 0.10, 0.08, 0.07, 0.06, 0.06, 0.05, 0.05],
        [0.11, 0.09, 0.07, 0.07, 0.06, 0.05, 0.05, 0.04],
        [0.10, 0.08, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04],
    ]);
    let frequency_table = ndarray::arr2(&[
        [0.02, 0.04, 0.08, 0.12, 0.08, 0.04, 0.02, 0.01],
        [0.01, 0.03, 0.06, 0.10, 0.06, 0.03, 0.01, 0.00],
        [0.01, 0.02, 0.05, 0.08, 0.05, 0.02, 0.01, 0.01],
    ]);

    let mut wind_rose = WindRose::new(
        wind_directions,
        wind_speeds.clone(),
        ti_table,
        Some(frequency_table),
        None,
        false,
        None,
        None,
    )?;

    // Method 1: Custom value function (value = 100/ws)
    println!("Method 1: Custom value function (value = 100/wind_speed)");
    let custom_value_func = |_wd: f64, ws: f64| -> f64 { 100.0 / ws };
    wind_rose.assign_value_using_wd_ws_function(&custom_value_func, false);

    if let Some(ref vt) = wind_rose.value_table {
        println!("  WS (m/s) | Value ($/MWh)");
        println!("  ---------|-------------");
        for (j, &ws) in wind_speeds.iter().enumerate() {
            println!("  {:8.1} | {:11.2}", ws, vt[[0, j]]);
        }
    }

    // Method 2: Piecewise linear
    println!("\nMethod 2: Piecewise linear (SPP market approximation)");
    wind_rose.assign_value_piecewise_linear(1.425, 4.5, 0.0, -0.135, false, false);

    if let Some(ref vt) = wind_rose.value_table {
        println!("  WS (m/s) | Value ($/MWh)");
        println!("  ---------|-------------");
        for (j, &ws) in wind_speeds.iter().enumerate() {
            println!("  {:8.1} | {:11.2}", ws, vt[[0, j]]);
        }
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
