//! Example: Generate the value of a wind rose
//!
//! This example shows how the value (e.g. the electricity price) of a wind
//! rose or a time series can be assigned: with a custom function of the wind
//! speed and with the piecewise linear function of the wind speed.
//!
//! Scenario aligned 1:1 with
//! examples_wind_data/003_generate_value.py.

use florus::types::{Array1, Array2};
use florus::wind_data::{TimeSeries, WindRose};

fn main() -> anyhow::Result<()> {
    // A wind rose with three wind directions and eight wind speeds
    let wind_directions = Array1::from_vec(vec![250.0, 260.0, 270.0]);
    let wind_speeds = Array1::from_iter((3..11).map(|i| i as f64));

    let mut wind_rose = WindRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        Array2::from_elem((3, 8), 0.06),
        None,
        None,
        false,
        None,
        None,
    )?;

    // A custom function where value = 100 / wind_speed
    wind_rose.assign_value_using_wd_ws_function(|_wd: f64, ws: f64| 100.0 / ws, false);

    println!("Value from the custom function:");
    for (j, ws) in wind_speeds.iter().enumerate() {
        let values: Vec<String> = (0..3)
            .map(|i| format!("{:.6}", wind_rose.value_table.as_ref().unwrap()[[i, j]]))
            .collect();
        println!("  ws = {ws:.1}: {}", values.join(" "));
    }

    // The piecewise linear function of the wind speed (the default settings
    // approximate the normalized mean electricity price of the SPP market)
    wind_rose.assign_value_piecewise_linear(1.425, 4.5, 0.0, -0.135, false, false);

    println!("Value from the piecewise linear function:");
    for (j, ws) in wind_speeds.iter().enumerate() {
        let values: Vec<String> = (0..3)
            .map(|i| format!("{:.6}", wind_rose.value_table.as_ref().unwrap()[[i, j]]))
            .collect();
        println!("  ws = {ws:.1}: {}", values.join(" "));
    }

    // The equivalent usage with a time series
    let n = 100;
    let mut time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_iter((0..n).map(|i| 3.0 + 10.0 * i as f64 / (n - 1) as f64)),
        Array1::from_vec(vec![0.06; n]),
    )?;
    time_series.assign_value_using_wd_ws_function(|_wd: f64, ws: f64| 100.0 / ws, false);

    println!("Time series values (the first five time steps):");
    for i in 0..5 {
        println!(
            "  ws = {:.2}, value = {:.6}",
            time_series.wind_speeds[i], time_series.values[i]
        );
    }

    Ok(())
}
