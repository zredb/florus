//! Example: Generate the turbulence intensity of a wind rose
//!
//! This example shows how the turbulence intensity of a wind rose or a time
//! series can be assigned: with a custom function of the wind speed and with
//! the IEC 61400-1 method.
//!
//! Scenario aligned 1:1 with
//! examples_wind_data/002_generate_ti.py.

use florus::types::{Array1, Array2};
use florus::wind_data::{TIParams, TimeSeries, WindRose};

fn main() -> anyhow::Result<()> {
    // A wind rose with three wind directions and six wind speeds
    let wind_directions = Array1::from_vec(vec![250.0, 260.0, 270.0]);
    let wind_speeds = Array1::from_vec(vec![5.0, 6.0, 7.0, 8.0, 9.0, 10.0]);

    let mut wind_rose = WindRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        Array2::from_elem((3, 6), 0.06),
        None,
        None,
        false,
        None,
        None,
    )?;

    // A custom function where TI = 1 / wind_speed
    wind_rose.assign_ti_using_wd_ws_function(|_wd: f64, ws: f64| 1.0 / ws);

    println!("Turbulence intensity from the custom function:");
    for (j, ws) in wind_speeds.iter().enumerate() {
        let ti: Vec<String> = (0..3)
            .map(|i| format!("{:.6}", wind_rose.ti_table[[i, j]]))
            .collect();
        println!("  ws = {ws:.1}: {}", ti.join(" "));
    }

    // The IEC 61400-1 method: TI is a function of the wind speed and Iref is
    // the TI value at 15 m/s
    let iref = 0.07;
    wind_rose.assign_ti_using_iec_method(Some(florus::wind_data::TIParams {
        iref,
        ..Default::default()
    }));

    println!("Turbulence intensity from the IEC method (Iref = {iref:.2}):");
    for (j, ws) in wind_speeds.iter().enumerate() {
        let ti: Vec<String> = (0..3)
            .map(|i| format!("{:.6}", wind_rose.ti_table[[i, j]]))
            .collect();
        println!("  ws = {ws:.1}: {}", ti.join(" "));
    }

    // The equivalent usage with a time series
    let n = 100;
    let mut time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_iter((0..n).map(|i| 5.0 + 10.0 * i as f64 / (n - 1) as f64)),
        Array1::from_vec(vec![0.06; n]),
    )?;
    time_series.assign_ti_using_iec_method(Some(florus::wind_data::TIParams {
        iref,
        ..Default::default()
    }));

    println!("Time series turbulence intensity (the first five time steps):");
    for i in 0..5 {
        println!(
            "  ws = {:.2}, TI = {:.6}",
            time_series.wind_speeds[i], time_series.turbulence_intensities[i]
        );
    }

    Ok(())
}
