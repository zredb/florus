//! Example: TI Generation Methods
//!
//! Demonstrates generating turbulence intensity using different methods.
//!
//! Corresponds to: examples_wind_data/002_generate_ti.py

use florus::wind_data::{TIParams, TimeSeries, WindRose};
use florus::{FlorisModel, Result};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== TI Generation Methods ===\n");

    let wind_directions = Array1::from_vec(vec![0.0, 120.0, 240.0]);
    let wind_speeds = Array1::from_vec(vec![4.0, 6.0, 8.0, 10.0, 12.0, 14.0]);
    let ti_table = ndarray::arr2(&[
        [0.12, 0.10, 0.08, 0.07, 0.06, 0.06],
        [0.11, 0.09, 0.07, 0.07, 0.06, 0.05],
        [0.10, 0.08, 0.07, 0.06, 0.06, 0.05],
    ]);
    let frequency_table = ndarray::arr2(&[
        [0.03, 0.05, 0.10, 0.15, 0.10, 0.02],
        [0.02, 0.04, 0.08, 0.12, 0.08, 0.01],
        [0.01, 0.03, 0.06, 0.10, 0.06, 0.02],
    ]);

    let mut wind_rose = WindRose::new(
        wind_directions.clone(),
        wind_speeds.clone(),
        ti_table.clone(),
        Some(frequency_table),
        None,
        false,
        None,
        None,
    )?;

    // Method 1: Custom TI function (TI = 1/ws)
    println!("Method 1: Custom TI function (TI = 1/wind_speed)");
    let custom_ti_func = |_wd: f64, ws: f64| -> f64 { 1.0 / ws };
    wind_rose.assign_ti_using_wd_ws_function(&custom_ti_func);

    println!("  WD  |  WS  |  TI (custom)");
    println!("  ----|------|-------------");
    for (i, &wd) in wind_directions.iter().enumerate() {
        for (j, &ws) in wind_speeds.iter().enumerate() {
            println!(
                "  {:4.0} | {:4.1} | {:11.4}",
                wd,
                ws,
                wind_rose.ti_table[[i, j]]
            );
        }
    }

    // Method 2: IEC 61400-1 method
    println!("\nMethod 2: IEC 61400-1 method (Iref=0.07)");
    wind_rose.assign_ti_using_iec_method(Some(TIParams {
        iref: 0.07,
        offset: 3.8,
    }));

    println!("  WD  |  WS  |  TI (IEC)");
    println!("  ----|------|----------");
    for (i, &wd) in wind_directions.iter().enumerate() {
        for (j, &ws) in wind_speeds.iter().enumerate() {
            println!(
                "  {:4.0} | {:4.1} | {:9.4}",
                wd,
                ws,
                wind_rose.ti_table[[i, j]]
            );
        }
    }

    // Method 3: TimeSeries with IEC method
    println!("\nMethod 3: TimeSeries with IEC method");
    let mut ts = TimeSeries::new(
        Array1::from_vec(vec![270.0, 270.0, 270.0]),
        Array1::from_vec(vec![6.0, 8.0, 10.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06]),
    )?;
    ts.assign_ti_using_iec_method(Some(TIParams {
        iref: 0.07,
        offset: 3.8,
    }));

    println!("  WS  |  TI (IEC)");
    println!("  ----|----------");
    for i in 0..ts.wind_speeds.len() {
        println!(
            "  {:4.1} | {:9.4}",
            ts.wind_speeds[i], ts.turbulence_intensities[i]
        );
    }

    // Show the model with IEC TI
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(&ndarray::arr1(&[0.0, 500.0]), &ndarray::arr1(&[0.0, 0.0]))?;

    println!("\nFarm power with IEC TI values:");
    println!("{:<10} | {:<10} | {:<10}", "WS (m/s)", "Farm (kW)", "TI");
    println!("{}", "-".repeat(36));
    for j in 0..wind_speeds.len() {
        let ws_val = wind_speeds[j];
        let ti_val = wind_rose.ti_table[[2, j]];
        fmodel.set_wind_conditions(
            Array1::from_vec(vec![ws_val]),
            Array1::from_vec(vec![240.0]),
            Array1::from_vec(vec![ti_val]),
        )?;
        fmodel.run()?;
        let p = fmodel.get_turbine_powers()?;
        println!(
            "{:<10.1} | {:<10.1} | {:<10.4}",
            ws_val,
            p.row(0).sum() / 1000.0,
            ti_val
        );
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
