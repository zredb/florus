//! Example: The wind data objects
//!
//! FLORIS provides a set of wind data objects to hold the ambient wind
//! conditions: `TimeSeries`, `WindRose` and `WindTIRose`. This example
//! shows how they are initialized, how values are broadcast, how they are
//! converted into each other and how they are applied to a model.
//!
//! Scenario aligned 1:1 with examples/003_wind_data_objects.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array3};
use florus::wind_data::{TimeSeries, WindRose, WindTIRose};

fn main() -> anyhow::Result<()> {
    // ---------------------------------------------------------------
    // Initializing
    // ---------------------------------------------------------------

    // A time series of 100 steps
    let n = 100;
    let wind_speeds: Vec<f64> = (0..n).map(|i| 8.0 + 0.1 * i as f64).collect();
    let wind_directions: Vec<f64> = (0..n).map(|i| 270.0 + 0.2 * i as f64).collect();
    let turbulence_intensities: Vec<f64> = (0..n).map(|i| 0.06 + 0.0001 * i as f64).collect();
    let values: Vec<f64> = (0..n).map(|i| 25.0 + 0.05 * i as f64).collect();

    let mut time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(turbulence_intensities.clone()),
    )?;
    time_series.values = Array1::from_vec(values.clone());

    println!("TimeSeries with {n} time steps");
    println!(
        "  wind speeds: {:.2} .. {:.2} m/s",
        wind_speeds[0],
        wind_speeds[n - 1]
    );
    println!(
        "  wind directions: {:.2} .. {:.2} deg",
        wind_directions[0],
        wind_directions[n - 1]
    );
    println!(
        "  turbulence intensities: {:.4} .. {:.4}",
        turbulence_intensities[0],
        turbulence_intensities[n - 1]
    );
    println!("  values: {:.2} .. {:.2}", values[0], values[n - 1]);

    // ---------------------------------------------------------------
    // Converting between the objects
    // ---------------------------------------------------------------

    // A time series can be converted into a wind rose with a 5 degree and a
    // 1 m/s resolution
    let wind_rose = time_series.to_wind_rose(5.0, 1.0);
    println!("WindRose converted from the time series:");
    println!(
        "  {} wind directions, {} wind speeds",
        wind_rose.wind_directions.len(),
        wind_rose.wind_speeds.len()
    );

    // The wind rose can be converted back into a time series
    let time_series_2 = wind_rose.to_time_series();
    println!(
        "TimeSeries converted from the wind rose: {} steps",
        time_series_2.wind_speeds.len()
    );

    // A wind rose with a turbulence intensity dimension
    let wind_ti_rose = WindTIRose::new(
        Array1::from_vec(vec![270.0, 280.0]),
        Array1::from_vec(vec![8.0, 9.0]),
        Array1::from_vec(vec![0.06, 0.08]),
        Array3::zeros((2, 2, 2)),
        None,
        None,
    )?;
    println!(
        "WindTIRose: {} directions, {} speeds, {} turbulence intensities",
        wind_ti_rose.wind_directions.len(),
        wind_ti_rose.wind_speeds.len(),
        wind_ti_rose.turbulence_intensities.len()
    );

    // ---------------------------------------------------------------
    // Applying the wind data to a model
    // ---------------------------------------------------------------

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 630.0, 1260.0]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // With the time series
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Turbine powers with the time series [kW] (the first three steps):");
    for i in 0..3 {
        println!(
            "  {:.1}, {:.1}, {:.1}",
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0,
            powers[[i, 2]] / 1000.0
        );
    }

    // With the wind rose
    fmodel.set_wind_data(&wind_rose)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    println!("Turbine powers with the wind rose [kW] (the first three conditions):");
    for i in 0..3 {
        println!(
            "  {:.1}, {:.1}, {:.1}",
            powers[[i, 0]] / 1000.0,
            powers[[i, 1]] / 1000.0,
            powers[[i, 2]] / 1000.0
        );
    }

    Ok(())
}
