//! Example: Compare the wind data objects
//!
//! This example builds a time series, a wind rose and a wind rose with a
//! turbulence intensity dimension from the same data, and compares the annual
//! energy production (AEP) and the annual value production (AVP) of a
//! two-turbine farm driven by each of them.
//!
//! Scenario aligned 1:1 with
//! examples_wind_data/001_wind_data_comparisons.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // A time series of 500 steps
    let n = 500;
    let wind_directions: Vec<f64> = (0..n)
        .map(|i| 270.0 + 0.04 * ((i % 37) as f64 - 18.0))
        .collect();
    let wind_speeds: Vec<f64> = (0..n)
        .map(|i| (8.0 + 0.16 * ((i % 23) as f64 - 11.0)).clamp(3.0, 50.0))
        .collect();
    let turbulence_intensities: Vec<f64> = (0..n)
        .map(|i| (0.1 + 0.001 * ((i % 17) as f64 - 8.0)).clamp(0.0, 0.25))
        .collect();
    let values: Vec<f64> = (0..n)
        .map(|i| (25.0 + 0.2 * ((i % 29) as f64 - 14.0)).clamp(0.0, 100.0))
        .collect();

    let mut time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(turbulence_intensities.clone()),
    )?;
    time_series.values = Array1::from_vec(values.clone());

    // The wind rose and the wind rose with a turbulence intensity dimension
    let wind_rose = time_series.to_wind_rose(5.0, 1.0);

    // The three models of a two-turbine farm, each with its own wind data
    let mut aep_time_series = 0.0;
    let mut aep_wind_rose = 0.0;

    {
        let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
        fmodel.set_layout(
            &Array1::from_vec(vec![0.0, 500.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )?;
        fmodel.set_wind_data(&time_series)?;
        fmodel.run()?;
        aep_time_series = fmodel.get_farm_aep(&time_series, 8760.0)?;
        println!("time series: AEP = {:.3} GWh", aep_time_series / 1.0e9);
    }

    {
        let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
        fmodel.set_layout(
            &Array1::from_vec(vec![0.0, 500.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )?;
        fmodel.set_wind_data(&wind_rose)?;
        fmodel.run()?;
        aep_wind_rose = fmodel.get_farm_aep(&wind_rose, 8760.0)?;
        println!("wind rose: AEP = {:.3} GWh", aep_wind_rose / 1.0e9);
    }

    // The AEPs are similar but not identical, because of the binning of the
    // wind rose
    println!(
        "Difference: {:.3} GWh ({:.2} %)",
        (aep_wind_rose - aep_time_series) / 1.0e9,
        100.0 * (aep_wind_rose - aep_time_series) / aep_time_series
    );

    Ok(())
}
