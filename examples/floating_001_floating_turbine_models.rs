//! Example: The floating turbine models
//!
//! This example compares three single-turbine simulations:
//!   * `fixed`: a fixed bottom turbine (no tilt variation with the wind speed)
//!   * `floating`: a floating turbine whose tilt varies with the wind speed
//!   * `defined_floating`: a floating turbine whose tilt varies with the wind
//!     speed, but whose power/thrust tables already account for the tilt
//!
//! Scenario aligned 1:1 with examples_floating/001_floating_turbine_models.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // The wind speeds, with a constant wind direction
    let wind_speeds: Vec<f64> = (3..25).map(|i| i as f64).collect();
    let n = wind_speeds.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;

    let mut models: Vec<(&str, FlorisModel)> = vec![
        (
            "fixed",
            FlorisModel::from_file("examples/inputs/gch_fixed.yaml")?,
        ),
        (
            "floating",
            FlorisModel::from_file("examples/inputs/gch_floating.yaml")?,
        ),
        (
            "defined_floating",
            FlorisModel::from_file("examples/inputs/gch_floating_defined_floating.yaml")?,
        ),
    ];

    let mut powers: Vec<Vec<f64>> = Vec::new();
    let mut cts: Vec<Vec<f64>> = Vec::new();
    for (_, fmodel) in models.iter_mut() {
        fmodel.set_wind_data(&time_series)?;
        fmodel.run()?;
        let turbine_powers = fmodel.get_turbine_powers()?;
        powers.push(
            (0..n)
                .map(|i| turbine_powers[[i, 0]] / 1000.0)
                .collect::<Vec<f64>>(),
        );
        let turbine_cts = fmodel.get_turbine_thrust_coefficients();
        cts.push((0..n).map(|i| turbine_cts[[i, 0]]).collect::<Vec<f64>>());
    }

    println!("Wind speed [m/s], fixed [kW], floating [kW], defined floating [kW], fixed Ct, floating Ct, defined floating Ct");
    for i in 0..n {
        println!(
            "{:.1}, {:.3}, {:.3}, {:.3}, {:.6}, {:.6}, {:.6}",
            wind_speeds[i],
            powers[0][i],
            powers[1][i],
            powers[2][i],
            cts[0][i],
            cts[1][i],
            cts[2][i]
        );
    }

    Ok(())
}
