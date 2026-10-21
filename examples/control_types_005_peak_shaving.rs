//! Example: Peak shaving
//!
//! This example compares the simple and the peak-shaving operation models:
//! with a high turbulence intensity, the peak-shaving model reduces the
//! thrust coefficient (and therefore the wake) of the upstream turbine in the
//! region around the rated wind speed.
//!
//! Scenario aligned 1:1 with
//! examples_control_types/005_peak_shaving.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;

    // A sweep of the wind speeds with a high turbulence intensity, so that
    // the peak shaving is engaged
    let wind_speeds: Vec<f64> = (0..100).map(|i| 30.0 * i as f64 / 99.0).collect();
    let n = wind_speeds.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.10; n]),
    )?;
    fmodel.set_wind_data(&time_series)?;

    // "Normal" operation with the simple turbine operation model
    fmodel.set_operation_model("simple")?;
    fmodel.run()?;
    let powers_base = fmodel.get_turbine_powers()? / 1000.0;
    let thrust_base = fmodel.get_turbine_thrust_coefficients();

    // Switch to the peak-shaving operation model
    fmodel.set_operation_model("peak-shaving")?;
    fmodel.run()?;
    let powers_peak = fmodel.get_turbine_powers()? / 1000.0;
    let thrust_peak = fmodel.get_turbine_thrust_coefficients();

    println!("Wind speed [m/s], Ct without peak shaving, Ct with peak shaving, power without [kW], power with [kW]");
    for i in 0..n {
        println!(
            "{:.2}, {:.4}, {:.4}, {:.1}, {:.1}",
            wind_speeds[i],
            thrust_base[[i, 0]],
            thrust_peak[[i, 0]],
            powers_base[[i, 0]],
            powers_peak[[i, 0]]
        );
    }

    Ok(())
}
