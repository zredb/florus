//! Example: Getting the power of the farm
//!
//! This example demonstrates the different ways of getting the power out of
//! a FLORIS model: the turbine powers, the farm power, the expected farm
//! power and the annual energy production (AEP), with and without wakes.
//!
//! Scenario aligned 1:1 with examples/005_getting_power.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A three-turbine layout
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 126.0 * 5.0, 126.0 * 10.0]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // A time series in which the wind speed and the turbulence intensity are
    // constant but the wind direction sweeps from 250 to 290 degrees
    let wind_directions: Vec<f64> = (250..290).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![9.9; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.set_wind_data(&time_series)?;

    fmodel.run()?;

    let turbine_powers = fmodel.get_turbine_powers()?;
    let farm_power = fmodel.get_farm_power()?;

    // The same without wakes, which is useful to compute the wake losses
    fmodel.run_no_wake()?;
    let farm_power_no_wake = fmodel.get_farm_power()?;

    println!(
        "Turbine power has shape ({}, {})",
        turbine_powers.shape()[0],
        turbine_powers.shape()[1]
    );
    println!("Farm power has shape ({},)", farm_power.shape()[0]);

    println!("Wind direction [deg], turbine 0 [kW], turbine 1 [kW], turbine 2 [kW], farm power [kW], farm power without wakes [kW], wake losses [%]");
    for i in 0..n {
        let wake_losses = 100.0 * (farm_power_no_wake[i] - farm_power[i]) / farm_power_no_wake[i];
        println!(
            "{:.0}, {:.1}, {:.1}, {:.1}, {:.1}, {:.1}, {:.2}",
            wind_directions[i],
            turbine_powers[[i, 0]] / 1000.0,
            turbine_powers[[i, 1]] / 1000.0,
            turbine_powers[[i, 2]] / 1000.0,
            farm_power[i] / 1000.0,
            farm_power_no_wake[i] / 1000.0,
            wake_losses
        );
    }

    Ok(())
}
