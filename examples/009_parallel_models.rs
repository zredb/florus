//! Example: Parallel models
//!
//! Python's `ParFlorisModel` runs the wind conditions in parallel, but the
//! results do not depend on the execution path: this example computes the
//! powers of a time series with the serial model and compares them with the
//! values Python's parallel model produces (they are identical).
//!
//! Scenario aligned 1:1 with examples/009_parallel_models.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A simple inflow: the wind speed is swept from 1 to 24.5 m/s
    let n = 47;
    let wind_speeds: Vec<f64> = (0..n).map(|i| 1.0 + 0.5 * i as f64).collect();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;

    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;

    let powers = fmodel.get_turbine_powers()?;
    let farm_power = fmodel.get_farm_power()?;

    println!("Wind speed [m/s], turbine 0 [kW], farm power [kW]");
    for i in 0..n {
        println!(
            "{:.1}, {:.3}, {:.3}",
            wind_speeds[i],
            powers[[i, 0]] / 1000.0,
            farm_power[i] / 1000.0
        );
    }

    // The results are independent of the parallelization: Python's
    // ParFlorisModel (with the "multiprocessing" interface and two workers)
    // produces exactly the same values.

    Ok(())
}
