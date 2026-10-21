//! Example: The uncertain model with parallelization
//!
//! The uncertain model can be solved with the wind conditions split over
//! several workers (Python's `ParFlorisModel`). This example demonstrates the
//! usage and shows that the result obtained with and without the
//! parallelization is the same.
//!
//! Scenario aligned 1:1 with
//! examples_uncertain/003_uncertain_model_with_parallelization.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::uncertain::UncertainFlorisModel;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // The nominal model and the two uncertain models (serial and parallel)
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut ufmodel = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml")?,
        3.0,
        None,
    )?;
    let mut pufmodel = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml")?,
        3.0,
        None,
    )?;

    // The wind direction is swept while the wind speed and the turbulence
    // intensity are held constant
    let wind_directions: Vec<f64> = (240..300).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;

    // A two turbine farm
    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, d * 6.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);

    fmodel.set_layout(&layout_x, &layout_y)?;
    fmodel.set_wind_data(&time_series)?;
    ufmodel.set_layout(&layout_x, &layout_y)?;
    ufmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;
    pufmodel.set_layout(&layout_x, &layout_y)?;
    pufmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;

    // Run the models; the parallel model splits the conditions over two workers
    fmodel.run()?;
    ufmodel.run()?;
    pufmodel.run_parallel(2)?;

    let nominal = fmodel.get_turbine_powers()?;
    let serial = ufmodel.get_turbine_powers()?;
    let parallel = pufmodel.get_turbine_powers()?;

    println!("Wind direction [deg], nominal t0 [kW], nominal t1 [kW], serial t0 [kW], serial t1 [kW], parallel t0 [kW], parallel t1 [kW]");
    for i in 0..n {
        println!(
            "{:.0}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}",
            wind_directions[i],
            nominal[[i, 0]] / 1000.0,
            nominal[[i, 1]] / 1000.0,
            serial[[i, 0]] / 1000.0,
            serial[[i, 1]] / 1000.0,
            parallel[[i, 0]] / 1000.0,
            parallel[[i, 1]] / 1000.0
        );
    }

    // The parallel and the serial results are identical
    let difference = (&serial - &parallel).mapv(|v| v.abs());
    println!(
        "The maximum difference between the serial and the parallel results: {} W",
        difference.iter().copied().fold(0.0, f64::max)
    );

    Ok(())
}
