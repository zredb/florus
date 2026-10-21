//! Example: The uncertain model
//!
//! This example compares the nominal farm power with the power of a farm
//! whose wind direction has an uncertainty: the wind direction of every
//! condition is perturbed with the sample points of a Gaussian distribution
//! and the results are averaged with the Gaussian weights.
//!
//! Scenario aligned 1:1 with examples/008_uncertain_models.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::uncertain::UncertainFlorisModel;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // The nominal model and the two uncertain models (a standard deviation of
    // 3 and 5 degrees)
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut ufmodel_3 = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml")?,
        3.0,
        None,
    )?;
    let mut ufmodel_5 = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml")?,
        5.0,
        None,
    )?;

    // An inflow where the wind direction is swept while the wind speed and the
    // turbulence intensity are held constant
    let wind_directions: Vec<f64> = (240..300).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;

    // A two turbine farm with the inflow
    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, d * 6.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);

    fmodel.set_layout(&layout_x, &layout_y)?;
    fmodel.set_wind_data(&time_series)?;
    ufmodel_3.set_layout(&layout_x, &layout_y)?;
    ufmodel_5.set_layout(&layout_x, &layout_y)?;
    ufmodel_3.set_wind_conditions(
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;
    ufmodel_5.set_wind_conditions(
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;

    // Run the models
    fmodel.run()?;
    ufmodel_3.run()?;
    ufmodel_5.run()?;

    let powers_nom = fmodel.get_turbine_powers()?;
    let powers_3 = ufmodel_3.get_turbine_powers()?;
    let powers_5 = ufmodel_5.get_turbine_powers()?;

    println!("Wind direction [deg], nominal t0 [kW], nominal t1 [kW], uncertain(3 deg) t0 [kW], uncertain(3 deg) t1 [kW], uncertain(5 deg) t0 [kW], uncertain(5 deg) t1 [kW]");
    for i in 0..n {
        println!(
            "{:.0}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}",
            wind_directions[i],
            powers_nom[[i, 0]] / 1000.0,
            powers_nom[[i, 1]] / 1000.0,
            powers_3[[i, 0]] / 1000.0,
            powers_3[[i, 1]] / 1000.0,
            powers_5[[i, 0]] / 1000.0,
            powers_5[[i, 1]] / 1000.0
        );
    }

    Ok(())
}
