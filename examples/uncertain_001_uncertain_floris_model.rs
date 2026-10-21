//! Example: The uncertain model parameters
//!
//! This example demonstrates the parameters of the uncertain model: the
//! resolutions of the inputs, the wind direction sample points, the standard
//! deviation of the wind direction distribution and the verbosity.
//!
//! Scenario aligned 1:1 with
//! examples_uncertain/001_uncertain_model_params.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::uncertain::UncertainFlorisModel;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // The resolution parameters define the precision of the wind direction,
    // wind speed, turbulence intensity and control parameters
    let wd_resolution = 1.0; // Degree
    let ws_resolution = 1.0; // m/s
    let ti_resolution = 0.01; // Decimal fraction
    let yaw_resolution = 1.0; // Degree
    let power_setpoint_resolution = 100.0; // kW

    // The wind direction sample points
    let wd_sample_points = vec![-6.0, -3.0, 0.0, 3.0, 6.0];

    // The standard deviation of the wind direction distribution
    let wd_std = 3.0;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut ufmodel = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml")?,
        wd_std,
        Some(wd_sample_points.clone()),
    )?;
    ufmodel.verbose = true;

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

    fmodel.run()?;
    ufmodel.run()?;

    println!("Resolutions: wd = {wd_resolution} deg, ws = {ws_resolution} m/s, ti = {ti_resolution}, yaw = {yaw_resolution} deg, power setpoint = {power_setpoint_resolution} kW");
    println!("Wind direction sample points: {wd_sample_points:?}");
    println!("Wind direction standard deviation: {wd_std} deg");

    let nominal = fmodel.get_turbine_powers()?;
    let uncertain = ufmodel.get_turbine_powers()?;
    let farm_nominal = fmodel.get_farm_power()?;
    let farm_uncertain = ufmodel.get_farm_power()?;

    println!("Wind direction [deg], nominal t0 [kW], nominal t1 [kW], uncertain t0 [kW], uncertain t1 [kW], nominal farm [kW], uncertain farm [kW]");
    for i in 0..n {
        println!(
            "{:.0}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}",
            wind_directions[i],
            nominal[[i, 0]] / 1000.0,
            nominal[[i, 1]] / 1000.0,
            uncertain[[i, 0]] / 1000.0,
            uncertain[[i, 1]] / 1000.0,
            farm_nominal[i] / 1000.0,
            farm_uncertain[i] / 1000.0
        );
    }

    Ok(())
}
