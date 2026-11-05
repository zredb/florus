//! Example: Compare the TurbOPark implementations
//!
//! This example compares the original TurbOPark model with the TurbOPark-Gauss
//! model: the wake of a constant thrust coefficient turbine is sampled over a
//! sweep of the wind direction and along a row of ten turbines.
//!
//! Scenario aligned 1:1 with
//! examples_turbopark/001_compare_turbopark_implementations.py.

use florus::core::turbines::turbine_type::TurbineType;
use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // A constant thrust coefficient turbine (not realistic, but the power is
    // not used here)
    let const_ct_turbine = TurbineType::build(
        "ConstantCT",
        &[0.0, 30.0],
        Some(&[0.0, 1.0]),
        None,
        None,
        Some(&[0.75, 0.75]),
        None,
        100.0, // hub height
        120.0, // rotor diameter
        8.0,   // TSR
        1.225, // reference air density
        0.0,   // reference tilt
        1.88,  // cosine loss exponent yaw
        1.88,  // cosine loss exponent tilt
    )?;

    // The TurbOPark-Gauss model and the original TurbOPark model
    let mut fmodel_new = FlorisModel::from_file("examples/inputs/turboparkgauss_cubature.yaml")?;
    let mut fmodel_orig = FlorisModel::from_file("examples/inputs/turbopark_cubature.yaml")?;

    fmodel_new.set_turbine_type(&const_ct_turbine)?;
    fmodel_orig.set_turbine_type(&const_ct_turbine)?;

    // The wake of a single turbine, sampled over a sweep of the wind direction
    let wd_array: Vec<f64> = (0..900).map(|i| 225.0 + 0.1 * i as f64).collect();
    let n = wd_array.len();
    let wind_data = TimeSeries::new(
        Array1::from_vec(wd_array.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;

    fmodel_orig.set_layout(
        &Array1::from_vec(vec![0.0, 600.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    fmodel_orig.set_wind_data(&wind_data)?;
    fmodel_orig.run()?;
    let orig_vels = fmodel_orig.turbine_average_velocities();

    fmodel_new.set_layout(
        &Array1::from_vec(vec![0.0, 600.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    fmodel_new.set_wind_data(&wind_data)?;
    fmodel_new.run()?;
    let new_vels = fmodel_new.turbine_average_velocities();

    let u0 = 8.0;
    println!("Wind direction [deg], TurbOPark normalized waked speed, TurbOPark-Gauss normalized waked speed");
    for i in 0..n {
        println!(
            "{:.1}, {:.6}, {:.6}",
            wd_array[i],
            orig_vels[[i, 1]] / u0,
            new_vels[[i, 1]] / u0
        );
    }

    // The velocities along a row of ten turbines aligned with the flow
    let layout_x = Array1::from_iter((0..10).map(|i| 540.0 * i as f64));
    let layout_y = Array1::from_iter((0..10).map(|_| 0.0));
    let wind_data_row = TimeSeries::new(
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![0.06]),
    )?;

    fmodel_orig.set_layout(&layout_x, &layout_y)?;
    fmodel_orig.set_wind_data(&wind_data_row)?;
    fmodel_orig.run()?;
    let orig_row = fmodel_orig.turbine_average_velocities();

    fmodel_new.set_layout(&layout_x, &layout_y)?;
    fmodel_new.set_wind_data(&wind_data_row)?;
    fmodel_new.run()?;
    let new_row = fmodel_new.turbine_average_velocities();

    println!("Row of ten turbines, normalized velocity at every turbine:");
    println!(
        "  TurbOPark:      {:?}",
        orig_row.row(0).mapv(|v| v / u0).to_vec()
    );
    println!(
        "  TurbOPark-Gauss: {:?}",
        new_row.row(0).mapv(|v| v / u0).to_vec()
    );

    Ok(())
}
