//! Example: Setting yaw angles and disabling turbines
//!
//! This example sets the yaw angle of the most upstream turbine of a
//! two-turbine and a three-turbine farm (the middle turbine of the latter is
//! disabled) and sweeps the yaw angle from -20 to 20 degrees.
//!
//! Scenario aligned 1:1 with
//! examples_control_types/003_setting_yaw_and_disabling.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};
use florus::wind_data::TimeSeries;

fn main() -> anyhow::Result<()> {
    // Two models: a two-turbine and a three-turbine layout
    let mut fmodel_2 = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut fmodel_3 = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Change to the mixed model turbine
    fmodel_2.set_operation_model("mixed")?;
    fmodel_3.set_operation_model("mixed")?;

    fmodel_2.set_layout(
        &Array1::from_vec(vec![0.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    fmodel_3.set_layout(
        &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    // Constant wind conditions
    let n = 50;
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel_2.set_wind_data(&time_series)?;
    fmodel_3.set_wind_data(&time_series)?;

    // The yaw angle of the most upstream turbine is swept from -20 to 20
    // degrees, the other turbines are set to 0
    let upstream_yaw_angles: Vec<f64> = (0..n)
        .map(|i| -20.0 + 40.0 * i as f64 / (n - 1) as f64)
        .collect();
    let mut yaw_angles_2 = Array2::zeros((n, 2));
    let mut yaw_angles_3 = Array2::zeros((n, 3));
    for (i, yaw) in upstream_yaw_angles.iter().enumerate() {
        yaw_angles_2[[i, 0]] = *yaw;
        yaw_angles_3[[i, 0]] = *yaw;
    }

    // The middle turbine of the three-turbine case is disabled
    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((n, 3), false);
    for i in 0..n {
        disable_turbines[[i, 1]] = true;
    }

    fmodel_2.set_operation(Some(yaw_angles_2), None, None, None, None, None)?;
    fmodel_3.set_operation(
        Some(yaw_angles_3),
        None,
        None,
        None,
        None,
        Some(disable_turbines),
    )?;

    fmodel_2.run()?;
    fmodel_3.run()?;

    let powers_2 = fmodel_2.get_turbine_powers()?;
    let powers_3 = fmodel_3.get_turbine_powers()?;

    println!("Yaw [deg], 2-turbine upstream [kW], 2-turbine downstream [kW], 3-turbine upstream [kW], 3-turbine downstream [kW]");
    for i in 0..n {
        println!(
            "{:.2}, {:.1}, {:.1}, {:.1}, {:.1}",
            upstream_yaw_angles[i],
            powers_2[[i, 0]] / 1000.0,
            powers_2[[i, 1]] / 1000.0,
            powers_3[[i, 0]] / 1000.0,
            powers_3[[i, 2]] / 1000.0
        );
    }

    Ok(())
}
