//! Example: The reference turbines of the turbine library
//!
//! This example compares the power and the thrust coefficient of every
//! turbine definition of the FLORIS turbine library over a sweep of the wind
//! speed, and the effect of the yaw angle on the power of one turbine.
//!
//! Scenario aligned 1:1 with
//! examples_turbine/001_reference_turbines.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    // The wind speeds and the yaw angles of the comparison
    let ws_array: Vec<f64> = (1..150).map(|i| 0.1 + 0.2 * i as f64).collect();
    let n = ws_array.len();
    let wd_array = vec![270.0; n];
    let turbulence_intensities = vec![0.06; n];

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A single turbine
    fmodel.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;

    // The available turbine definitions of the library
    let turbines = fmodel.available_turbine_types(true);

    println!("Turbine, wind speed [m/s], power [kW], thrust coefficient");
    for turbine in &turbines {
        // Set the turbine type and reset the operation settings
        fmodel.set_turbine_types(&[turbine.clone()])?;
        fmodel.reset_operation();

        // The reference wind height follows the hub height of the turbine
        fmodel.assign_hub_height_to_ref_height()?;

        fmodel.set_wind_conditions(
            Array1::from_vec(ws_array.clone()),
            Array1::from_vec(wd_array.clone()),
            Array1::from_vec(turbulence_intensities.clone()),
        )?;
        fmodel.run()?;

        let powers = fmodel.get_turbine_powers()?;
        let cts = fmodel.get_turbine_thrust_coefficients();
        for i in 0..n {
            println!(
                "{}, {:.1}, {:.3}, {:.6}",
                turbine,
                ws_array[i],
                powers[[i, 0]] / 1000.0,
                cts[[i, 0]]
            );
        }
    }

    // The effect of the yaw angle on the power of one turbine
    let n_yaw = 60;
    let yaw_angles: Vec<f64> = (0..n_yaw)
        .map(|i| -30.0 + 60.0 * i as f64 / (n_yaw - 1) as f64)
        .collect();
    fmodel.set_turbine_types(&["nrel_5MW".to_string()])?;
    fmodel.reset_operation();
    fmodel.assign_hub_height_to_ref_height()?;

    // The wind speed at which the yaw effect is tested
    let wind_speed_to_test_yaw = 11.0;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![wind_speed_to_test_yaw; n_yaw]),
        Array1::from_vec(vec![270.0; n_yaw]),
        Array1::from_vec(vec![0.06; n_yaw]),
    )?;
    fmodel.set_yaw_angles(Array2::from_shape_vec(
        (n_yaw, 1),
        yaw_angles.iter().copied().collect::<Vec<f64>>(),
    )?)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;

    println!("Yaw angle [deg], power [kW]");
    for i in 0..n_yaw {
        println!("{:.3}, {:.3}", yaw_angles[i], powers[[i, 0]] / 1000.0);
    }

    Ok(())
}
