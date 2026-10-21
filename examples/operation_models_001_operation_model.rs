//! Example: Compare the yaw loss of the operation models
//!
//! This example compares the power of a yawed turbine for the different
//! operation models ("simple", "cosine-loss", "controller-dependent" and
//! "unified-momentum") at three wind speeds.
//!
//! Scenario aligned 1:1 with
//! examples_operation_models/001_compare_yaw_loss.py.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    // The yaw angle sweep
    let n = 101;
    let yaw_max = 30.0;
    let yaw_angles: Vec<f64> = (0..n)
        .map(|i| -yaw_max + 2.0 * yaw_max * i as f64 / (n - 1) as f64)
        .collect();

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![1.0; n]),
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel.set_yaw_angles(Array2::from_shape_vec(
        (n, 1),
        yaw_angles.iter().copied().collect::<Vec<f64>>(),
    )?)?;

    let operation_models = [
        "simple",
        "cosine-loss",
        "controller-dependent",
        "unified-momentum",
    ];
    let wind_speeds = [8.0, 11.5, 15.0];

    for operation_model in operation_models {
        for wind_speed in wind_speeds {
            println!("Evaluating model: {operation_model}   wind speed: {wind_speed} m/s");
            fmodel.set_wind_conditions(
                Array1::from_vec(vec![wind_speed; n]),
                Array1::from_vec(vec![270.0; n]),
                Array1::from_vec(vec![0.06; n]),
            )?;
            fmodel.set_operation_model(operation_model)?;
            fmodel.run()?;

            let powers = fmodel.get_turbine_powers()?;
            // The power at the extreme yaw angles and at zero yaw
            println!(
                "  yaw = {:+.1} deg: {:.3} kW, yaw = 0 deg: {:.3} kW, yaw = {:+.1} deg: {:.3} kW",
                yaw_angles[0],
                powers[[0, 0]] / 1000.0,
                powers[[n / 2, 0]] / 1000.0,
                yaw_angles[n - 1],
                powers[[n - 1, 0]] / 1000.0
            );
        }
    }

    Ok(())
}
