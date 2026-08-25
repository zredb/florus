//! Example: Compare Yaw Loss Under Different Operation Models
//!
//! Demonstrates how different operation models alter power loss to yaw.
//!
//! Corresponds to: examples_operation_models/001_compare_yaw_loss.py

use florus::{FlorisModel, Result};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Compare Yaw Loss Under Different Operation Models ===\n");

    let n_yaw = 51;
    let yaw_max = 30.0;
    let wind_speeds = [8.0, 11.5, 15.0];
    let op_models = ["simple", "cosine-loss", "controller-dependent", "unified-momentum"];

    let yaw_angles: Vec<f64> = (0..n_yaw)
        .map(|i| -yaw_max + 2.0 * yaw_max * i as f64 / (n_yaw - 1) as f64)
        .collect();

    for &ws in &wind_speeds {
        println!("Wind speed: {} m/s", ws);
        println!("{:<20} | {:<15} | {:<15}", "Operation Model", "Max Power (kW)", "Power at 30deg");
        println!("{}", "-".repeat(55));

        for op_model in &op_models {
            let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
            fmodel.set_layout(
                &ndarray::arr1(&[0.0]),
                &ndarray::arr1(&[0.0]),
            )?;
            fmodel.set_operation_model(*op_model)?;

            let wind_directions = Array1::from_vec(vec![270.0; n_yaw]);
            let wind_speeds_arr = Array1::from_vec(vec![ws; n_yaw]);
            let turbulence_intensities = Array1::from_vec(vec![0.06; n_yaw]);

            fmodel.set_wind_conditions(wind_speeds_arr, wind_directions, turbulence_intensities)?;

            fmodel.set_yaw_angles(ndarray::Array2::from_shape_vec(
                (n_yaw, 1),
                yaw_angles.clone(),
            )?)?;

            fmodel.run()?;

            let powers = fmodel.get_turbine_powers() / 1000.0;

            let max_power = powers.column(0).fold(0.0_f64, |a, &b| a.max(b));
            let power_at_30 = powers[[n_yaw - 1, 0]];

            println!("{:<20} | {:<15.2} | {:<15.2}", op_model, max_power, power_at_30);
        }
        println!();
    }

    println!("=== Example Complete ===");
    Ok(())
}
