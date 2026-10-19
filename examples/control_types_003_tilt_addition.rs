//! Example: Setting Yaw Angles and Disabling Turbine
//!
//! Demonstrates combining yaw control with turbine disabling.
//! A 3-turbine case with the middle turbine disabled produces the same
//! results as a 2-turbine case.
//!
//! Corresponds to: examples_control_types/003_setting_yaw_and_disabling.py

use florus::wind_data::TimeSeries;
use florus::{FlorisModel, Result};
use ndarray::Array1;

fn main() -> Result<()> {
    println!("=== Yaw + Disable Example ===\n");

    let n = 50;
    let upstream_yaw: Vec<f64> = (0..n)
        .map(|i| -20.0 + 40.0 * i as f64 / (n - 1) as f64)
        .collect();

    // 2-turbine model
    let mut fmodel_2 = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_2.set_layout(&ndarray::arr1(&[0.0, 1000.0]), &ndarray::arr1(&[0.0, 0.0]))?;
    fmodel_2.set_operation_model("mixed")?;
    let ts2 = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel_2.set_wind_data(&ts2)?;

    // 3-turbine model (middle turbine will be disabled)
    let mut fmodel_3 = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_3.set_layout(
        &ndarray::arr1(&[0.0, 500.0, 1000.0]),
        &ndarray::arr1(&[0.0, 0.0, 0.0]),
    )?;
    fmodel_3.set_operation_model("mixed")?;
    let ts3 = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )?;
    fmodel_3.set_wind_data(&ts3)?;

    // Set yaw angles
    let yaw_2: Vec<f64> = upstream_yaw.iter().flat_map(|&yaw| [yaw, 0.0]).collect();
    let yaw_3: Vec<f64> = upstream_yaw
        .iter()
        .flat_map(|&yaw| [yaw, 0.0, 0.0])
        .collect();
    fmodel_2.set_yaw_angles(ndarray::Array2::from_shape_vec((n, 2), yaw_2)?)?;
    fmodel_3.set_yaw_angles(ndarray::Array2::from_shape_vec((n, 3), yaw_3)?)?;

    // Disable middle turbine in 3-turbine case
    let mut disable = ndarray::Array2::from_elem((n, 3), false);
    for i in 0..n {
        disable[[i, 1]] = true;
    }
    fmodel_3.set_operation(None, None, None, None, None, Some(disable))?;

    fmodel_2.run()?;
    fmodel_3.run()?;

    let p2 = fmodel_2.get_turbine_powers()? / 1000.0;
    let p3 = fmodel_3.get_turbine_powers()? / 1000.0;

    println!("Upstream turbine power comparison:");
    println!(
        "{:<12} | {:<15} | {:<15} | {:<12}",
        "Yaw (deg)", "2T Power (kW)", "3T Power (kW)", "Match?"
    );
    println!("{}", "-".repeat(58));
    for i in [0, n / 4, n / 2, 3 * n / 4, n - 1] {
        let diff = (p2[[i, 0]] - p3[[i, 0]]).abs();
        println!(
            "{:<12.1} | {:<15.2} | {:<15.2} | {:<12}",
            upstream_yaw[i],
            p2[[i, 0]],
            p3[[i, 0]],
            if diff < 0.01 { "Yes" } else { "No" }
        );
    }

    println!("\nDownstream turbine power comparison:");
    println!(
        "{:<12} | {:<15} | {:<15} | {:<12}",
        "Yaw (deg)", "2T Power (kW)", "3T Power (kW)", "Match?"
    );
    println!("{}", "-".repeat(58));
    for i in [0, n / 4, n / 2, 3 * n / 4, n - 1] {
        let diff = (p2[[i, 1]] - p3[[i, 2]]).abs();
        println!(
            "{:<12.1} | {:<15.2} | {:<15.2} | {:<12}",
            upstream_yaw[i],
            p2[[i, 1]],
            p3[[i, 2]],
            if diff < 0.01 { "Yes" } else { "No" }
        );
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
