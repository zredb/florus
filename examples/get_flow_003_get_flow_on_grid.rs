//! Example: Extract the turbulence intensity at points
//!
//! This example demonstrates the use of the `sample_ti_at_points` method: the
//! turbulence intensity is extracted at user-specified locations in the flow.
//! The TI is sampled at a single x, y location and four different heights
//! over a sweep of wind directions, which mimics the TI measurements of a met
//! mast.
//!
//! Scenario aligned 1:1 with
//! examples_get_flow/003_extract_turbulence_intensity_at_points.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    // The met mast location (0 to 3, as in the Python example)
    let met_mast_option = 0;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A two-turbine farm
    let d = 126.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 3.0 * d]),
        &Array1::from_vec(vec![0.0, 3.0 * d]),
    )?;

    // The wind direction is swept over 360 degrees
    let wd_array: Vec<f64> = (0..360).map(|i| i as f64).collect();
    let n = wd_array.len();
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(wd_array.clone()),
        Array1::from_vec(vec![0.06; n]),
    )?;

    // The met mast in between the turbines
    let (points_x, points_y) = match met_mast_option {
        0 => (vec![3.0 * d; 4], vec![0.0; 4]),
        1 => (vec![200.0; 4], vec![200.0; 4]),
        2 => (vec![20.0; 4], vec![20.0; 4]),
        _ => (vec![305.0; 4], vec![158.0; 4]),
    };
    let points_z = [30.0, 90.0, 150.0, 250.0];

    let ti_at_points = fmodel.sample_ti_at_points(&points_x, &points_y, &points_z)?;

    println!("Wind direction [deg], TI at z=30 m, z=90 m, z=150 m, z=250 m");
    for i in 0..n {
        println!(
            "{:.0}, {:.6}, {:.6}, {:.6}, {:.6}",
            wd_array[i],
            ti_at_points[[i, 0]],
            ti_at_points[[i, 1]],
            ti_at_points[[i, 2]],
            ti_at_points[[i, 3]]
        );
    }

    Ok(())
}
