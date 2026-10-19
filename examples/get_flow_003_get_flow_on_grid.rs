//! Example: Extract Turbulence Intensity at Points
//!
//! Demonstrates sample_ti_at_points to extract TI at user-specified locations.
//!
//! Corresponds to: examples_get_flow/003_extract_turbulence_intensity_at_points.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Extract Turbulence Intensity at Points ===\n");

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let d = 126.0;
    fmodel.set_layout(&ndarray::arr1(&[0.0, 3.0 * d]), &ndarray::arr1(&[0.0, 0.0]))?;

    let points_x = vec![3.0 * d; 4];
    let points_y = vec![0.0; 4];
    let points_z = vec![30.0, 90.0, 150.0, 250.0];

    println!(
        "Met Mast at ({:.0}, 0) with heights: {:?} m\n",
        3.0 * d,
        points_z
    );

    println!("TI at Met Mast:");
    println!(
        "{:<10} | {:<10} | {:<10} | {:<10} | {:<10}",
        "WD (deg)", "z=30m", "z=90m", "z=150m", "z=250m"
    );
    println!("{}", "-".repeat(58));

    for wd in (0..360).step_by(30) {
        fmodel.set_wind_conditions(
            ndarray::arr1(&[8.0]),
            ndarray::arr1(&[wd as f64]),
            ndarray::arr1(&[0.06]),
        )?;
        fmodel.run()?;

        let ti_at_points = fmodel.sample_ti_at_points(&points_x, &points_y, &points_z)?;

        // Returns (n_findex, n_points); the model has a single condition.
        println!(
            "{:<10} | {:<10.4} | {:<10.4} | {:<10.4} | {:<10.4}",
            wd,
            ti_at_points[[0, 0]],
            ti_at_points[[0, 1]],
            ti_at_points[[0, 2]],
            ti_at_points[[0, 3]]
        );
    }

    println!("\n=== Example Complete ===");
    Ok(())
}
