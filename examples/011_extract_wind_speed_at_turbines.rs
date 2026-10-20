//! Example: Extract Wind Speed at Turbines
//!
//! This example demonstrates how to extract the wind speed at the turbine
//! points from the FLORUS model. Both the u velocities and the turbine
//! average velocities are grabbed from the model, then the turbine average
//! is recalculated from the u velocities to show that they are equivalent.
//!
//! Scenario aligned 1:1 with examples_get_flow/001_extract_wind_speed_at_turbines.py
//! (4-turbine 2 x 2 layout, 1 findex, turbine_grid_points = 3).

use florus::{FlorisModel, Result};
use ndarray::Array;

fn main() -> Result<()> {
    // Initialize the FLORUS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Create a 4-turbine layout (Python: layout_x=[0, 0, 500, 500],
    // layout_y=[0, 300, 0, 300])
    fmodel.set_layout(
        &Array::from_vec(vec![0.0, 0.0, 500.0, 500.0]),
        &Array::from_vec(vec![0.0, 300.0, 0.0, 300.0]),
    )?;

    // Calculate wake
    fmodel.run()?;

    // Collect the wind speed at all the turbine points
    let u_points = fmodel.core().flow_field.u.clone();

    println!("U points is 1 findex x 4 turbines x 3 x 3 points (turbine_grid_points=3)");
    let s = u_points.shape();
    println!("({}, {}, {}, {})", s[0], s[1], s[2], s[3]);

    println!("turbine_average_velocities is 1 findex x 4 turbines");
    let avg = fmodel.turbine_average_velocities();
    println!("[[{}]]", format_row(&avg));

    // Show that one is equivalent to the other following averaging
    println!(
        "turbine_average_velocities is determined by taking the cube root of mean \
         of the cubed value across the points "
    );
    println!("turbine_average_velocities: [[{}]]", format_row(&avg));
    let recomputed = recompute_cubic_mean(&u_points);
    println!("Recomputed:       [[{}]]", format_row(&recomputed));

    Ok(())
}

fn format_row(avg: &ndarray::Array2<f64>) -> String {
    // numpy-style float printing: 8 decimals, trailing zeros trimmed
    // ("4.78455500" -> "4.784555", "7.97368581" kept).
    avg.iter()
        .map(|v| {
            let mut s = format!("{v:.8}");
            if s.contains('.') {
                while s.ends_with('0') {
                    s.pop();
                }
                if s.ends_with('.') {
                    s.pop();
                }
            }
            s
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn recompute_cubic_mean(u: &ndarray::Array4<f64>) -> ndarray::Array2<f64> {
    let u_cubed = u.mapv(|x: f64| x.powi(3));
    let mean_axis2 = u_cubed.mean_axis(ndarray::Axis(2)).unwrap();
    let mean = mean_axis2.mean_axis(ndarray::Axis(2)).unwrap();
    mean.mapv(|x: f64| x.cbrt())
}
