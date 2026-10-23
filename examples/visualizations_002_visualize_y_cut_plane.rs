//! Example: Visualize the y cut plane
//!
//! This example visualizes the y cut plane of a three-turbine farm with the
//! wind direction along the row.
//!
//! Scenario aligned 1:1 with
//! examples_visualizations/002_visualize_y_cut_plane.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::visualization::flow_visualization;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A three turbine layout with the wind direction along the row
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    fmodel.run()?;

    // The y plane
    let y_plane = fmodel.calculate_y_plane(200, 100, 0.0, None, None)?;

    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();
    flow_visualization::visualize_cut_plane(
        &y_plane,
        format!("{output_dir}/002_y_cut_plane.png"),
        Some(3.0),
        Some(9.0),
        "coolwarm",
        false,
        "Y Cut Plane",
    )?;

    // The velocity range of the plane
    let u = &y_plane.data.u;
    let min = u.iter().copied().fold(f64::INFINITY, f64::min);
    let max = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("Y cut plane velocity range: {min:.3} .. {max:.3} m/s");

    Ok(())
}
