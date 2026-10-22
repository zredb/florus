//! Example: Visualize the cross plane
//!
//! This example visualizes the cross plane downstream of a single turbine.
//!
//! Scenario aligned 1:1 with
//! examples_visualizations/003_visualize_cross_plane.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::visualization::flow_visualization;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A single turbine
    fmodel.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    fmodel.run()?;

    // The cross plane 500 m downstream of the turbine
    let cross_plane = fmodel.calculate_cross_plane(100, 100, 500.0)?;

    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();
    flow_visualization::visualize_cut_plane(
        &cross_plane,
        format!("{output_dir}/003_cross_plane.png"),
        Some(3.0),
        Some(9.0),
        "coolwarm",
        false,
        "Cross Plane",
    )?;

    // The velocity range of the plane
    let u = &cross_plane.data.u;
    let min = u.iter().copied().fold(f64::INFINITY, f64::min);
    let max = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("Cross plane velocity range: {min:.3} .. {max:.3} m/s");

    Ok(())
}
