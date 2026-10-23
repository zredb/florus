//! Example: Visualize the flow by sweeping the turbines
//!
//! Some wake models do not have a visualization method included; for these
//! cases a slower version can be used, which scans a turbine model to produce
//! the horizontal flow field.
//!
//! Scenario aligned 1:1 with
//! examples_visualizations/005_visualize_flow_by_sweeping_turbines.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::visualization::flow_visualization;

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A two turbine layout
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0]),
        &Array1::from_vec(vec![0.0, 0.0]),
    )?;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    fmodel.run()?;

    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();

    // The horizontal plane, scanned over the turbines
    flow_visualization::visualize_horizontal_plane(
        &fmodel,
        format!("{output_dir}/005_horizontal_plane_scan.png"),
        Some(3.0),
        Some(9.0),
        Some("Horizontal (turbine scan method)"),
    )?;

    // The velocity range of the plane
    let plane = fmodel.calculate_horizontal_plane(90.0, 100, 50, None, None)?;
    let u = &plane.data.u;
    let min = u.iter().copied().fold(f64::INFINITY, f64::min);
    let max = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("Horizontal plane velocity range: {min:.3} .. {max:.3} m/s");

    Ok(())
}
