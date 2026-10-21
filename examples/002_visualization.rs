//! Example: The layout and flow visualizations
//!
//! FLORIS includes two modules for visualization:
//!   1) `layout_visualization`: the layout of the wind farm
//!   2) `flow_visualization`: the flow field
//!
//! This example produces the layout visualizations (the turbine points, the
//! labels and the waking directions) and the horizontal flow plane at hub
//! height with the turbine rotors.
//!
//! Scenario aligned 1:1 with examples/002_visualizations.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::visualization::{flow_visualization, layout_visualization};

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A farm with 8 turbines placed irregularly
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 500.0, 0.0, 128.0, 1000.0, 900.0, 1500.0, 1250.0]),
        &Array1::from_vec(vec![0.0, 300.0, 750.0, 1400.0, 0.0, 567.0, 888.0, 1450.0]),
    )?;

    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();

    // The turbine points
    layout_visualization::plot_turbine_points(
        &fmodel,
        format!("{output_dir}/002_turbine_points.png"),
        None,
        "k",
        10,
    )?;

    // The turbine points and their labels
    layout_visualization::plot_turbine_labels(
        &fmodel,
        format!("{output_dir}/002_turbine_labels.png"),
        None,
        None,
        false,
    )?;

    // The turbine points, the labels and the waking directions
    layout_visualization::plot_waking_directions(
        &fmodel,
        format!("{output_dir}/002_waking_directions.png"),
        None,
        Some(2),
    )?;

    // The flow field: the model is set to a single condition
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![290.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    fmodel.run()?;

    let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 100)?;
    flow_visualization::visualize_cut_plane(
        &horizontal_plane,
        format!("{output_dir}/002_horizontal_plane.png"),
        None,
        None,
        "viridis",
        true,
        "Horizontal flow at hub height",
    )?;

    // The turbine rotors and the labels on top of the flow field
    let turbine_names = ["T1", "T2", "T3", "T4", "T9", "T10", "T75", "T78"];
    flow_visualization::visualize_cut_plane_with_rotors_and_labels(
        &horizontal_plane,
        &fmodel,
        format!("{output_dir}/002_horizontal_plane_with_rotors.png"),
        None,
        None,
        "viridis",
        &turbine_names,
        true,
        true,
        "Horizontal flow with the turbine rotors and labels",
    )?;

    println!("The visualizations were written to {output_dir}:");
    println!("  002_turbine_points.png");
    println!("  002_turbine_labels.png");
    println!("  002_waking_directions.png");
    println!("  002_horizontal_plane.png");
    println!("  002_horizontal_plane_with_rotors.png");

    Ok(())
}
