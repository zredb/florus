//! Example: The layout visualizations
//!
//! This example produces the layout visualizations of a five-turbine farm: the
//! turbine points, the turbine labels, the waking directions and the flow
//! field at hub height.
//!
//! Scenario aligned 1:1 with
//! examples_visualizations/001_layout_visualizations.py.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::visualization::{flow_visualization, layout_visualization};

fn main() -> anyhow::Result<()> {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // A five-turbine layout with the wind from the north west
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 0.0, 1000.0, 1000.0, 1000.0]),
        &Array1::from_vec(vec![0.0, 500.0, 0.0, 500.0, 1000.0]),
    )?;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![300.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    fmodel.run()?;

    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();

    // The turbine points
    layout_visualization::plot_turbine_points(
        &fmodel,
        format!("{output_dir}/001_turbine_points.png"),
        None,
        "k",
        10,
    )?;

    // The turbine points and their labels
    let turbine_names: Vec<String> = vec![
        "T10".to_string(),
        "T11".to_string(),
        "T12".to_string(),
        "T13".to_string(),
        "T22".to_string(),
    ];
    layout_visualization::plot_turbine_labels(
        &fmodel,
        format!("{output_dir}/001_turbine_labels.png"),
        Some(&turbine_names),
        None,
        false,
    )?;

    // The waking directions
    layout_visualization::plot_waking_directions(
        &fmodel,
        format!("{output_dir}/001_waking_directions.png"),
        None,
        Some(2),
    )?;

    // The flow field at hub height
    let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 100)?;
    flow_visualization::visualize_cut_plane(
        &horizontal_plane,
        format!("{output_dir}/001_horizontal_plane.png"),
        Some(1.0),
        Some(8.0),
        "coolwarm",
        false,
        "Horizontal plane at hub height",
    )?;

    println!("The visualizations were written to {output_dir}:");
    println!("  001_turbine_points.png");
    println!("  001_turbine_labels.png");
    println!("  001_waking_directions.png");
    println!("  001_horizontal_plane.png");

    Ok(())
}
