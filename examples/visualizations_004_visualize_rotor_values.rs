//! Example: Visualize the values at the rotors
//!
//! This example plots the velocity values at the resolution of the rotor
//! grids of the turbines.
//!
//! Scenario aligned 1:1 with
//! examples_visualizations/004_visualize_rotor_values.py.

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

    // The values at the rotor of every turbine
    let output_dir = "output/examples_visualizations";
    std::fs::create_dir_all(output_dir).ok();
    flow_visualization::plot_rotor_values(
        &fmodel.flow_field().u,
        0,
        1,
        2,
        format!("{output_dir}/004_rotor_values.png"),
        "viridis",
    )?;

    // The velocity range at the rotors
    let u = &fmodel.flow_field().u;
    let min = u.iter().copied().fold(f64::INFINITY, f64::min);
    let max = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("Rotor velocity range: {min:.3} .. {max:.3} m/s");

    Ok(())
}
