//! Example: Flow Visualization
//!
//! Demonstrate the use of flow visualization functions in FLORUS
//! Corresponds to examples_visualizations in Python FLORIS

use florus::{FlorisModel, Result};
use florus::visualization::flow_visualization;
use ndarray::Array;
use std::path::Path;

fn main() -> Result<()> {
    println!("=== FLORUS Flow Visualization Example ===\n");

    // Initialize FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 5-turbine layout with wind direction from northwest
    let layout_x = Array::from_vec(vec![0.0, 630.0, 1260.0]);
    let layout_y = Array::from_vec(vec![0.0, 0.0, 0.0]);
    fmodel.set_layout(&layout_x, &layout_y)?;
    
    let wind_speeds = Array::from_vec(vec![8.0]);
    let wind_directions = Array::from_vec(vec![300.0]);
    fmodel.set(
        Some(wind_speeds),
        Some(wind_directions),
        None,  // reference_wind_height
        None,  // turbulence_intensities (use default)
        None,  // air_density
        None,  // wind_shear
        None,  // wind_veer
        None,  // turbine_layout
        None,  // yaw_angles
        None,  // tilt_angles
        None,  // power_setpoints
        None,  // awc_modes
        None,  // awc_amplitudes
        None,  // awc_frequencies
        None,  // correct_cp_ct_for_tilt
    )?;

    // Create output directory
    let output_dir = "examples/outputs/visualization";
    std::fs::create_dir_all(output_dir)?;

    // Example 1: Visualize horizontal plane with turbine points
    println!("Example 1: Flow visualization with turbine points");
    visualize_horizontal_plane_with_turbines(&fmodel, &format!("{}/01_horizontal_flow.png", output_dir))?;

    // Example 2: Visualize with yawed turbine
    println!("\nExample 2: Flow visualization with yawed turbine");
    visualize_yawed_flow(&mut fmodel, &format!("{}/02_yawed_flow.png", output_dir))?;

    // Example 3: Visualize with turbine rotors and labels
    println!("\nExample 3: Flow visualization with turbine rotors and labels");
    visualize_flow_with_rotors_and_labels(&fmodel, &format!("{}/03_flow_with_rotors_labels.png", output_dir))?;

    println!("\n=== Visualization complete! ===");
    println!("Output files saved to: {}", output_dir);

    Ok(())
}

/// Visualize horizontal plane with turbine points
fn visualize_horizontal_plane_with_turbines<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
) -> Result<()> {
    // Calculate horizontal plane
    let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 200)?;

    // Visualize the cut plane
    flow_visualization::visualize_cut_plane(
        &horizontal_plane,
        &output_path,
        Some(1.0),  // min_speed
        Some(8.0),  // max_speed
        "coolwarm",
        false,      // color_bar
        "",         // title
    )?;

    println!("  Saved: {}", output_path.as_ref().display());
    Ok(())
}

/// Visualize flow with yawed turbine
fn visualize_yawed_flow<P: AsRef<Path>>(
    fmodel: &mut FlorisModel,
    output_path: P,
) -> Result<()> {
    // Set yaw angles for second turbine
    let yaw_angles = Array::from_shape_vec((1, 3), vec![0.0, 30.0, 0.0])?;
    fmodel.set_yaw_angles(yaw_angles)?;

    // Calculate horizontal plane
    let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 200)?;

    // Visualize with turbine rotors showing yaw
    flow_visualization::visualize_cut_plane_with_rotors(
        &horizontal_plane,
        fmodel,
        &output_path,
        Some(1.0),  // min_speed
        Some(8.0),  // max_speed
        "coolwarm",
        true,       // show_rotors
        false,      // color_bar
        "Flow visualization with yawed turbine",
    )?;

    // Reset yaw angles
    let yaw_angles_reset = Array::from_shape_vec((1, 3), vec![0.0, 0.0, 0.0])?;
    fmodel.set_yaw_angles(yaw_angles_reset)?;

    println!("  Saved: {}", output_path.as_ref().display());
    Ok(())
}

/// Visualize flow with turbine rotors and labels
fn visualize_flow_with_rotors_and_labels<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
) -> Result<()> {
    // Set turbine names
    let turbine_names = vec!["T1", "T2", "T3"];

    // Calculate horizontal plane
    let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 200)?;

    // Visualize with rotors and labels
    flow_visualization::visualize_cut_plane_with_rotors_and_labels(
        &horizontal_plane,
        fmodel,
        &output_path,
        Some(1.0),  // min_speed
        Some(8.0),  // max_speed
        "coolwarm",
        &turbine_names,
        true,       // show_rotors
        false,      // color_bar
        "Horizontal Flow with Turbine Rotors and labels",
    )?;

    println!("  Saved: {}", output_path.as_ref().display());
    Ok(())
}
