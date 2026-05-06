//! Example: Y-Plane and Cross-Plane Visualization
//!
//! Demonstrate visualizing y-plane and cross-plane cuts through the flow field.

use florus::{FlorisModel, Result};
use florus::visualization::flow_visualization;
use ndarray::Array;
use std::path::Path;

fn main() -> Result<()> {
    println!("=== FLORUS Y-Plane and Cross-Plane Visualization ===\n");

    // Initialize FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 3-turbine layout with wind direction along the row
    let layout_x = Array::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = Array::from_vec(vec![0.0, 0.0, 0.0]);
    fmodel.set_layout(&layout_x, &layout_y)?;
    
    let wind_speeds = Array::from_vec(vec![8.0]);
    let wind_directions = Array::from_vec(vec![270.0]);
    let turbulence_intensities = Array::from_vec(vec![0.06]);
    fmodel.set(
        Some(wind_speeds),
        Some(wind_directions),
        None,  // wind_shear
        None,  // wind_veer
        None,  // reference_wind_height
        Some(turbulence_intensities),
        None,  // air_density
        None,  // layout_x
        None,  // layout_y
        None,  // yaw_angles
        None,  // power_setpoints
        None,  // awc_modes
        None,  // awc_amplitudes
        None,  // awc_frequencies
        None,  // disable_turbines
    )?;

    // Create output directory
    let output_dir = "examples/outputs/visualization";
    std::fs::create_dir_all(output_dir)?;

    // Example 1: Y-plane visualization
    println!("Example 1: Y-plane visualization");
    visualize_y_plane(&fmodel, &format!("{}/04_y_plane.png", output_dir))?;

    // Example 2: Cross-plane visualization
    println!("\nExample 2: Cross-plane visualization");
    visualize_cross_plane(&fmodel, &format!("{}/05_cross_plane.png", output_dir))?;

    println!("\n=== Visualization complete! ===");
    println!("Output files saved to: {}", output_dir);

    Ok(())
}

/// Visualize y-plane (vertical plane along wind direction)
fn visualize_y_plane<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
) -> Result<()> {
    // Calculate y-plane at crossstream_dist=0.0
    let y_plane = fmodel.calculate_y_plane(200, 100, 0.0)?;

    // Visualize the cut plane
    flow_visualization::visualize_cut_plane(
        &y_plane,
        &output_path,
        Some(3.0),  // min_speed
        Some(9.0),  // max_speed
        "coolwarm",
        false,      // color_bar
        "Y Cut Plane",
    )?;

    println!("  Saved: {}", output_path.as_ref().display());
    Ok(())
}

/// Visualize cross-plane (vertical plane across wind direction)
fn visualize_cross_plane<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
) -> Result<()> {
    // Calculate cross-plane at downstream_dist=500.0
    let cross_plane = fmodel.calculate_cross_plane(100, 100, 500.0)?;

    // Visualize the cut plane
    flow_visualization::visualize_cut_plane(
        &cross_plane,
        &output_path,
        Some(3.0),  // min_speed
        Some(9.0),  // max_speed
        "coolwarm",
        false,      // color_bar
        "Cross Plane",
    )?;

    println!("  Saved: {}", output_path.as_ref().display());
    Ok(())
}
