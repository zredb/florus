//! Example: Visualize Cross Plane (from examples_visualizations folder)
//!
//! Demonstrate visualizing a plane cut vertically through the flow field across the wind direction.
//! This corresponds to `examples_visualizations/003_visualize_cross_plane.py` in Python FLORIS.
//!
//! Note: Full cross-plane calculation requires implementing calculate_cross_plane() method.
//! This example uses horizontal plane as a placeholder to demonstrate the concept.

use florus::{FlorisModel, Result};
use florus::visualization::*;

fn main() -> Result<()> {
    println!("=== Cross Plane Visualization Example ===\n");
    println!("This example demonstrates cross-stream plane visualization concepts.");
    println!("Note: Full cross-plane implementation requires calculate_cross_plane() method.\n");

    // Initialize FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 1 turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0]),
        &ndarray::arr1(&[0.0]),
    )?;
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    // Create output directory
    std::fs::create_dir_all("output/examples_visualizations")?;

    println!("Running simulation...");
    fmodel.run()?;

    // Since calculate_cross_plane() is not yet implemented, we'll use horizontal plane
    // as a demonstration of flow visualization capabilities
    println!("\nNote: Full cross-plane visualization requires implementing:");
    println!("  - calculate_cross_plane() method in FlorisModel");
    println!("  - Vertical plane extraction perpendicular to wind direction");
    println!("  - Y-Z coordinate system at specified downstream distance\n");

    // Demonstrate with horizontal plane instead
    println!("Creating horizontal plane visualization as demonstration...");
    match visualize_horizontal_plane(
        &fmodel,
        "output/examples_visualizations/cross_plane_demo.png",
        Some(3.0),   // Min speed
        Some(9.0),   // Max speed
        Some("Horizontal Plane (Cross-plane not yet implemented)"),
    ) {
        Ok(_) => {
            println!("   ✓ Saved to output/examples_visualizations/cross_plane_demo.png");
            println!("   ℹ This shows horizontal plane at hub height");
            println!("   ℹ True cross-plane would show vertical slice at downstream_dist=500m");
        }
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    println!("\n=== Cross Plane Example Complete ===");
    println!("\nTo fully implement cross-plane visualization:");
    println!("  1. Add calculate_cross_plane() to FlorisModel");
    println!("  2. Extract vertical slice at specified downstream_dist");
    println!("  3. Interpolate flow field onto Y-Z grid perpendicular to wind");
    println!("  4. Visualize with appropriate axis labels and aspect ratio");

    Ok(())
}
