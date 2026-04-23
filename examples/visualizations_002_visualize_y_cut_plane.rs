//! Example: Visualize Y Cut Plane (from examples_visualizations folder)
//!
//! Demonstrate visualizing a plane cut vertically through the flow field along the wind direction.
//! This corresponds to `examples_visualizations/002_visualize_y_cut_plane.py` in Python FLORIS.
//!
//! Note: Full Y-plane calculation requires implementing calculate_y_plane() method.
//! This example uses horizontal plane as a placeholder to demonstrate the concept.

use florus::{FlorisModel, Result};
use florus::visualization::*;

fn main() -> Result<()> {
    println!("=== Y Cut Plane Visualization Example ===\n");
    println!("This example demonstrates vertical plane visualization concepts.");
    println!("Note: Full Y-plane implementation requires calculate_y_plane() method.\n");

    // Initialize FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 3 turbine layout with wind direction along the row
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0, 1000.0]),
        &ndarray::arr1(&[0.0, 0.0, 0.0]),
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

    // Since calculate_y_plane() is not yet implemented, we'll use horizontal plane
    // as a demonstration of flow visualization capabilities
    println!("\nNote: Full Y-plane visualization requires implementating:");
    println!("  - calculate_y_plane() method in FlorisModel");
    println!("  - Vertical plane extraction from 3D flow field");
    println!("  - Y-Z coordinate system visualization\n");

    // Demonstrate with horizontal plane instead
    println!("Creating horizontal plane visualization as demonstration...");
    match visualize_horizontal_plane(
        &fmodel,
        "output/examples_visualizations/y_plane_demo.png",
        Some(3.0),   // Min speed
        Some(9.0),   // Max speed
        Some("Horizontal Plane (Y-plane not yet implemented)"),
    ) {
        Ok(_) => {
            println!("   ✓ Saved to output/examples_visualizations/y_plane_demo.png");
            println!("   ℹ This shows horizontal plane at hub height");
            println!("   ℹ True Y-plane would show vertical slice along wind direction");
        }
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    println!("\n=== Y Cut Plane Example Complete ===");
    println!("\nTo fully implement Y-plane visualization:");
    println!("  1. Add calculate_y_plane() to FlorisModel");
    println!("  2. Extract vertical slice at specified crossstream_dist");
    println!("  3. Interpolate flow field onto Y-Z grid");
    println!("  4. Visualize with appropriate axis labels");

    Ok(())
}
