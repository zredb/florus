//! Example: Visualize Flow by Sweeping Turbines (from examples_visualizations folder)
//!
//! Demonstrate the use calculate_horizontal_plane_with_turbines.
//! This corresponds to `examples_visualizations/005_visualize_flow_by_sweeping_turbines.py` in Python FLORIS.
//!
//! Note: Full implementation requires calculate_horizontal_plane_with_turbines() method.
//! This example uses standard horizontal plane calculation as a placeholder.

use florus::{FlorisModel, Result};
use florus::visualization::*;

fn main() -> Result<()> {
    println!("=== Flow Visualization by Sweeping Turbines Example ===\n");
    println!("This example demonstrates turbine scanning for flow visualization.");
    println!("Note: Full implementation requires calculate_horizontal_plane_with_turbines().\n");

    // Initialize FLORIS
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set a 2 turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
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

    println!("\nNote: The 'calculate_horizontal_plane_with_turbines' method is not yet implemented.");
    println!("This method would:");
    println!("  - Scan across the flow field using turbine models");
    println!("  - Calculate velocities at each point by applying wake models");
    println!("  - Useful for wake models without native visualization support\n");

    // Use standard horizontal plane as demonstration
    println!("Creating standard horizontal plane visualization as demonstration...");
    match visualize_horizontal_plane(
        &fmodel,
        "output/examples_visualizations/sweep_turbines_demo.png",
        Some(3.0),   // Min speed
        Some(9.0),   // Max speed
        Some("Horizontal Plane (Turbine Sweep Method Not Implemented)"),
    ) {
        Ok(_) => {
            println!("   ✓ Saved to output/examples_visualizations/sweep_turbines_demo.png");
            println!("   ℹ This shows standard horizontal plane calculation");
            println!("   ℹ True turbine sweep would scan points using wake model directly");
        }
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    println!("\n=== Sweeping Turbines Example Complete ===");
    println!("\nTo fully implement turbine sweep visualization:");
    println!("  1. Add calculate_horizontal_plane_with_turbines() to flow_visualization module");
    println!("  2. Create grid of points across horizontal plane");
    println!("  3. For each point, calculate velocity by applying wake models");
    println!("  4. Return CutPlane object with interpolated values");
    println!("\nUse case:");
    println!("  - Some wake models may not have native visualization methods");
    println!("  - Turbine sweep provides a generic way to visualize any wake model");
    println!("  - Slower but more flexible than direct flow field extraction");

    Ok(())
}
