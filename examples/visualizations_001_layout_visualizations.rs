//! Example: Layout Visualizations (from examples_visualizations folder)
//!
//! Demonstrate the use of all the functions within the layout_visualization module.
//! This corresponds to `examples_visualizations/001_layout_visualizations.py` in Python FLORIS.
//!
//! Note: Since Plotters doesn't support multi-subplot layouts like matplotlib,
//! each visualization is saved as a separate file instead of subplots.

use florus::{FlorisModel, Result};
use florus::visualization::*;

fn main() -> Result<()> {
    println!("=== Layout Visualizations Example ===\n");
    println!("This example demonstrates all layout visualization functions.");
    println!("Each plot is saved as a separate file (Plotters limitation).\n");

    // Initialize FLORIS with the given input file
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Change to 5-turbine layout with a wind direction from northwest
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 0.0, 1000.0, 1000.0, 1000.0]),
        &ndarray::arr1(&[0.0, 500.0, 0.0, 500.0, 1000.0]),
    )?;
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[300.0]),
        ndarray::arr1(&[0.06]),
    )?;

    // Create output directory
    std::fs::create_dir_all("output/examples_visualizations")?;

    let min_ws = 1.0;
    let max_ws = 8.0;

    // Plot 1: Flow visualization and turbine points
    println!("1. Creating flow visualization with turbine points...");
    fmodel.run()?;
    match visualize_horizontal_plane(
        &fmodel,
        "output/examples_visualizations/plot_1_flow_and_points.png",
        Some(min_ws),
        Some(max_ws),
        Some("Flow visualization and turbine points"),
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_1_flow_and_points.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Plot 2: Show turbine names with a red bounding box
    println!("\n2. Creating turbine labels with bounding box...");
    let turbine_names = vec![
        "T10".to_string(), "T11".to_string(), "T12".to_string(),
        "T13".to_string(), "T22".to_string(),
    ];
    match plot_turbine_labels(
        &fmodel,
        "output/examples_visualizations/plot_2_labels_with_bbox.png",
        Some(&turbine_names),
        None,
        true,  // show_bbox
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_2_labels_with_bbox.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Plot 3: Flow visualization with yawed turbine
    println!("\n3. Creating flow visualization with yawed turbine...");
    fmodel.set_yaw_angles(ndarray::arr2(&[[0.0, 30.0, 0.0, 0.0, 0.0]]))?;
    fmodel.run()?;
    
    // First create the flow visualization
    match visualize_horizontal_plane(
        &fmodel,
        "output/examples_visualizations/plot_3_flow_with_yaw.png",
        Some(min_ws),
        Some(max_ws),
        Some("Flow visualization with yawed turbine"),
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_3_flow_with_yaw.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }
    
    // Then overlay rotor orientations
    match plot_turbine_rotors(
        &fmodel,
        "output/examples_visualizations/plot_3_rotors_overlay.png",
        "white",
        Some(300.0),
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_3_rotors_overlay.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Plot 4: Show turbine names and wake direction
    println!("\n4. Creating layout with wake directions...");
    fmodel.set_yaw_angles(ndarray::arr2(&[[0.0, 0.0, 0.0, 0.0, 0.0]]))?;
    match plot_waking_directions(
        &fmodel,
        "output/examples_visualizations/plot_4_wake_directions.png",
        None,  // No distance limit
        None,  // No connection limit
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_4_wake_directions.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Plot 5: Plot a subset and limit wake line distance
    println!("\n5. Creating subset layout with limited wake lines...");
    match plot_waking_directions(
        &fmodel,
        "output/examples_visualizations/plot_5_subset_limited.png",
        Some(7.0),  // Limit to 7D distance
        None,       // No connection limit
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_5_subset_limited.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Plot 6: Plot with a shaded region
    println!("\n6. Creating layout with shaded region...");
    let region_points = vec![
        (0.0, 0.0),
        (300.0, 0.0),
        (300.0, 1000.0),
        (0.0, 700.0),
    ];
    match shade_region(
        &region_points,
        "output/examples_visualizations/plot_6_shaded_region.png",
        false,      // Don't show vertex points
        "black",    // Region color
        0.3,        // Alpha (transparency)
        "black",    // Point color (not used since show_points=false)
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_6_shaded_region.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Plot 7: Plot farm terrain (hub heights as proxy for terrain)
    println!("\n7. Creating farm terrain visualization...");
    // Note: In Rust version, we need to modify hub_heights through the core
    // For this example, we'll just use the existing hub heights
    match plot_farm_terrain(
        &fmodel,
        "output/examples_visualizations/plot_7_farm_terrain.png",
    ) {
        Ok(_) => println!("   ✓ Saved to output/examples_visualizations/plot_7_farm_terrain.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    println!("\n=== Layout Visualizations Complete ===");
    println!("\nGenerated files:");
    println!("  1. plot_1_flow_and_points.png - Flow field with turbine points");
    println!("  2. plot_2_labels_with_bbox.png - Turbine labels with bounding box");
    println!("  3. plot_3_flow_with_yaw.png - Flow field with yawed turbine");
    println!("     plot_3_rotors_overlay.png - Rotor orientations (separate file)");
    println!("  4. plot_4_wake_directions.png - Wake directions between turbines");
    println!("  5. plot_5_subset_limited.png - Subset with limited wake lines");
    println!("  6. plot_6_shaded_region.png - Shaded polygonal region");
    println!("  7. plot_7_farm_terrain.png - Farm terrain (hub heights)");
    println!("\nNote: Unlike Python's matplotlib subplot system, Plotters saves");
    println!("each visualization as a separate file.");

    Ok(())
}
