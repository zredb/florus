//! Example 002: Visualizations
//!
//! This example demonstrates the use of the flow and layout visualizations in FLORUS.
//! First, an example wind farm layout is plotted, with the turbine names and the directions
//! and distances between turbines shown in different configurations by subplot.
//! Next, the horizontal flow field at hub height is plotted for a single wind condition.
//!
//! FLORUS includes two modules for visualization:
//!   1) flow_visualization: for visualizing the flow field
//!   2) layout_visualization: for visualizing the layout of the wind farm
//! The two modules can be used together to visualize the flow field and the layout
//! of the wind farm.
//!
//! Corresponds to `002_visualizations.py` in Python FLORIS.

use florus::{FlorisModel, Result};
use florus::visualization::*;

fn main() -> Result<()> {
    println!("=== FLORUS Visualization Example ===\n");

    // Load a FlorisModel
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set the farm layout to have 8 turbines irregularly placed
    let layout_x = vec![0.0, 500.0, 0.0, 128.0, 1000.0, 900.0, 1500.0, 1250.0];
    let layout_y = vec![0.0, 300.0, 750.0, 1400.0, 0.0, 567.0, 888.0, 1450.0];
    fmodel.set_layout(
        &ndarray::arr1(&layout_x),
        &ndarray::arr1(&layout_y),
    )?;

    println!("Layout visualization contains functions for visualizing the layout:");
    println!("  - plot_turbine_points");
    println!("  - plot_turbine_labels");
    println!("  - plot_turbine_rotors");
    println!("  - plot_waking_directions");
    println!("Each can be overlaid to provide further information about the layout.");
    println!("\nCreating layout visualization subplots...\n");

    // Create output directory
    std::fs::create_dir_all("output")?;

    // Subplot 1: Turbine Points
    println!("1. Creating turbine points plot...");
    match plot_turbine_points(
        &fmodel,
        "output/layout_subplot_1_points.png",
        None,        // All turbines
        "black",     // Color
        10,          // Marker size
    ) {
        Ok(_) => println!("   ✓ Saved to output/layout_subplot_1_points.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Subplot 2: Turbine Points and Labels
    println!("\n2. Creating turbine points and labels plot...");
    match plot_turbine_labels(
        &fmodel,
        "output/layout_subplot_2_labeled.png",
        None,        // Use default indices
        None,        // Use default offset
        false,       // No bounding box
    ) {
        Ok(_) => println!("   ✓ Saved to output/layout_subplot_2_labeled.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Subplot 3: Turbine Points, Labels, and Waking Directions
    println!("\n3. Creating turbine points, labels, and waking directions plot...");
    match plot_waking_directions(
        &fmodel,
        "output/layout_subplot_3_waking.png",
        None,        // No distance limit
        Some(2),     // Limit to 2 nearest neighbors
    ) {
        Ok(_) => println!("   ✓ Saved to output/layout_subplot_3_waking.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Subplot 4: Use Provided Turbine Names
    println!("\n4. Creating plot with custom turbine names...");
    let turbine_names = vec![
        "T1".to_string(), "T2".to_string(), "T3".to_string(), "T4".to_string(),
        "T9".to_string(), "T10".to_string(), "T75".to_string(), "T78".to_string(),
    ];
    match plot_turbine_labels(
        &fmodel,
        "output/layout_subplot_4_custom_names.png",
        Some(&turbine_names),  // Use custom names
        None,                  // Use default offset
        false,                 // No bounding box
    ) {
        Ok(_) => println!("   ✓ Saved to output/layout_subplot_4_custom_names.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Visualizations of the flow field
    println!("\n--- Flow Field Visualization ---");
    println!("Visualizations of the flow field are made by using calculate plane methods.");
    println!("In this example we show the horizontal plane at hub height.\n");

    // For flow visualizations, set to a single condition (n_findex = 1)
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),     // Wind speed
        ndarray::arr1(&[290.0]),   // Wind direction
        ndarray::arr1(&[0.06]),    // Turbulence intensity
    )?;

    println!("Running simulation for flow visualization...");
    fmodel.run()?;

    // Plot the flow field
    println!("\n5. Creating horizontal flow field visualization...");
    match visualize_horizontal_plane(
        &fmodel,
        "output/horizontal_flow_with_rotors.png",
        Some(3.0),   // Min wind speed for color scale
        Some(12.0),  // Max wind speed for color scale
        Some("Horizontal Flow with Turbine Rotors and Labels"),
    ) {
        Ok(_) => println!("   ✓ Saved to output/horizontal_flow_with_rotors.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    // Plot turbine rotors with yaw angles
    println!("\n6. Setting yaw angles and plotting rotors on flow field...");
    fmodel.set_yaw_angles(ndarray::arr2(&[[0.0, 10.0, 5.0, 15.0, 0.0, 20.0, 10.0, 5.0]]))?;
    
    match plot_turbine_rotors(
        &fmodel,
        "output/flow_field_with_rotors.png",
        "black",
        Some(290.0),  // Use current wind direction
    ) {
        Ok(_) => println!("   ✓ Saved to output/flow_field_with_rotors.png"),
        Err(e) => eprintln!("   ✗ Error: {}", e),
    }

    println!("\n=== Visualization Complete ===");
    println!("All visualizations saved to 'output/' directory");
    println!("\nGenerated files:");
    println!("  Layout visualizations:");
    println!("    - layout_subplot_1_points.png");
    println!("    - layout_subplot_2_labeled.png");
    println!("    - layout_subplot_3_waking.png");
    println!("    - layout_subplot_4_custom_names.png");
    println!("  Flow field visualizations:");
    println!("    - horizontal_flow_with_rotors.png");
    println!("    - flow_field_with_rotors.png");

    Ok(())
}
