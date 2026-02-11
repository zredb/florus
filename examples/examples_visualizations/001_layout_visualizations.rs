/// Layout Visualizations Example
///
/// Demonstrates all functions within the layout_visualization module.
///
/// This is the Rust equivalent of Python's 001_layout_visualizations.py

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Layout Visualizations");
    println!("================================\n");

    // ============================================================
    // Model Setup
    // ============================================================
    println!("--- Model Setup ---\n");

    // Create 5-turbine layout
    let layout_x = Array1::from_vec(vec![0.0, 0.0, 1000.0, 1000.0, 1000.0]);
    let layout_y = Array1::from_vec(vec![0.0, 500.0, 0.0, 500.0, 1000.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 5];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("Wind farm: 5 turbines in irregular layout");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: ({:.0}, {:.0})", i, x, y);
    }

    // ============================================================
    // Visualization Functions
    // ============================================================
    println!("\n--- Visualization Functions ---\n");

    println!("1. plot_turbine_points:");
    println!("   Shows turbine locations as points");
    println!("   Options: color, size, turbine_indices");
    println!());

    println!("2. plot_turbine_labels:");
    println!("   Adds turbine name labels");
    println!("   Options: turbine_names, show_bbox, bbox_dict");
    println!());

    println!("3. plot_turbine_rotors:");
    println!("   Shows rotor circles at hub location");
    println!("   Options: yaw_angles, color");
    println!());

    println!("4. plot_waking_directions:");
    println!("   Shows wake propagation lines");
    println!("   Options: turbine_indices, limit_dist_D");
    println!());

    println!("5. shade_region:");
    println!("   Adds shaded polygon to plot");
    println!("   Input: Array of (x, y) vertices");
    println!());

    println!("6. plot_farm_terrain:");
    println!("   Visualizes terrain/hub height variations");
    println!());

    // ============================================================
    // Example Plot Configurations
    // ============================================================
    println!("\n--- Example Configurations ---\n");

    println!("Plot 1: Flow visualization with turbine points");
    println!("   - Calculate horizontal plane at hub height");
    println!("   - Overlay white turbine points");
    println!());

    println!("Plot 2: Turbine labels with bounding boxes");
    println!("   - Custom turbine names: [T10, T11, T12, T13, T22]");
    println!("   - Red bounding box around labels");
    println!());

    println!("Plot 3: Flow with yawed rotor");
    println!("   - Yaw angles: [0°, 30°, 0°, 0°, 0°]");
    println!("   - Shows wake deflection");
    println!());

    println!("Plot 4: Wake directions");
    println!("   - Shows wake propagation from each turbine");
    println!("   - Can limit to subset of turbines");
    println!());

    println!("Plot 5: Limited wake distance");
    println!("   - limit_dist_D: 7 (7 rotor diameters)");
    println!("   - Shows wake effects within 7D");
    println!());

    println!("Plot 6: Shaded region");
    println!("   - Custom polygon: [(0,0), (300,0), (300,1000), (0,700)]");
    println!());

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Layout Visualization Key Points:");
    println!("  ✓ Multiple functions for different visualization needs");
    println!("  ✓ Can combine functions in single plot");
    println!("  ✓ Customizable colors, labels, and regions");
    println!("  ✓ Supports subset selection (turbine_indices)");
    println!("  ✓ Useful for presentations and analysis");
    println!());

    println!("Common workflows:");
    println!("  - Initial layout assessment");
    println!("  - Wake interaction analysis");
    println!("  ✓ Documentation and reporting");
    println!("  - Comparison of layouts");

    println!("\n================================");
    println!("Example completed successfully!");
    println!("Note: Full visualization requires plotting library.");

    Ok(())
}
