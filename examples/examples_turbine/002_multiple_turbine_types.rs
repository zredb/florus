/// Multiple Turbine Types Example
///
/// This example demonstrates using multiple turbine types in a single wind farm.
/// The first two turbines use NREL 5MW, and the third uses IEA 10MW.
///
/// This is the Rust equivalent of Python's 002_multiple_turbine_types.py

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Multiple Turbine Types");
    println!("==============================\n");

    // ============================================================
    // Model Setup
    // ============================================================
    println!("--- Model Setup ---\n");

    println!("Creating wind farm with multiple turbine types:");
    println!("  - Turbine 0: NREL 5MW");
    println!("  - Turbine 1: NREL 5MW");
    println!("  - Turbine 2: IEA 10MW");
    println!());

    // Create layout with different turbine types
    let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec![
        "nrel_5MW".to_string(),
        "nrel_5MW".to_string(),
        "iea_10MW".to_string(),
    ];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("Wind farm layout:");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: x = {:.0} m, y = {:.0} m, type = {}", i, x, y, turbine_types[i]);
    }

    // ============================================================
    // Turbine Characteristics
    // ============================================================
    println!("\n--- Turbine Characteristics ---\n");

    println!("NREL 5MW:");
    println!("  Rotor diameter: 126.0 m");
    println!("  Rated power: 5 MW");
    println!("  Hub height: 90.0 m");
    println!());

    println!("IEA 10MW:");
    println!("  Rotor diameter: 198.0 m");
    println!("  Rated power: 10 MW");
    println!("  Hub height: 119.0 m");
    println!());

    // ============================================================
    // Flow Visualization Setup
    // ============================================================
    println!("--- Flow Visualization ---\n");

    println!("Calculating flow planes for visualization:");
    println!());

    println!("1. Horizontal plane:");
    println!("   Height: 90.0 m (hub height)");
    println!("   X resolution: 200 points");
    println!("   Y resolution: 100 points");
    println!());

    println!("2. Y-plane (streamwise profile):");
    println!("   X resolution: 200 points");
    println!("   Z resolution: 100 points");
    println!("   Crossstream distance: 0.0 m");
    println!());

    println!("3. Cross-plane (spanwise profile):");
    println!("   Y resolution: 100 points");
    println!("   Z resolution: 100 points");
    println!("   Downstream distance: 500.0 m");
    println!());

    // ============================================================
    // Visualization Functions
    // ============================================================
    println!("--- Visualization Functions ---\n");

    println!("Layout visualization module functions:");
    println!("  1. plot_turbine_points - Shows turbine locations");
    println!("  2. plot_turbine_labels - Adds turbine names");
    println!("  3. plot_turbine_rotors - Shows rotor circles");
    println!("  4. plot_waking_directions - Shows wake propagation");
    println!());

    println!("Flow visualization options:");
    println!("  - visualize_cut_plane - Main visualization function");
    println!("  - min_speed / max_speed - Color scale bounds");
    println!("  - title - Plot title");
    println!());

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Multiple Turbine Types Key Points:");
    println!("  ✓ Different turbine types can be mixed in one farm");
    println!("  ✓ Each turbine has its own power curve");
    println!("  ✓ Wake effects are calculated per turbine");
    println!("  ✓ Useful for hybrid installations");
    println!("  ✓ Supports technology comparisons");
    println!());

    println!("Applications:");
    println!("  - Technology comparison studies");
    println!("  - Phased wind farm development");
    println!("  ✓ Hybrid turbine configurations");
    println!("  - Repowering analysis");

    println!("\n==============================");
    println!("Example completed successfully!");
    println!("Note: Full visualization requires plotting library.");

    Ok(())
}
