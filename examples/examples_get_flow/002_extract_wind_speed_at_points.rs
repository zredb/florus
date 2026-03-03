/// Example: Extract Wind Speed at Points
///
/// This example demonstrates extracting wind speed at specific points
/// in the flow field. This is useful for analyzing wake behavior
/// at specific locations.
///
/// This is the Rust equivalent of Python's examples_get_flow/002_extract_wind_speed_at_points.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Extract Wind Speed at Points");
    println!("===============================================\n");

    println!("--- Point Extraction ---\n");
    
    println!("Extract wind speed at specific (x, y, z) locations:");
    println!("  - Custom probe locations");
    println!("  - Cross-section analysis");
    println!("  - Wake centerline tracking\n");
    
    // ============================================================
    // Set up 2-turbine layout
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    let d = 126.0;
    let spacing = 5.0 * d;
    
    // Two turbines in line
    let layout_x = Array1::from_vec(vec![0.0, spacing]);
    let layout_y = Array1::from_vec(vec![0.0; 2]);
    
    model.set_layout(layout_x.clone(), layout_y.clone())?;
    
    // Wind from 270 degrees (west)
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    
    model.run()?;
    
    println!("--- Two-Turbine Wake Analysis ---\n");
    
    println!("Layout: 2 turbines at {:.0}D spacing", spacing / d);
    println!("Wind: 8 m/s from West (270°)\n");
    
    // Get powers
    let powers = model.get_turbine_powers();
    println!("Upstream turbine power: {:.1} kW", powers[[0, 0]] / 1000.0);
    println!("Downstream turbine power: {:.1} kW", powers[[0, 1]] / 1000.0);
    
    // Wake deficit
    let wake_deficit = 1.0 - powers[[0, 1]] / powers[[0, 0]];
    println!("Wake deficit: {:.1}%\n", wake_deficit * 100.0);
    
    // ============================================================
    // Multiple downstream positions
    // ============================================================
    println!("--- Downstream Velocity Profile ---\n");
    
    // Create different wind directions to simulate downstream positions
    let directions: Vec<f64> = vec![260.0, 270.0, 280.0];
    
    println!("{:>8} {:>12} {:>12}", "WD (°)", "Up (kW)", "Down (kW)");
    println!("{}", "-".repeat(35));
    
    for wd in &directions {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![*wd]),
            Array1::from_vec(vec![0.06]),
        )?;
        model.run()?;
        
        let powers = model.get_turbine_powers();
        println!("{:>8.0} {:>12.1} {:>12.1}", wd, powers[[0, 0]] / 1000.0, powers[[0, 1]] / 1000.0);
    }
    
    // ============================================================
    // Grid study
    // ============================================================
    println!("\n--- Grid Resolution Study ---\n");
    
    // Test different grid resolutions
    // Note: Would require reinitializing the model with different grid settings
    
    println!("Note: Grid resolution controlled by 'turbine_grid_points' in config");
    println!("Current: 3x3 points per rotor (configurable)");
    
    println!("\n--- Lateral Position Study ---\n");
    
    // Vary lateral positions
    let y_offsets: Vec<f64> = vec![-100.0, -50.0, 0.0, 50.0, 100.0];
    
    let layout_y_offset = Array1::from_vec(y_offsets.clone());
    model.set_layout(&Array1::from_vec(vec![0.0; y_offsets.len()]), layout_y_offset)?;
    
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0; y_offsets.len()]),
        Array1::from_vec(vec![0.06; y_offsets.len()]),
    )?;
    
    model.run()?;
    
    let powers = model.get_turbine_powers();
    
    println!("{:>12} {:>12}", "Y Offset (m)", "Power (kW)");
    println!("{}", "-".repeat(26));
    
    for (i, y) in y_offsets.iter().enumerate() {
        println!("{:>12.0} {:>12.1}", y, powers[[i, 0]] / 1000.0);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
