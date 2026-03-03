/// Example: Compare Farm Power with Neighbor
///
/// This example demonstrates comparing power production between different
/// wind farm configurations or with neighboring farms.
///
/// This is the Rust equivalent of Python's 010_compare_farm_power_with_neighbor.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 010: Compare Farm Power with Neighbor");
    println!("====================================================\n");

    println!("--- Farm Comparison ---\n");
    
    println!("This example demonstrates:");
    println!("  - Single farm vs multi-farm comparison");
    println!("  - Layout optimization impact");
    println!("  - Neighbor farm wake effects\n");
    
    // ============================================================
    // Single farm baseline
    // ============================================================
    println!("--- Single Farm Baseline ---\n");
    
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    let d = 126.0;
    let spacing = 5.0 * d;
    
    // 4 turbines in a row
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing, 3.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0; 4]);
    
    model.set_layout(layout_x.clone(), layout_y.clone())?;
    
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    
    model.run()?;
    
    let single_farm_power = model.get_farm_power()[[0]] / 1000.0;
    
    println!("Single farm (4 turbines at 5D spacing):");
    println!("  Total power: {:.1} kW\n", single_farm_power);
    
    // ============================================================
    // Spacing comparison
    // ============================================================
    println!("--- Spacing Comparison ---\n");
    
    let spacings = vec![3.0, 5.0, 7.0, 10.0, 15.0];
    
    println!("{:>12} {:>12} {:>12}", "Spacing (D)", "Farm (kW)", "Per Turb (kW)");
    println!("{}", "-".repeat(40));
    
    for sp in &spacings {
        let layout_x = Array1::from_vec(vec![
            0.0, 
            sp * d, 
            2.0 * sp * d, 
            3.0 * sp * d
        ]);
        let layout_y = Array1::from_vec(vec![0.0; 4]);
        
        model.set_layout(layout_x, layout_y)?;
        
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let farm_power = model.get_farm_power()[[0]] / 1000.0;
        let per_turb = farm_power / 4.0;
        
        println!("{:>12.1} {:>12.1} {:>12.1}", sp, farm_power, per_turb);
    }
    
    // ============================================================
    // Wind direction comparison
    // ============================================================
    println!("\n--- Wind Direction Comparison ---\n");
    
    let directions: Vec<f64> = (180..340).step_by(20).map(|d| d as f64).collect();
    
    // Reset to 5D spacing
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing, 3.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0; 4]);
    
    model.set_layout(layout_x.clone(), layout_y.clone())?;
    
    println!("{:>8} {:>12} {:>12}", "WD (°)", "Farm (kW)", "Loss %");
    println!("{}", "-".repeat(35));
    
    // Calculate freestream power (single turbine)
    model.set_layout(
        Array1::from_vec(vec![0.0]),
        Array1::from_vec(vec![0.0]),
    )?;
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    model.run()?;
    let freestream = model.get_turbine_powers()[[0, 0]] / 1000.0;
    let max_power = freestream * 4.0;
    
    // Restore 4-turbine layout
    model.set_layout(layout_x, layout_y)?;
    
    for wd in &directions {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![*wd]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let farm_power = model.get_farm_power()[[0]] / 1000.0;
        let loss = (1.0 - farm_power / max_power) * 100.0;
        
        println!("{:>8.0} {:>12.1} {:>12.1}%", wd, farm_power, loss);
    }
    
    // ============================================================
    // Layout pattern comparison
    // ============================================================
    println!("\n--- Layout Pattern Comparison ---\n");
    
    // Linear layout
    let linear_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing, 3.0 * spacing]);
    let linear_y = Array1::from_vec(vec![0.0; 4]);
    
    // Staggered layout
    let stagger_x = Array1::from_vec(vec![0.0, spacing, 0.0, spacing]);
    let stagger_y = Array1::from_vec(vec![0.0, 0.0, spacing, spacing]);
    
    // Grid layout (2x2)
    let grid_x = Array1::from_vec(vec![0.0, spacing, 0.0, spacing]);
    let grid_y = Array1::from_vec(vec![0.0, 0.0, spacing, spacing]);
    
    let patterns = vec![
        ("Linear", linear_x.clone(), linear_y.clone()),
        ("Staggered", stagger_x.clone(), stagger_y.clone()),
        ("Grid", grid_x.clone(), grid_y.clone()),
    ];
    
    println!("4-turbine layouts at 5D spacing:\n");
    println!("{:>12} {:>12}", "Pattern", "Farm (kW)");
    println!("{}", "-".repeat(26));
    
    for (name, x, y) in &patterns {
        model.set_layout(x.clone(), y.clone())?;
        
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let power = model.get_farm_power()[[0]] / 1000.0;
        println!("{:>12} {:>12.1}", name, power);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
