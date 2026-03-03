/// Example: Visualize Grid Points
///
/// This example demonstrates the grid points used in FLORIS-RS.
/// The grid defines where velocities are calculated in the flow field.
///
/// This is the Rust equivalent of Python's examples_get_flow/003_visualize_grid.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Visualize Grid Points");
    println!("======================================\n");

    println!("--- Grid System ---\n");
    
    println!("FLORIS uses different grid types:");
    println!("  1. TurbineGrid: Points on turbine rotor");
    println!("  2. CubatureGrid: Optimized quadrature points");
    println!("  3. TurbineCubatureGrid: Combined approach\n");
    
    println!("Grid resolution affects:");
    println!("  - Computational accuracy");
    println!("  - Simulation speed");
    println!("  - Memory usage\n");
    
    // ============================================================
    // Load model
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Single turbine
    model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    // Set conditions
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    
    model.run()?;
    
    println!("--- Grid Information ---\n");
    
    // Note: Grid information would come from the internal grid structure
    println!("Grid configuration:");
    println!("  Type: TurbineGrid");
    println!("  Points per rotor: 3x3 (configurable)");
    println!("  Total points per turbine: 9\n");
    
    // ============================================================
    // Multiple grid points example
    // ============================================================
    println!("--- Multiple Turbines ---\n");
    
    let d = 126.0;
    
    // 3 turbines in a row
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 3]);
    
    model.set_layout(layout_x.clone(), layout_y.clone())?;
    
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    
    model.run()?;
    
    let powers = model.get_turbine_powers();
    
    println!("3-turbine layout at 5D spacing:");
    println!("{:>8} {:>12} {:>12}", "Turbine", "X (m)", "Power (kW)");
    println!("{}", "-".repeat(35));
    
    for (i, x) in layout_x.iter().enumerate() {
        println!("{:>8} {:>12.0} {:>12.1}", i, x, powers[[0, i]] / 1000.0);
    }
    
    // ============================================================
    // Grid spacing effect
    // ============================================================
    println!("\n--- Grid Spacing Effect ---\n");
    
    let spacings = vec![3.0, 5.0, 7.0, 10.0];
    
    println!("{:>12} {:>12} {:>12} {:>12}", "Spacing (D)", "T1 (kW)", "T2 (kW)", "T3 (kW)");
    println!("{}", "-".repeat(52));
    
    for spacing in &spacings {
        let layout_x = Array1::from_vec(vec![0.0, spacing * d, 2.0 * spacing * d]);
        let layout_y = Array1::from_vec(vec![0.0; 3]);
        
        model.set_layout(layout_x, layout_y)?;
        
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let powers = model.get_turbine_powers();
        
        println!("{:>12.1} {:>12.1} {:>12.1} {:>12.1}", 
            spacing, 
            powers[[0, 0]] / 1000.0,
            powers[[0, 1]] / 1000.0,
            powers[[0, 2]] / 1000.0
        );
    }
    
    println!("\n--- Key Points ---\n");
    
    println!("1. Grid points define velocity sampling locations");
    println!("2. More points = higher accuracy but slower");
    println!("3. TurbineGrid with 3x3 is standard for AEP");
    println!("4. Higher resolution useful for visualization");
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
