/// Example: Heterogeneous Inflow in 2D and 3D
///
/// This example demonstrates 2D and 3D heterogeneous inflow in FLORIS-RS.
/// 2D: Speed varies in horizontal (x) direction
/// 3D: Speed varies in vertical (z) direction
///
/// This is the Rust equivalent of Python's examples_heterogeneous/004_heterogeneous_2d_and_3d.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Heterogeneous Inflow in 2D and 3D");
    println!("==================================================\n");

    println!("--- 2D vs 3D Heterogeneous Inflow ---\n");
    
    println!("2D Heterogeneous Inflow:");
    println!("  - Speed varies in x-direction (horizontal)");
    println!("  - Typically models terrain speedups");
    println!("  - Single slice through the domain\n");
    
    println!("3D Heterogeneous Inflow:");
    println!("  - Speed varies in x and z directions");
    println!("  - Models vertical speed profiles");
    println!("  - More comprehensive but slower\n");
    
    // ============================================================
    // Test 2D vs 3D effects
    // ============================================================
    println!("--- Model Comparison ---\n");
    
    let mut base_model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    let mut hetero_model = florus::FlorisModel::from_file("examples/inputs/gch_heterogeneous_inflow.yaml")?;
    
    // Single turbine
    base_model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    hetero_model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    let wind_speeds = vec![5.0, 8.0, 10.0, 12.0, 15.0];
    
    println!("{:>8} {:>12} {:>12} {:>12}", "WS (m/s)", "Base (kW)", "Hetero (kW)", "Ratio");
    println!("{}", "-".repeat(48));
    
    for ws in &wind_speeds {
        base_model.set_wind_conditions(
            Array1::from_vec(vec![*ws]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        base_model.run()?;
        let base_power = base_model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        hetero_model.set_wind_conditions(
            Array1::from_vec(vec![*ws]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        hetero_model.run()?;
        let hetero_power = hetero_model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        let ratio = hetero_power / base_power;
        
        println!("{:>8.1} {:>12.1} {:>12.1} {:>12.3}", ws, base_power, hetero_power, ratio);
    }
    
    // ============================================================
    // Multi-turbine array
    // ============================================================
    println!("\n--- Multi-Turbine Array (5D spacing) ---\n");
    
    let d = 126.0;
    let spacing = 5.0 * d;
    
    let layout_x = Array1::from_vec(vec![0.0, spacing, 2.0 * spacing]);
    let layout_y = Array1::from_vec(vec![0.0; 3]);
    
    base_model.set_layout(layout_x.clone(), layout_y.clone())?;
    hetero_model.set_layout(layout_x.clone(), layout_y.clone())?;
    
    println!("Layout: 3 turbines at {:.0}D spacing\n", spacing / d);
    
    let ws = 8.0;
    
    // Base model
    base_model.set_wind_conditions(
        Array1::from_vec(vec![ws]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    base_model.run()?;
    let base_powers = base_model.get_turbine_powers();
    
    // Heterogeneous model
    hetero_model.set_wind_conditions(
        Array1::from_vec(vec![ws]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    hetero_model.run()?;
    let hetero_powers = hetero_model.get_turbine_powers();
    
    println!("{:>8} {:>12} {:>12}", "Turbine", "Base (kW)", "Hetero (kW)");
    println!("{}", "-".repeat(35));
    
    for i in 0..3 {
        println!("{:>8} {:>12.1} {:>12.1}", i, base_powers[[0, i]] / 1000.0, hetero_powers[[0, i]] / 1000.0);
    }
    
    let base_total: f64 = (0..3).map(|i| base_powers[[0, i]]).sum::<f64>() / 1000.0;
    let hetero_total: f64 = (0..3).map(|i| hetero_powers[[0, i]]).sum::<f64>() / 1000.0;
    
    println!("{:>8} {:>12.1} {:>12.1}", "Total", base_total, hetero_total);
    
    println!("\n--- Key Observations ---\n");
    
    println!("1. 2D heterogeneous inflow affects all turbines equally in x-direction");
    println!("2. 3D heterogeneous inflow can model vertical speed profiles");
    println!("3. Speedups are typically 5-20% in complex terrain");
    println!("4. Effect varies with wind direction");
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
