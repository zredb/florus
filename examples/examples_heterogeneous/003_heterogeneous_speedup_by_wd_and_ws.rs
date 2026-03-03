/// Example: Heterogeneous Speedup by Wind Direction and Wind Speed
///
/// This example demonstrates how heterogeneous speedups can vary with
/// both wind direction and wind speed in FLORIS-RS.
///
/// This is the Rust equivalent of Python's examples_heterogeneous/003_heterogeneous_speedup_by_wd_and_ws.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Heterogeneous Speedup by Wind Direction and Speed");
    println!("==================================================================\n");

    println!("--- Heterogeneous Speedups ---\n");
    
    println!("Heterogeneous speedups can vary with:");
    println!("  - Wind direction (different terrain effects)");
    println!("  - Wind speed (atmospheric stability effects)");
    println!("  - Spatial location\n");
    
    // ============================================================
    // Load model with heterogeneous config
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch_heterogeneous_inflow.yaml")?;
    
    // Single turbine
    model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    println!("--- Speedup by Wind Speed ---\n");
    
    // Test different wind speeds
    let wind_speeds = vec![5.0, 8.0, 10.0, 12.0, 15.0, 20.0];
    
    // Without heterogeneous
    let mut base_model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    base_model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    println!("{:>8} {:>12} {:>12} {:>12}", "WS (m/s)", "Homog (kW)", "Hetero (kW)", "Speedup");
    println!("{}", "-".repeat(48));
    
    for ws in &wind_speeds {
        // Base model (homogeneous)
        base_model.set_wind_conditions(
            Array1::from_vec(vec![*ws]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        base_model.run()?;
        let base_power = base_model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        // Heterogeneous model
        model.set_wind_conditions(
            Array1::from_vec(vec![*ws]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        model.run()?;
        let hetero_power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        let speedup = hetero_power / base_power;
        
        println!("{:>8.1} {:>12.1} {:>12.1} {:>12.3}", ws, base_power, hetero_power, speedup);
    }
    
    println!("\n--- Speedup by Wind Direction ---\n");
    
    let directions: Vec<f64> = (240..300).step_by(10).map(|d| d as f64).collect();
    let ws = 8.0;
    
    println!("{:>8} {:>12} {:>12} {:>12}", "WD (°)", "Homog (kW)", "Hetero (kW)", "Speedup");
    println!("{}", "-".repeat(48));
    
    for wd in &directions {
        // Base model (homogeneous)
        base_model.set_wind_conditions(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![*wd]),
            Array1::from_vec(vec![0.06]),
        )?;
        base_model.run()?;
        let base_power = base_model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        // Heterogeneous model
        model.set_wind_conditions(
            Array1::from_vec(vec![ws]),
            Array1::from_vec(vec![*wd]),
            Array1::from_vec(vec![0.06]),
        )?;
        model.run()?;
        let hetero_power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        let speedup = hetero_power / base_power;
        
        println!("{:>8.0} {:>12.1} {:>12.1} {:>12.3}", wd, base_power, hetero_power, speedup);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
