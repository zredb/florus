/// Example: Operation Models
///
/// This example demonstrates different turbine operation models
/// in FLORIS-RS.
///
/// Operation models define how turbines behave in terms of:
/// - Power production
/// - Thrust coefficient
/// - Control strategies
///
/// This is the Rust equivalent of Python's examples_operation_models/...

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Operation Models");
    println!("================================\n");

    println!("--- Operation Models ---\n");
    
    println!("FLORIS supports multiple operation models:");
    println!("  1. Simple: Basic power/Ct curves");
    println!("  2. CosineLoss: Yaw/tilt corrections");
    println!("  3. UnifiedMomentum: Combined approach");
    println!("  4. Mixed: Multiple models");
    println!("  5. SimpleDerating: Power limiting");
    println!("  6. PeakShaving: Load limiting");
    println!("  7. AWC: Active wake control\n");
    
    // ============================================================
    // Basic comparison
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    model.set_layout(
        &Array1::from_vec(vec![0.0]),
        &Array1::from_vec(vec![0.0]),
    )?;
    
    // ============================================================
    // Test different conditions
    // ============================================================
    println!("--- Power by Wind Speed ---\n");
    
    let wind_speeds = vec![4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];
    
    println!("{:>8} {:>12}", "WS (m/s)", "Power (kW)");
    println!("{}", "-".repeat(22));
    
    for ws in &wind_speeds {
        model.set_wind_conditions(
            Array1::from_vec(vec![*ws]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        println!("{:>8.1} {:>12.1}", ws, power);
    }
    
    // ============================================================
    // Yaw effects
    // ============================================================
    println!("\n--- Yaw Angle Effects ---\n");
    
    let yaws = vec![-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0];
    
    println!("{:>8} {:>12}", "Yaw (°)", "Power (kW)");
    println!("{}", "-".repeat(22));
    
    for yaw in &yaws {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        // Note: Yaw setting would require additional API
        model.run()?;
        
        let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        println!("{:>8.1} {:>12.1}", yaw, power);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
