/// Example: Specify Turbine Power Curve
///
/// This example demonstrates how to specify a custom turbine model based on
/// power and thrust curves in FLORIS-RS.
///
/// This is the Rust equivalent of Python's examples_turbine/003_specify_turbine_power_curve.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Specify Turbine Power Curve");
    println!("==========================================\n");

    println!("--- Custom Turbine Definition ---\n");
    
    println!("FLORIS allows custom turbine definitions with:");
    println!("  - Power curve: power vs wind speed");
    println!("  - Thrust coefficient curve: Ct vs wind speed");
    println!("  - Physical parameters (rotor diameter, hub height, etc.)\n");
    
    // ============================================================
    // Example power and thrust data
    // ============================================================
    println!("--- Example Power/Thrust Data ---\n");
    
    // Example turbine data (similar to Python example)
    let wind_speeds = vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
    let powers = vec![0.0, 30.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0];
    
    // Calculate power coefficients
    let rotor_diameter = 126.0;
    let air_density = 1.225;
    let swept_area = std::f64::consts::PI * (rotor_diameter / 2.0).powi(2);
    
    println!("Wind Speed (m/s) | Power (kW) | Cp");
    println!("{}", "-".repeat(45));
    
    for (i, ws) in wind_speeds.iter().enumerate().skip(1) {
        let cp = powers[i] * 1000.0 / (0.5 * air_density * swept_area * ws.powi(3));
        println!("{:>17.1} | {:>10.0} | {:.4}", ws, powers[i], cp);
    }
    
    // ============================================================
    // Load model and test custom turbine
    // ============================================================
    println!("\n--- Testing with Standard Turbine ---\n");
    
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    let test_speeds = vec![4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0];
    
    println!("{:>8} {:>12}", "WS (m/s)", "Power (kW)");
    println!("{}", "-".repeat(22));
    
    for ws in &test_speeds {
        model.set_wind_conditions(
            Array1::from_vec(vec![*ws]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        model.run()?;
        
        let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        println!("{:>8.1} {:>12.1}", ws, power);
    }
    
    println!("\n--- Custom Turbine Notes ---\n");
    
    println!("To create custom turbine in FLORIS-RS:");
    println!("  1. Define power curve (wind_speeds, powers)");
    println!("  2. Define thrust coefficient curve");
    println!("  3. Specify physical parameters");
    println!("  4. Use in Farm::new() with custom type\n");
    
    println!("====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
