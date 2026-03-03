/// Example: Floating Offshore Wind Turbine Overview
///
/// This example provides an overview of floating offshore wind turbines in FLORIS-RS.
/// Floating turbines differ from fixed-bottom turbines in that they have:
/// - Dynamic tilt angles based on wind speed
/// - Tilt-dependent power and thrust coefficients
/// - Additional wake deflection due to tilt
///
/// This is the Rust equivalent of Python's examples_floating/001_floating_overview.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Floating Offshore Wind Turbine Overview");
    println!("============================================================\n");

    println!("--- Floating Turbine Concepts ---\n");
    
    println!("1. FLOATING TURBINE BASICS:");
    println!("   - Floating turbines are mounted on floating platforms");
    println!("   - They experience dynamic tilt (pitch) angles");
    println!("   - Platform motion affects power and wake behavior\n");
    
    println!("2. TILT ANGLE EFFECTS:");
    println!("   - Higher tilt at low wind speeds (peak power tracking)");
    println!("   - Lower tilt at high wind speeds (structural stability)");
    println!("   - Tilt affects both power capture and wake deflection\n");
    
    println!("3. WAKE IMPLICATIONS:");
    println!("   - Tilt-induced vertical wake deflection");
    println!("   - Platform stability affects wake coherence");
    println!("   - Different from fixed-bottom turbine wakes\n");
    
    // ============================================================
    // Load floating turbine model
    // ============================================================
    println!("--- Loading Floating Turbine Model ---\n");
    
    let mut model = florus::FlorisModel::from_file("examples/inputs_floating/gch_floating.yaml")?;
    
    println!("Model loaded successfully!");
    println!("  Configuration: Gauss-Curl-Hybrid with floating turbine");
    
    // Set up a single turbine
    model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    // Test different wind speeds
    let wind_speeds = vec![5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0];
    let wind_directions = Array1::from_vec(vec![270.0; wind_speeds.len()]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; wind_speeds.len()]);
    
    println!("\n--- Wind Speed Sweep ---\n");
    println!("{:>8} {:>12} {:>12}", "WS (m/s)", "Power (kW)", "Ct");
    println!("{}", "-".repeat(35));
    
    for ws in &wind_speeds {
        let ws_array = Array1::from_vec(vec![*ws; wind_speeds.len()]);
        model.set_wind_conditions(
            ws_array.clone(),
            wind_directions.clone(),
            turbulence_intensities.clone(),
        )?;
        model.run()?;
        
        let powers = model.get_turbine_powers();
        let cts = model.get_turbine_thrust_coefficients();
        
        println!("{:>8.1} {:>12.1} {:>12.4}", ws, powers[[0, 0]] / 1000.0, cts[[0, 0]]);
    }
    
    println!("\n--- Floating vs Fixed Comparison ---\n");
    
    // Compare with fixed turbine
    let mut fixed_model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fixed_model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    println!("{:>8} {:>12} {:>12} {:>12}", "WS (m/s)", "Float (kW)", "Fixed (kW)", "Diff %");
    println!("{}", "-".repeat(50));
    
    for ws in &wind_speeds {
        let ws_array = Array1::from_vec(vec![*ws; wind_speeds.len()]);
        
        model.set_wind_conditions(
            ws_array.clone(),
            wind_directions.clone(),
            turbulence_intensities.clone(),
        )?;
        model.run()?;
        let float_power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        fixed_model.set_wind_conditions(
            ws_array.clone(),
            wind_directions.clone(),
            turbulence_intensities.clone(),
        )?;
        fixed_model.run()?;
        let fixed_power = fixed_model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        let diff = (float_power - fixed_power) / fixed_power * 100.0;
        println!("{:>8.1} {:>12.1} {:>12.1} {:>12.1}%", ws, float_power, fixed_power, diff);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
