/// Example: Fixed Tilt Floating Turbine
///
/// This example demonstrates fixed tilt angles for floating turbines.
/// Fixed tilt is commonly used when platform motion is restricted or
/// for comparison studies.
///
/// This is the Rust equivalent of Python's examples_floating/002_fixed_tilt.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Fixed Tilt Floating Turbine");
    println!("==============================================\n");

    println!("--- Fixed Tilt Concept ---\n");
    
    println!("Fixed tilt configurations:");
    println!("  - FixedTilt5: 5 degree tilt (similar to onshore)");
    println!("  - FixedTilt15: 15 degree tilt (extreme platform tilt)");
    println!("  - DefinedFloating: Custom tilt schedule\n");
    
    // Test different tilt configurations
    let configs = vec![
        ("gch_floating.yaml", "Dynamic Tilt"),
        ("gch_floating_fixedtilt5.yaml", "Fixed 5°"),
        ("gch_floating_fixedtilt15.yaml", "Fixed 15°"),
    ];
    
    let wind_speeds = vec![5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0];
    let wind_directions = Array1::from_vec(vec![270.0; wind_speeds.len()]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; wind_speeds.len()]);
    
    println!("--- Power Output by Configuration ---\n");
    println!("{:>8}", "WS (m/s)");
    for (_, name) in &configs {
        print!(" {:>12}", name);
    }
    println!();
    println!("{}", "-".repeat(8 + 12 * configs.len() + 1));
    
    for ws in &wind_speeds {
        print!("{:>8.1}", ws);
        
        for (config, _) in &configs {
            let config_path = format!("examples/inputs_floating/{}", config);
            let ws_array = Array1::from_vec(vec![*ws; wind_speeds.len()]);
            
            if let Ok(mut model) = florus::FlorisModel::from_file(&config_path) {
                model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
                model.set_wind_conditions(
                    ws_array,
                    wind_directions.clone(),
                    turbulence_intensities.clone(),
                )?;
                if let Ok(()) = model.run() {
                    let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
                    print!(" {:>12.1}", power);
                } else {
                    print!(" {:>12}", "N/A");
                }
            } else {
                print!(" {:>12}", "Error");
            }
        }
        println!();
    }
    
    // ============================================================
    // Compare wake effects
    // ============================================================
    println!("\n--- Wake Effects with Tilt ---\n");
    
    // Create a 2-turbine layout
    let d = 126.0;
    let spacing = 5.0 * d;
    let layout_x = Array1::from_vec(vec![0.0, spacing]);
    let layout_y = Array1::from_vec(vec![0.0; 2]);
    
    println!("2-turbine layout with {:.0}D spacing:", spacing / d);
    println!("{:>8} {:>12} {:>12} {:>12}", "WS (m/s)", "Dyn. Tilt", "Fix 5°", "Fix 15°");
    println!("{}", "-".repeat(50));
    
    for ws in &wind_speeds {
        print!("{:>8.1}", ws);
        
        for (config, _) in &configs {
            let config_path = format!("examples/inputs_floating/{}", config);
            let ws_array = Array1::from_vec(vec![*ws; wind_speeds.len()]);
            
            if let Ok(mut model) = florus::FlorisModel::from_file(&config_path) {
                model.set_layout(layout_x.clone(), layout_y.clone())?;
                model.set_wind_conditions(
                    ws_array,
                    wind_directions.clone(),
                    turbulence_intensities.clone(),
                )?;
                if let Ok(()) = model.run() {
                    let powers = model.get_turbine_powers();
                    // Second turbine power (affected by wake)
                    let downstream_power = powers[[0, 1]] / 1000.0;
                    print!(" {:>12.1}", downstream_power);
                } else {
                    print!(" {:>12}", "N/A");
                }
            } else {
                print!(" {:>12}", "Error");
            }
        }
        println!();
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
