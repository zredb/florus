/// Example: Tilt-Driven Vertical Wake Deflection
///
/// This example demonstrates vertical wake deflection due to floating turbine tilt.
/// When a floating turbine tilts, its wake is deflected vertically, which can
/// affect downstream turbine placement and power production.
///
/// This is the Rust equivalent of Python's examples_floating/003_tilt_driven_vertical_wake_deflection.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Tilt-Driven Vertical Wake Deflection");
    println!("========================================================\n");

    println!("--- Vertical Wake Deflection ---\n");
    
    println!("Floating turbine tilt causes vertical wake deflection:");
    println!("  - Upward tilt pushes wake upward");
    println!("  - Downward tilt (rare) pushes wake downward");
    println!("  - Effect increases with tilt angle\n");
    
    println!("This effect is important for:");
    println!("  - Array layout optimization");
    println!("  - Vertical wake interaction");
    println!("  - Extreme tilt conditions\n");
    
    // ============================================================
    // Compare fixed vs floating at multiple wind speeds
    // ============================================================
    println!("--- Power Comparison: Fixed vs Floating ---\n");
    
    // Use Empirical Gauss model for tilt effects
    let configs = vec![
        ("examples/inputs_floating/emgauss_fixed.yaml", "Fixed Bottom"),
        ("examples/inputs_floating/emgauss_floating.yaml", "Floating"),
        ("examples/inputs_floating/emgauss_floating_fixedtilt5.yaml", "Fixed Tilt 5°"),
        ("examples/inputs_floating/emgauss_floating_fixedtilt15.yaml", "Fixed Tilt 15°"),
    ];
    
    let wind_speeds = vec![6.0, 8.0, 10.0, 12.0, 15.0, 18.0];
    let wind_directions = Array1::from_vec(vec![270.0; wind_speeds.len()]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; wind_speeds.len()]);
    
    // Single turbine
    println!("Single turbine power:");
    println!("{:>8}", "WS (m/s)");
    for (_, name) in &configs {
        print!(" {:>14}", name);
    }
    println!();
    println!("{}", "-".repeat14 * configs.len(8 + () + 1));
    
    for ws in &wind_speeds {
        print!("{:>8.1}", ws);
        
        for (config, _) in &configs {
            let ws_array = Array1::from_vec(vec![*ws; wind_speeds.len()]);
            
            if let Ok(mut model) = florus::FlorisModel::from_file(config) {
                model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
                model.set_wind_conditions(
                    ws_array,
                    wind_directions.clone(),
                    turbulence_intensities.clone(),
                )?;
                if model.run().is_ok() {
                    let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
                    print!(" {:>14.1}", power);
                } else {
                    print!(" {:>14}", "Error");
                }
            } else {
                print!(" {:>14}", "Error");
            }
        }
        println!();
    }
    
    // ============================================================
    // 2-turbine array
    // ============================================================
    println!("\n--- Two-Turbine Array (5D spacing) ---\n");
    
    let d = 126.0;
    let spacing = 5.0 * d;
    let layout_x = Array1::from_vec(vec![0.0, spacing]);
    let layout_y = Array1::from_vec(vec![0.0; 2]);
    
    println!("Layout: 2 turbines at {:.0}D spacing", spacing / d);
    println!("{:>8} {:>14} {:>14} {:>14} {:>14}", "WS", "Upstream", "Fixed", "Float", "15° Tilt");
    println!("{}", "-".repeat(65));
    
    for ws in &wind_speeds {
        let ws_array = Array1::from_vec(vec![*ws; wind_speeds.len()]);
        
        // Fixed
        if let Ok(mut model) = florus::FlorisModel::from_file("examples/inputs_floating/emgauss_fixed.yaml") {
            model.set_layout(layout_x.clone(), layout_y.clone())?;
            model.set_wind_conditions(ws_array.clone(), wind_directions.clone(), turbulence_intensities.clone())?;
            if model.run().is_ok() {
                let powers = model.get_turbine_powers();
                let upstream_fixed = powers[[0, 0]] / 1000.0;
                let downstream_fixed = powers[[0, 1]] / 1000.0;
                print!("{:>8.1} {:>14.1} {:>14.1}", ws, upstream_fixed, downstream_fixed);
            } else {
                print!("{:>8.1} {:>14} {:>14}", ws, "Error", "Error");
            }
        }
        
        // Floating
        if let Ok(mut model) = florus::FlorisModel::from_file("examples/inputs_floating/emgauss_floating.yaml") {
            model.set_layout(layout_x.clone(), layout_y.clone())?;
            model.set_wind_conditions(ws_array.clone(), wind_directions.clone(), turbulence_intensities.clone())?;
            if model.run().is_ok() {
                let power = model.get_turbine_powers()[[0, 1]] / 1000.0;
                print!(" {:>14.1}", power);
            } else {
                print!(" {:>14}", "Error");
            }
        }
        
        // Fixed Tilt 15
        if let Ok(mut model) = florus::FlorisModel::from_file("examples/inputs_floating/emgauss_floating_fixedtilt15.yaml") {
            model.set_layout(layout_x.clone(), layout_y.clone())?;
            model.set_wind_conditions(ws_array, wind_directions.clone(), turbulence_intensities.clone())?;
            if model.run().is_ok() {
                let power = model.get_turbine_powers()[[0, 1]] / 1000.0;
                print!(" {:>14.1}", power);
            } else {
                print!(" {:>14}", "Error");
            }
        }
        
        println!();
    }
    
    println!("\n--- Key Observations ---\n");
    println!("1. Floating turbines have different power curves due to tilt");
    println!("2. Downstream power varies with tilt angle");
    println!("3. Vertical wake deflection affects array performance");
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
