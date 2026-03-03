/// Example: Set Wind Conditions
///
/// This example demonstrates how to set wind conditions in FLORIS-RS,
/// including wind speed, wind direction, turbulence intensity, and more.
///
/// This is the Rust equivalent of Python's 004_set.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 004: Setting Wind Conditions");
    println!("==============================================\n");

    println!("--- Wind Condition Parameters ---\n");
    
    println!("Key wind conditions to set:");
    println!("  - Wind speed (m/s)");
    println!("  - Wind direction (degrees)");
    println!("  - Turbulence intensity");
    println!("  - Air density (optional)");
    println!("  - Wind shear exponent\n");
    
    // ============================================================
    // Basic setup
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    model.set_layout(
        &Array1::from_vec(vec![0.0]),
        &Array1::from_vec(vec![0.0]),
    )?;
    
    // ============================================================
    // Set wind speed
    // ============================================================
    println!("--- Setting Wind Speed ---\n");
    
    let wind_speeds = vec![4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
    
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
    // Set wind direction
    // ============================================================
    println!("\n--- Setting Wind Direction ---\n");
    
    let directions: Vec<f64> = (180..360).step_by(20).map(|d| d as f64).collect();
    
    println!("{:>8} {:>12}", "WD (°)", "Power (kW)");
    println!("{}", "-".repeat(22));
    
    for wd in &directions {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![*wd]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        println!("{:>8.0} {:>12.1}", wd, power);
    }
    
    // ============================================================
    // Set turbulence intensity
    // ============================================================
    println!("\n--- Setting Turbulence Intensity ---\n");
    
    let tis = vec![0.03, 0.05, 0.07, 0.10, 0.15, 0.20];
    
    println!("{:>8} {:>12}", "TI", "Power (kW)");
    println!("{}", "-".repeat(22));
    
    for ti in &tis {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![*ti]),
        )?;
        
        model.run()?;
        
        let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
        println!("{:>8.2} {:>12.1}", ti, power);
    }
    
    // ============================================================
    // Multiple conditions at once
    // ============================================================
    println!("\n--- Multiple Wind Conditions ---\n");
    
    // Multiple wind speeds and directions
    let ws_array = Array1::from_vec(vec![6.0, 8.0, 10.0, 12.0]);
    let wd_array = Array1::from_vec(vec![270.0, 280.0, 290.0, 300.0]);
    let ti_array = Array1::from_vec(vec![0.06, 0.07, 0.08, 0.09]);
    
    model.set_wind_conditions(ws_array, wd_array.clone(), ti_array)?;
    
    model.run()?;
    
    let powers = model.get_turbine_powers();
    
    println!("{:>8} {:>8} {:>8} {:>12}", "WS", "WD", "TI", "Power");
    println!("{}", "-".repeat(40));
    
    for i in 0..4 {
        println!("{:>8.1} {:>8.0} {:>8.2} {:>12.1}", 
            ws_array[i], wd_array[i], ti_array[i], powers[[i, 0]] / 1000.0);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
