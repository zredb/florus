/// Example: Wind Resource Grid
///
/// This example demonstrates Wind Resource Grid (WRG) usage in FLORIS-RS.
/// WRG files contain spatially-varying wind data across a wind farm.
///
/// This is the Rust equivalent of Python's examples_wind_resource_grid/002_wind_resource_grid.py
///
/// Note: This is a conceptual example.

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Wind Resource Grid");
    println!("====================================\n");

    println!("--- Wind Resource Grid ---\n");
    
    println!("Wind Resource Grid (WRG) provides:");
    println!("  - Spatially-varying wind conditions");
    println!("  - Multiple measurement locations");
    println!("  - Used for site characterization\n");
    
    println!("--- Current FLORIS-RS Support ---\n");
    
    println!("Currently, FLORIS-RS supports:");
    println!("  - TimeSeries: Time-varying conditions");
    println!("  - WindRose: Sector-averaged conditions");
    println!("  - WindTIRose: TI-averaged conditions\n");
    
    println!("WRG functionality is planned for future releases.\n");
    
    // ============================================================
    // Demonstrate WindRose
    // ============================================================
    println!("--- WindRose Demo ---\n");
    
    // Load wind rose
    let wind_rose = florus::wind_data::WindRose::from_csv_long(
        "examples/inputs/wind_rose.csv",
        "wd", "ws", "freq_val", 0.06
    )?;
    
    println!("WindRose statistics:");
    println!("  Directions: {} sectors", wind_rose.n_dir());
    println!("  Wind speeds: {} bins", wind_rose.n_ws());
    
    // Get wind rose data
    let ws = wind_rose.wind_speeds();
    let freq = wind_rose.frequency_table()?;
    
    println!("\n{:>8} {:>12}", "WS (m/s)", "Frequency %");
    println!("{}", "-".repeat(22));
    
    for (i, w) in ws.iter().enumerate() {
        if i < 10 { // First 10 speeds
            println!("{:>8.1} {:>12.2}", w, freq[[i, 0]] * 100.0);
        }
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
