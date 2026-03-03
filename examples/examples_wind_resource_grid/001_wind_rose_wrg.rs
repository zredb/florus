/// Example: WindRoseWRG
///
/// WindRoseWRG is a type of WindData object that combines WindRose
/// functionality with Wind Resource Grid data.
///
/// This is the Rust equivalent of Python's examples_wind_resource_grid/001_wind_rose_wrg.py
///
/// Note: This is a placeholder showing the WindRoseWRG concept.

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: WindRoseWRG");
    println!("============================\n");

    println!("--- WindRoseWRG Overview ---\n");
    
    println!("WindRoseWRG combines:");
    println!("  - WindRose: Direction-sector averaged data");
    println!("  - WRG: Spatial wind resource variation");
    println!("  - Used for spatially-varying AEP calculations\n");
    
    println!("--- Current Implementation ---\n");
    
    println!("Note: Full WindRoseWRG support in FLORIS-RS is future work.");
    println!("Currently, use WindRose for sector-averaged calculations.\n");
    
    // ============================================================
    // Demonstrate current WindRose functionality
    // ============================================================
    println!("--- Using WindRose ---\n");
    
    // Load wind rose from CSV
    let wind_rose = florus::wind_data::WindRose::from_csv_long(
        "examples/inputs/wind_rose.csv",
        "wd", "ws", "freq_val", 0.06
    )?;
    
    println!("Wind Rose loaded:");
    println!("  Number of directions: {}", wind_rose.n_dir());
    println!("  Number of wind speeds: {}", wind_rose.n_ws());
    
    println!("\n--- WRG Placeholder ---\n");
    
    println!("WindRoseWRG would allow:");
    println!("  1. Load WRG file");
    println!("  2. Extract wind data at each grid point");
    println!("  3. Calculate AEP at each location");
    println!("  4. Spatial optimization\n");
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
