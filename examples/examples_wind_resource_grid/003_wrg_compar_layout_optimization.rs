/// Example: Layout Optimization with WindRoseWRG Comparison
///
/// This example compares layout optimization using different wind data sources.
///
/// This is the Rust equivalent of Python's examples_wind_resource_grid/003_wrg_compar_layout_optimization.py
///
/// Note: This is a conceptual example for future WRG support.

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Layout Optimization with WindRoseWRG Comparison");
    println!("==============================================================\n");

    println!("--- Overview ---\n");
    
    println!("Layout optimization methods:");
    println!("  1. WindRose: Sector-averaged (standard)");
    println!("  2. WindRoseWRG: Spatially-varying (future)");
    println!("  3. TimeSeries: Time-varying (detailed)\n");
    
    println!("--- Current Implementation ---\n");
    
    println!("Current FLORIS-RS supports WindRose-based optimization.\n");
    
    // ============================================================
    // Show WindRose optimization example
    // ============================================================
    println!("--- WindRose Optimization ---\n");
    
    let wind_rose = florus::wind_data::WindRose::from_csv_long(
        "examples/inputs/wind_rose.csv",
        "wd", "ws", "freq_val", 0.06
    )?;
    
    println!("WindRose: {} directions, {} speeds", 
        wind_rose.n_dir(), wind_rose.n_ws());
    
    // AEP calculation would go here
    
    println!("\n--- Future: WindRoseWRG ---\n");
    
    println!("WindRoseWRG would allow:");
    println!("  - Optimization at each grid point");
    println!("  - Spatially-varying optimal layout");
    println!("  - More accurate site-specific AEP\n");
    
    println!("====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
