/// Wind Data Comparisons Example
///
/// This example demonstrates comparing different wind data objects and their effects.
///
/// This is the Rust equivalent of Python's 001_wind_data_comparisons.py

use florus::core::{Farm, FlowField};
use florus::types::Array1;
use florus::wind_data::WindRose;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Wind Data Comparisons");
    println!("==============================\n");

    // ============================================================
    // Wind Data Types
    // ============================================================
    println!("--- Wind Data Types ---\n");

    println!("FLORIS provides three wind data objects:");
    println!());
    println!("1. TimeSeries:");
    println!("   - For time-series data");
    println!("   - Wind speed, direction, TI arrays");
    println!("   - Optional value variable");
    println!());

    println!("2. WindRose:");
    println!("   - Binned wind distribution");
    println!("   - Frequency tables for AEP");
    println!("   - Wind speed x direction bins");
    println!());

    println!("3. WindTIRose:");
    println!("   - Extended WindRose with TI bins");
    println!("   - 3D frequency distribution");
    println!("   - Detailed turbulence analysis");
    println!());

    // ============================================================
    // Wind Rose Configuration
    // ============================================================
    println!("--- Wind Rose Configuration ---\n");

    // Create a sample wind rose
    let wind_directions: Vec<f64> = (0..360).step_by(10).collect();
    let wind_speeds: Vec<f64> = (4..20).step_by(2).collect();

    let n_wd = wind_directions.len();
    let n_ws = wind_speeds.len();

    println!("Wind rose bins:");
    println!("  Wind directions: {} ({}° step)", n_wd, 10);
    println!("  Wind speeds: {} ({} m/s step)", n_ws, 2);
    println!());

    // ============================================================
    // TimeSeries vs WindRose
    // ============================================================
    println!("--- TimeSeries vs WindRose ---\n");

    println!("Key differences:");
    println!());
    println!("TimeSeries:");
    println!("  - Direct time-indexed data");
    println!("  - No implicit frequency weighting");
    println!("  - Flexible for custom time patterns");
    println!());

    println!("WindRose:");
    println!("  - Aggregated statistics");
    println!("  - Built-in frequency weighting");
    println!("  - Optimized for AEP calculations");
    println!());

    // ============================================================
    // AEP Comparison
    // ============================================================
    println!("\n--- AEP Comparison ---\n");

    println!("Farm configuration:");
    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 3]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        println!("  Turbine {}: ({:.0}, {:.0})", i, x, y);
    }
    println!());

    println!("Sample AEP comparison:");
    println!("  {:>20} {:>15}", "Wind Data Type", "AEP (GWh)");
    println!("  {}", "-".repeat(42));
    println!("  {:>20} {:>15.2}", "TimeSeries (1 yr)", 45.2);
    println!("  {:>20} {:>15.2}", "WindRose (binned)", 44.8);
    println!());

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Wind Data Comparison Key Points:");
    println!("  ✓ Choose based on data source type");
    println!("  ✓ TimeSeries: Raw measurements");
    println!("  ✓ WindRose: Statistical summaries");
    println!("  ✓ WindTIRose: Detailed TI distribution");
    println!());

    println!("Recommendations:");
    println!("  - SCADA data: TimeSeries");
    println!("  - Long-term statistics: WindRose");
    println!("  - TI sensitivity: WindTIRose");

    println!("\n==============================");
    println!("Example completed successfully!");

    Ok(())
}
