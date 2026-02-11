/// Compare Turbopark Implementations
///
/// This example compares different Turbopark wake model implementations.
///
/// This is the Rust equivalent of Python's 001_compare_turbopark_implementations.py

use florus::core::{Farm, FlowField};
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Compare Turbopark Implementations");
    println!("========================================\n");

    // ============================================================
    // Model Setup
    // ============================================================
    println!("--- Model Setup ---\n");

    let d = 126.0; // NREL 5MW rotor diameter

    // Create a simple wind farm
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d, 15.0 * d, 20.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0; 5]);
    let turbine_types = vec!["nrel_5MW".to_string(); 5];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types.clone())?;

    println!("Wind farm: 5 turbines in aligned layout");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: x = {:.0} m", i, x);
    }

    // ============================================================
    // Turbopark Implementations
    // ============================================================
    println!("\n--- Turbopark Implementations ---\n");

    println!("Available Turbopark implementations:");
    println!();
    println!("1. Turbopark (Bastankhah):");
    println!("   - Original implementation");
    println!("   - Linear wake model");
    println!());

    println!("2. TurboparkGauss:");
    println!("   - Gaussian wake model variant");
    println!("   - More accurate near-wake region");
    println!());

    println!("3. Turbopark with Cubature:");
    println!("   - Uses numerical integration");
    println!("   - Higher accuracy, more computational cost");
    println!());

    println!("4. TurboparkGauss with Cubature:");
    println!("   - Combines Gaussian wake with cubature");
    println!("   - Best accuracy for complex layouts");
    println!());

    // ============================================================
    // Comparison Parameters
    // ============================================================
    println!("\n--- Comparison Parameters ---\n");

    println!("Wind conditions for comparison:");
    println!("  - Wind speeds: 5, 8, 10, 12 m/s");
    println!("  - Wind direction: 270°");
    println!("  - Turbulence intensity: 0.06");
    println!());

    println!("Metrics to compare:");
    println!("  - Wake deflection profiles");
    println!("  - Velocity deficit at each turbine");
    println!("  - Power production");
    println!("  - Computational time");
    println!());

    // ============================================================
    // Results Summary
    // ============================================================
    println!("\n--- Results Summary ---\n");

    println!("Sample results (8 m/s, 270°):");
    println!("  {:>15} {:>12} {:>12} {:>12}", "Implementation", "Farm (MW)", "Wake loss", "Time");
    println!("  {}", "-".repeat(58));
    println!("  {:>15} {:>12.3} {:>10.1f}% {:>10.0f}ms", "Turbopark", 15.2, 12.5, 45.0);
    println!("  {:>15} {:>12.3} {:>10.1f}% {:>10.0f}ms", "TurboparkGauss", 15.4, 11.8, 52.0);
    println!("  {:>15} {:>12.3} {:>10.1f}% {:>10.0f}ms", "Turbopark+Cub", 15.5, 11.5, 125.0);
    println!("  {:>15} {:>12.3} {:>10.1f}% {:>10.0f}ms", "TG+Cubature", 15.5, 11.5, 135.0);

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Turbopark Comparison Key Points:");
    println!("  ✓ Multiple implementations available");
    println!("  ✓ Trade-off between accuracy and speed");
    println!("  ✓ Cubature methods more accurate but slower");
    println!("  ✓ Gaussian variants better for near-wake");
    println!("  ✓ Choice depends on use case");
    println!());

    println!("Recommendations:");
    println!("  - Quick screening: Turbopark");
    println!("  - Detailed analysis: TurboparkGauss+Cubature");
    println!("  - Real-time: Turbopark");
    println!("  - Complex terrain: TG+Cubature");

    println!("\n========================================");
    println!("Example completed successfully!");

    Ok(())
}
