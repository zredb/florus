//! Example: Heterogeneous with Wind Rose
//!
//! This example demonstrates heterogeneous inflow with wind rose analysis.
//!
//! Corresponds to: examples_heterogeneous/004_het_with_wind_rose.rs

use florus::Result;

fn main() -> Result<()> {
    println!("=== Heterogeneous with Wind Rose ===\n");

    println!("Heterogeneous Inflow with Wind Rose:\\n");
    
    println!("Concept:");
    println!("  - Combine spatial wind variation with directional distribution");
    println!("  - Most realistic wind farm modeling approach");
    println!("  - Captures full complexity of wind resources\\n");

    println!("Components:");
    println!("  1. Wind Rose:");
    println!("     - Wind direction frequencies");
    println!("     - Wind speed distributions");
    println!("     - TI variations\\n");

    println!("  2. Heterogeneous Map:");
    println!("     - Spatial wind speed variation");
    println!("     - Direction-dependent patterns");
    println!("     - Terrain effects\\n");

    println!("  3. Integration:");
    println!("     - For each wind direction:");
    println!("       * Apply corresponding heterogeneous map");
    println!("       * Calculate farm power");
    println!("       * Weight by frequency");
    println!("     - Sum for total AEP\\n");

    println!("Benefits:");
    println!("  - Accurate AEP estimation");
    println!("  - Realistic wake modeling");
    println!("  - Better layout optimization");
    println!("  - Improved financial projections\\n");

    println!("Challenges:");
    println!("  - Computational cost (many scenarios)");
    println!("  - Data requirements (detailed maps)");
    println!("  - Model complexity");
    println!("  - Validation difficulty\\n");

    println!("Applications:");
    println!("  - Final design validation");
    println!("  - Financial modeling");
    println!("  - Permitting studies");
    println!("  - Research applications\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires WindRose and HeterogeneousMap integration.");
    Ok(())
}
