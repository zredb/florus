//! Example: Floating vs Fixed-Bottom Farm Comparison
//!
//! This example compares floating and fixed-bottom wind farm performance.
//!
//! Corresponds to: examples_floating/002_floating_vs_fixedbottom_farm.rs

use florus::Result;

fn main() -> Result<()> {
    println!("=== Floating vs Fixed-Bottom Farm ===\n");

    println!("Comparison of Floating and Fixed-Bottom Farms:\\n");
    
    println!("Fixed-Bottom Turbines:");
    println!("  Advantages:");
    println!("    - Proven technology");
    println!("    - Lower cost in shallow water");
    println!("    - Simpler installation");
    println!("    - Stable platform\\n");
    
    println!("  Limitations:");
    println!("    - Depth limit (~60m)");
    println!("    - Foundation costs increase with depth");
    println!("    - Limited site availability\\n");

    println!("Floating Turbines:");
    println!("  Advantages:");
    println!("    - Access to deep water sites");
    println!("    - Higher wind resources offshore");
    println!("    - Reduced visual impact");
    println!("    - Potential for larger farms\\n");
    
    println!("  Challenges:");
    println!("    - Higher initial cost");
    println!("    - Complex dynamics");
    println!("    - Mooring system design");
    println!("    - Maintenance access\\n");

    println!("Performance Differences:");
    println!("  - Floating: Platform motion affects power");
    println!("  - Floating: Tilt control for wake steering");
    println!("  - Floating: Wave-dependent CP/CT tables");
    println!("  - Fixed: More predictable performance\\n");

    println!("Economic Considerations:");
    println!("  - Break-even depth: ~60-80m");
    println!("  - Floating LCOE decreasing rapidly");
    println!("  - Site-specific optimization needed\\n");

    println!("=== Example Complete ===");
    println!("Note: Full comparison requires both turbine types configured.");
    Ok(())
}
