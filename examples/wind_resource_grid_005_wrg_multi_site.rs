//! Example: WRG Multi-Site Analysis
//!
//! This example demonstrates analyzing multiple sites with WRG data.
//!
//! Corresponds to: examples_wind_resource_grid/005_wrg_multi_site.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== WRG Multi-Site Analysis ===\n");

    println!("Multi-Site Wind Resource Analysis:\\n");
    
    println!("Scenario:");
    println!("  - Multiple potential wind farm sites");
    println!("  - Compare wind resources across sites");
    println!("  - Select optimal development location\\n");

    println!("Analysis Process:");
    println!("  1. Load WRG for each site");
    println!("  2. Extract key metrics:");
    println!("     - Mean wind speed");
    println!("     - Wind power density");
    println!("     - Capacity factor estimate");
    println!("     - Turbulence characteristics\\n");

    println!("  3. Compare sites:");
    println!("     - Rank by energy potential");
    println!("     - Assess variability");
    println!("     - Evaluate risks");
    println!("     - Consider constraints\\n");

    println!("Comparison Metrics:");
    println!("  - Annual energy production (AEP)");
    println!("  - Capacity factor (%)");
    println!("  - Net present value (NPV)");
    println!("  - Levelized cost of energy (LCOE)");
    println!("  - Risk-adjusted returns\\n");

    println!("Decision Factors:");
    println!("  - Wind resource quality");
    println!("  - Grid connection availability");
    println!("  - Land acquisition costs");
    println!("  - Environmental constraints");
    println!("  - Permitting complexity\\n");

    println!("Applications:");
    println!("  - Portfolio optimization");
    println!("  - Site screening");
    println!("  - Investment decisions");
    println!("  - Strategic planning\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires WRG parser and multi-site analysis tools.");
    Ok(())
}
