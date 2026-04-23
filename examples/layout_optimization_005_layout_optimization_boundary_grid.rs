//! Example: Layout Optimization - Boundary Grid
//!
//! This example demonstrates boundary-constrained grid layout optimization.
//!
//! Corresponds to: examples_layout_optimization/005_layout_optimization_boundary_grid.py

use florus::Result;

fn main() -> Result<()> {
    println!("=== Layout Optimization - Boundary Grid ===\n");

    println!("Boundary-Constrained Grid Layout:\\n");
    
    println!("Concept:");
    println!("  - Define farm boundary polygon");
    println!("  - Create grid within boundary");
    println!("  - Place turbines at valid grid points");
    println!("  - Optimize subset selection\\n");

    println!("Boundary Definition:");
    println!("  - Polygon vertices (x, y coordinates)");
    println!("  - Can be convex or concave");
    println!("  - Multiple polygons (exclusion zones)");
    println!("  - Buffer zones from edges\\n");

    println!("Grid Generation:");
    println!("  1. Calculate bounding box of boundary");
    println!("  2. Create regular grid");
    println!("  3. Test each point for boundary inclusion");
    println!("  4. Keep only interior points");
    println!("  5. Apply minimum spacing filter\\n");

    println!("Optimization:");
    println!("  - Select N turbines from M grid points");
    println!("  - Maximize AEP subject to constraints");
    println!("  - Combinatorial problem");
    println!("  - Use greedy or genetic algorithms\\n");

    println!("Applications:");
    println!("  - Irregular land parcels");
    println!("  - Farms with exclusion zones");
    println!("  - Property boundary compliance");
    println!("  - Environmental restrictions\\n");

    println!("=== Example Complete ===");
    println!("Note: Requires boundary handling and grid generation.");
    Ok(())
}
