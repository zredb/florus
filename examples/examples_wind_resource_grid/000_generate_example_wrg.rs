/// Example: Generate Example WRG File
///
/// This example demonstrates how to create a Wind Resource Grid (WRG) file.
/// WRG files contain wind resource data at multiple grid points across a farm.
///
/// This is the Rust equivalent of Python's examples_wind_resource_grid/000_generate_example_wrg.py
///
/// Note: This is a placeholder showing the WRG concept.

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Generate Example WRG File");
    println!("============================================\n");

    println!("--- Wind Resource Grid (WRG) ---\n");
    
    println!("WRG files contain:");
    println!("  - Wind speed distributions at multiple locations");
    println!("  - Wind direction frequencies");
    println!("  - Turbulence intensity data");
    println!("  - Used for spatially-varying wind resources\n");
    
    println!("--- WRG File Format ---\n");
    
    println!("Top line: Nx Ny Xmin Ymin cell_size");
    println!("  Nx, Ny: Number of grid points in x and y");
    println!("  Xmin, Ymin: Starting coordinates");
    println!("  cell_size: Grid spacing\n");
    
    println!("Per grid point:");
    println!("  - x, y coordinates");
    println!("  - Weibull parameters (A, k)");
    println!("  - Mean wind direction");
    println!("  - Turbulence intensity\n");
    
    // ============================================================
    // Example WRG structure
    // ============================================================
    println!("--- Example WRG Structure ---\n");
    
    let nx = 3;
    let ny = 3;
    let xmin = 0.0;
    let ymin = 0.0;
    let cell_size = 1000.0;
    
    println!("Grid dimensions: {} x {}", nx, ny);
    println!("Domain: {}m x {}m", nx as f64 * cell_size, ny as f64 * cell_size);
    println!();
    
    println!("Example grid points:");
    println!("{:>6} {:>10} {:>10} {:>10} {:>10}", "i", "j", "x (m)", "y (m)", "A (m/s)");
    println!("{}", "-".repeat(50));
    
    // Generate example points
    for i in 0..nx {
        for j in 0..ny {
            let x = xmin + i as f64 * cell_size;
            let y = ymin + j as f64 * cell_size;
            let a = 8.0 + (i as f64 * 0.5); // Speed increases with x
            println!("{:>6} {:>10} {:>10.0} {:>10.0} {:>10.1}", i, j, x, y, a);
        }
    }
    
    println!("\n--- WRG for FLORIS ---\n");
    
    println!("Note: Full WRG support in FLORIS-RS is future work.");
    println!("Current workaround: Use WindRose with multiple sectors.");
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
