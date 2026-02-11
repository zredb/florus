/// Generate Gridded Layout Example
///
/// This example demonstrates how to generate a gridded wind farm layout
/// using meshgrid-style coordinate generation.
///
/// This is the Rust equivalent of Python's 004_generate_gridded_layout.py

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS: Generate Gridded Layout");
    println!("================================\n");

    // ============================================================
    // Grid Parameters
    // ============================================================
    println!("--- Grid Parameters ---\n");

    let d = 126.0; // NREL 5MW rotor diameter
    let n_rows = 5;
    let n_cols = 4;
    let spacing_x = 5.0 * d; // 630 m
    let spacing_y = 7.0 * d; // 882 m

    println!("Grid configuration:");
    println!("  Rows: {}", n_rows);
    println!("  Columns: {}", n_cols);
    println!("  Total turbines: {}", n_rows * n_cols);
    println!());

    println!("Spacing:");
    println!("  X spacing: {:.0} m ({:.1}D)", spacing_x, spacing_x / d);
    println!("  Y spacing: {:.0} m ({:.1}D)", spacing_y, spacing_y / d);
    println!());

    // ============================================================
    // Generate Grid
    // ============================================================
    println!("--- Generating Grid ---\n");

    // Generate x coordinates for each column
    // In Python: x = spacing_x * np.arange(0, n_cols, 1)
    let x_coords: Vec<f64> = (0..n_cols).map(|i| spacing_x * (i as f64)).collect();

    // Generate y coordinates for each row
    // In Python: y = spacing_y * np.arange(0, n_rows, 1)
    let y_coords: Vec<f64> = (0..n_rows).map(|i| spacing_y * (i as f64)).collect();

    println!("X coordinates (columns):");
    for (i, &x) in x_coords.iter().enumerate() {
        println!("  Col {}: {:.0} m", i, x);
    }
    println!());

    println!("Y coordinates (rows):");
    for (i, &y) in y_coords.iter().enumerate() {
        println!("  Row {}: {:.0} m", i, y);
    }
    println!());

    // Create meshgrid-like layout
    // In Python: X, Y = np.meshgrid(x, y); layout_x = X.flatten(); layout_y = Y.flatten()
    let mut layout_x = Vec::with_capacity(n_rows * n_cols);
    let mut layout_y = Vec::with_capacity(n_rows * n_cols);

    for &y in &y_coords {
        for &x in &x_coords {
            layout_x.push(x);
            layout_y.push(y);
        }
    }

    let turbine_types = vec!["nrel_5MW".to_string(); n_rows * n_cols];

    println!("Combined layout (row-major order):");
    for (i, (x, y)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
        let row = i / n_cols;
        let col = i % n_cols;
        println!("  Turbine {} (R{}, C{}): ({:.0}, {:.0})", i, row, col, x, y);
    }

    // ============================================================
    // Create Farm
    // ============================================================
    println!("\n--- Creating Farm ---\n");

    let farm = Farm::new(
        Array1::from_vec(layout_x.clone()),
        Array1::from_vec(layout_y.clone()),
        turbine_types.clone(),
    )?;

    println!("Farm created with {} turbines", farm.n_turbines());
    println!("Grid area: {:.0} x {:.0} m", spacing_x * (n_cols - 1), spacing_y * (n_rows - 1));

    // ============================================================
    // Alternative Grid Configurations
    // ============================================================
    println!("\n--- Alternative Configurations ---\n");

    println!("1. Staggered grid (hexagonal pattern):");
    let spacing_x_s = 5.0 * d;
    let spacing_y_s = 3.5 * d;
    println!("   Offset odd rows by half X spacing");
    println!("   Spacing: {:.0} x {:.0} m", spacing_x_s, spacing_y_s);
    println!());

    println!("2. Circular grid:");
    println!("   Turbines arranged in concentric circles");
    println!("   Useful for single-location studies");
    println!());

    println!("3. Random within boundary:");
    println!("   Random turbine placement within defined polygon");
    println!("   Requires minimum spacing constraints");
    println!());

    // ============================================================
    // Summary
    // ============================================================
    println!("\n--- Summary ---\n");

    println!("Gridded Layout Key Points:");
    println!("  ✓ Meshgrid approach for systematic layouts");
    println!("  ✓ Easy to control row/column count");
    println!("  ✓ Consistent spacing in both directions");
    println!("  ✓ Can be staggered for hexagonal patterns");
    println!("  ✓ Foundation for optimization studies");
    println!());

    println!("Applications:");
    println!("  - Baseline farm design comparisons");
    println!("  - Spacing sensitivity studies");
    println!("  - Grid convergence analysis");
    println!("  - Rapid layout prototyping");

    println!("\n================================");
    println!("Example completed successfully!");

    Ok(())
}
