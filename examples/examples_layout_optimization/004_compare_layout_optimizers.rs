/// Example: Compare Layout Optimizers
///
/// This example compares different layout optimization approaches
/// available in FLORIS-RS.
///
/// This is the Rust equivalent of Python's examples_layout_optimization/004_compare_layout_optimizers.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Compare Layout Optimizers");
    println!("========================================\n");

    println!("--- Layout Optimization Methods ---\n");
    
    println!("FLORIS-RS provides multiple optimization methods:");
    println!("  1. Random Search: Fast, simple exploration");
    println!("  2. SciPy-based: Gradient-based optimization");
    println!("  3. PyOptSparse: Advanced constrained optimization\n");
    
    // ============================================================
    // Simple example
    // ============================================================
    println!("--- Simple Optimization Demo ---\n");
    
    let d = 126.0;
    let n_turbines = 4;
    
    println!("Optimizing {} turbine layout\n", n_turbines);
    
    // Test different spacings
    let spacings = vec![3.0, 5.0, 7.0, 10.0];
    
    println!("{:>12} {:>12} {:>12}", "Spacing (D)", "Farm (kW)", "AEP (GWh)");
    println!("{}", "-".repeat(40));
    
    for spacing in &spacings {
        // Create grid layout
        let n = ((n_turbines as f64).sqrt()) as usize;
        let mut x_vals = Vec::new();
        let mut y_vals = Vec::new();
        
        for i in 0..n {
            for j in 0..n {
                if x_vals.len() < n_turbines {
                    x_vals.push(i as f64 * spacing * d);
                    y_vals.push(j as f64 * spacing * d);
                }
            }
        }
        
        let layout_x = Array1::from_vec(x_vals);
        let layout_y = Array1::from_vec(y_vals);
        
        // Run model
        let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
        model.set_layout(layout_x, layout_y)?;
        
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let power = model.get_farm_power()[[0]] / 1000.0;
        let aep = power * 8760.0 / 1_000_000.0; // Approximate annual
        
        println!("{:>12.1} {:>12.1} {:>12.2}", spacing, power, aep);
    }
    
    println!("\n--- Optimization Methods Comparison ---\n");
    
    println!("Method          | Pros                    | Cons");
    println!("----------------|-------------------------|------------------");
    println!("Random Search   | Simple, parallelizable  | May miss optimum");
    println!("SciPy           | Gradient-based, fast    | Needs derivatives");
    println!("PyOptSparse     | Handles constraints     | Complex setup\n");
    
    println!("====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
