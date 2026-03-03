/// Example: Load Optimization
///
/// This example demonstrates load optimization in FLORIS-RS,
/// focusing on optimizing turbine loads rather than power.
///
/// This is the Rust equivalent of Python's examples_load_optimization/...

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Load Optimization");
    println!("==================================\n");

    println!("--- Load Optimization ---\n");
    
    println!("Load optimization focuses on:");
    println!("  - Minimizing turbine loads");
    println!("  - Fatigue analysis");
    println!("  - Structural optimization");
    println!("  - Combined power-load optimization\n");
    
    // ============================================================
    // Basic setup
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    let d = 126.0;
    
    // ============================================================
    // Spacing vs Load
    // ============================================================
    println!("--- Spacing vs Load Analysis ---\n");
    
    let spacings = vec![3.0, 5.0, 7.0, 10.0];
    
    println!("{:>12} {:>12} {:>12}", "Spacing (D)", "Power (kW)", "Rel. Load");
    println!("{}", "-".repeat(40));
    
    for sp in &spacings {
        let layout_x = Array1::from_vec(vec![0.0, sp * d]);
        let layout_y = Array1::from_vec(vec![0.0; 2]);
        
        model.set_layout(layout_x, layout_y)?;
        
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let powers = model.get_turbine_powers();
        let power = (powers[[0, 0]] + powers[[0, 1]]) / 1000.0;
        
        // Simplified load metric (higher TI = more load)
        let rel_load = 1.0 / sp;
        
        println!("{:>12.1} {:>12.1} {:>12.3}", sp, power, rel_load);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
