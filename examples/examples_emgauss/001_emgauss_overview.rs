/// Example: Empirical Gauss Model Overview
///
/// This example provides an overview of the Empirical Gauss wake model
/// in FLORIS-RS.
///
/// The Empirical Gauss model combines Gaussian wake profiles with
/// empirical parameters for improved accuracy.
///
/// This is the Rust equivalent of Python's examples_emgauss/001_emgauss_overview.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Empirical Gauss Model Overview");
    println!("==============================================\n");

    println!("--- Empirical Gauss Model ---\n");
    
    println!("The Empirical Gauss model:");
    println!("  - Uses Gaussian wake deficit profiles");
    println!("  - Applies empirical expansion rates");
    println!("  - Good for offshore conditions");
    println!("  - Includes tilt-driven vertical deflection\n");
    
    // ============================================================
    // Compare with GCH
    // ============================================================
    println!("--- Model Comparison ---\n");
    
    let configs = vec![
        ("examples/inputs/gch.yaml", "Gauss-Curl-Hybrid"),
        ("examples/inputs/emgauss.yaml", "Empirical Gauss"),
    ];
    
    let wind_speeds = vec![5.0, 8.0, 10.0, 12.0, 15.0];
    
    println!("{:>8}", "WS (m/s)");
    for (_, name) in &configs {
        print!(" {:>16}", name);
    }
    println!();
    println!("{}", "-".repeat(8 + 18 * configs.len() + 1));
    
    for ws in &wind_speeds {
        print!("{:>8.1}", ws);
        
        for (config, _) in &configs {
            let mut model = florus::FlorisModel::from_file(config)?;
            
            model.set_layout(
                &Array1::from_vec(vec![0.0]),
                &Array1::from_vec(vec![0.0]),
            )?;
            
            model.set_wind_conditions(
                Array1::from_vec(vec![*ws]),
                Array1::from_vec(vec![270.0]),
                Array1::from_vec(vec![0.06]),
            )?;
            
            model.run()?;
            
            let power = model.get_turbine_powers()[[0, 0]] / 1000.0;
            print!(" {:>16.1}", power);
        }
        println!();
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
