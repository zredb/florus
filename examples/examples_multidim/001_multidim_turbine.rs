/// Example: Multi-Dimensional Turbine Parameters
///
/// This example demonstrates multi-dimensional turbine parameterization
/// in FLORIS-RS, including TI-dependent power/thrust curves.
///
/// This is the Rust equivalent of Python's examples_multidim/...

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Multi-Dimensional Turbine Parameters");
    println!("==================================================\n");

    println!("--- Multi-Dimensional Parameters ---\n");
    
    println!("Multi-dimensional turbine parameters allow:");
    println!("  - TI-dependent power curves");
    println!("  - TI-dependent thrust coefficients");
    println!("  - More accurate simulations");
    println!("  - Better load estimation\n");
    
    // ============================================================
    // Compare standard vs multi-dim
    // ============================================================
    println!("--- Standard vs Multi-Dim ---\n");
    
    // Standard turbine
    let mut standard_model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    standard_model.set_layout(
        &Array1::from_vec(vec![0.0]),
        &Array1::from_vec(vec![0.0]),
    )?;
    
    // Test different TI values
    let tis = vec![0.03, 0.06, 0.10, 0.15, 0.20];
    
    println!("{:>8} {:>16} {:>16}", "TI", "Standard (kW)", "Multi-Dim (kW)");
    println!("{}", "-".repeat(45));
    
    for ti in &tis {
        standard_model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![*ti]),
        )?;
        
        standard_model.run()?;
        
        let std_power = standard_model.get_turbine_powers()[[0, 0]] / 1000.0;
        
        // Try multi-dim if available
        let multi_power = if let Ok(mut mdm) = florus::FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct_TI.yaml") {
            mdm.set_layout(
                &Array1::from_vec(vec![0.0]),
                &Array1::from_vec(vec![0.0]),
            )?;
            mdm.set_wind_conditions(
                Array1::from_vec(vec![8.0]),
                Array1::from_vec(vec![270.0]),
                Array1::from_vec(vec![*ti]),
            )?;
            if mdm.run().is_ok() {
                Some(mdm.get_turbine_powers()[[0, 0]] / 1000.0)
            } else {
                None
            }
        } else {
            None
        };
        
        match multi_power {
            Some(mp) => println!("{:>8.2} {:>16.1} {:>16.1}", ti, std_power, mp),
            None => println!("{:>8.2} {:>16.1} {:>16}", ti, std_power, "N/A"),
        }
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
