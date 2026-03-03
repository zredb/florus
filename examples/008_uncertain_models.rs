/// Example: Uncertain Models
///
/// This example demonstrates using uncertainty quantification in FLORIS-RS,
/// including handling uncertainty in wind direction and other parameters.
///
/// This is the Rust equivalent of Python's 008_uncertain_models.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example 008: Uncertain Models");
    println!("====================================\n");

    println!("--- Uncertainty in FLORIS ---\n");
    
    println!("FLORIS supports uncertainty quantification:");
    println!("  - Wind direction uncertainty");
    println!("  - Turbulence intensity variation");
    println!("  - Wake model uncertainty");
    println!("  - Parameter sensitivity analysis\n");
    
    // ============================================================
    // Basic model
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    let d = 126.0;
    let spacing = 5.0 * d;
    
    // 2-turbine layout
    model.set_layout(
        &Array1::from_vec(vec![0.0, spacing]),
        &Array1::from_vec(vec![0.0; 2]),
    )?;
    
    // ============================================================
    // Single condition (deterministic)
    // ============================================================
    println!("--- Deterministic Simulation ---\n");
    
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    
    model.run()?;
    
    let powers = model.get_turbine_powers();
    
    println!("Single condition at 270°:");
    println!("  Upstream: {:.1} kW", powers[[0, 0]] / 1000.0);
    println!("  Downstream: {:.1} kW", powers[[0, 1]] / 1000.0);
    println!("  Wake loss: {:.1}%\n", (1.0 - powers[[0, 1]] / powers[[0, 0]]) * 100.0);
    
    // ============================================================
    // Wind direction sweep (simulating uncertainty)
    // ============================================================
    println!("--- Wind Direction Uncertainty ---\n");
    
    // Simulate wind direction distribution around 270°
    let base_direction = 270.0;
    let std_dev = 10.0; // Standard deviation in degrees
    
    println!("Simulating uncertainty around {}° (std: {}°)\n", base_direction, std_dev);
    
    // Generate sample directions (in practice, would use statistical distribution)
    let directions: Vec<f64> = vec![260.0, 265.0, 270.0, 275.0, 280.0];
    let weights: Vec<f64> = vec![0.1, 0.2, 0.4, 0.2, 0.1]; // Weights summing to 1
    
    let mut weighted_power = 0.0;
    
    println!("{:>8} {:>8} {:>12}", "WD (°)", "Weight", "Farm (kW)");
    println!("{}", "-".repeat(32));
    
    for (i, wd) in directions.iter().enumerate() {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![*wd]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let farm_power = model.get_farm_power()[[0]] / 1000.0;
        println!("{:>8.0} {:>8.2} {:>12.1}", wd, weights[i], farm_power);
        
        weighted_power += farm_power * weights[i];
    }
    
    println!("  ----------------");
    println!("  Weighted avg: {:.1} kW\n", weighted_power);
    
    // ============================================================
    // Turbulence intensity sensitivity
    // ============================================================
    println!("--- Turbulence Intensity Sensitivity ---\n");
    
    let tis = vec![0.04, 0.06, 0.08, 0.10, 0.15];
    
    println!("{:>8} {:>12} {:>12}", "TI", "Up (kW)", "Down (kW)");
    println!("{}", "-".repeat(35));
    
    for ti in &tis {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![*ti]),
        )?;
        
        model.run()?;
        
        let powers = model.get_turbine_powers();
        println!("{:>8.2} {:>12.1} {:>12.1}", 
            ti, powers[[0, 0]] / 1000.0, powers[[0, 1]] / 1000.0);
    }
    
    // ============================================================
    // Power uncertainty estimation
    // ============================================================
    println!("\n--- Power Uncertainty Estimation ---\n");
    
    // Calculate mean and std of power
    let mut powers: Vec<f64> = Vec::new();
    
    for wd in 260..=280 {
        model.set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![wd as f64]),
            Array1::from_vec(vec![0.06]),
        )?;
        
        model.run()?;
        
        let farm_power = model.get_farm_power()[[0]];
        powers.push(farm_power);
    }
    
    let n = powers.len() as f64;
    let mean: f64 = powers.iter().sum::<f64>() / n;
    let variance: f64 = powers.iter().map(|p| (p - mean).powi(2)).sum::<f64>() / n;
    let std = variance.sqrt();
    
    println!("Power statistics over 260-280°:");
    println!("  Mean: {:.1} kW", mean / 1000.0);
    println!("  Std: {:.1} kW", std / 1000.0);
    println!("  CV: {:.1}%\n", std / mean * 100.0);
    
    println!("====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
