/// Example: Get Started with Flow Field Analysis
///
/// This example demonstrates basic flow field analysis in FLORIS-RS.
/// The flow field contains detailed velocity information at each grid point
/// which is useful for visualization and analysis.
///
/// This is the Rust equivalent of Python's examples_get_flow/001_get_started.py

use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Get Started with Flow Field Analysis");
    println!("========================================================\n");

    println!("--- Flow Field Overview ---\n");
    
    println!("Flow field analysis provides:");
    println!("  - Velocity at each grid point");
    println!("  - Rotor-averaged velocities");
    println!("  - Wake deficit visualization");
    println!("  - Power and thrust calculations\n");
    
    // ============================================================
    // Load model and set up
    // ============================================================
    let mut model = florus::FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Single turbine
    model.set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))?;
    
    // Set wind conditions
    model.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;
    
    // Run
    model.run()?;
    
    println!("--- Flow Field Results ---\n");
    
    // Get velocities (would need access to internal flow field)
    let powers = model.get_turbine_powers();
    let cts = model.get_turbine_thrust_coefficients();
    let ais = model.get_turbine_ais();
    
    println!("Turbine Power: {:.1} kW", powers[[0, 0]] / 1000.0);
    println!("Thrust Coefficient: {:.4}", cts[[0, 0]]);
    println!("Axial Induction: {:.4}", ais[[0, 0]]);
    
    // ============================================================
    // Multiple wind directions
    // ============================================================
    println!("\n--- Wind Direction Sweep ---\n");
    
    let wind_directions: Vec<f64> = (250..290).map(|d| d as f64).collect();
    let wind_speeds = Array1::from_vec(vec![8.0; wind_directions.len()]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; wind_directions.len()]);
    
    model.set_wind_conditions(
        wind_speeds,
        Array1::from_vec(wind_directions.clone()),
        turbulence_intensities,
    )?;
    
    model.run()?;
    
    let powers = model.get_turbine_powers();
    
    println!("{:>8} {:>12}", "WD (°)", "Power (kW)");
    println!("{}", "-".repeat(22));
    
    for (i, wd) in wind_directions.iter().enumerate().step_by(5) {
        println!("{:>8.0} {:>12.1}", wd, powers[[i, 0]] / 1000.0);
    }
    
    // ============================================================
    // Multiple wind speeds
    // ============================================================
    println!("\n--- Wind Speed Sweep ---\n");
    
    let wind_speeds: Vec<f64> = (4..20).map(|d| d as f64).collect();
    let wind_directions = Array1::from_vec(vec![270.0; wind_speeds.len()]);
    let turbulence_intensities = Array1::from_vec(vec![0.06; wind_speeds.len()]);
    
    model.set_wind_conditions(
        Array1::from_vec(wind_speeds.clone()),
        wind_directions,
        turbulence_intensities,
    )?;
    
    model.run()?;
    
    let powers = model.get_turbine_powers();
    let cts = model.get_turbine_thrust_coefficients();
    
    println!("{:>8} {:>12} {:>12}", "WS (m/s)", "Power (kW)", "Ct");
    println!("{}", "-".repeat(35));
    
    for (i, ws) in wind_speeds.iter().enumerate() {
        println!("{:>8.1} {:>12.1} {:>12.4}", ws, powers[[i, 0]] / 1000.0, cts[[i, 0]]);
    }
    
    println!("\n====================================");
    println!("Example completed successfully!");
    
    Ok(())
}
