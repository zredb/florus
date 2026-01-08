/// Example: Multi-Turbine Wind Farm
///
/// Demonstrates wake interactions in a larger wind farm with 5 turbines

use florus::core::Farm;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Multi-Turbine Wind Farm\n");
    
    // Create a 5-turbine wind farm in a line
    // Spacing: 7D (7 * 126m = 882m) is typical for offshore wind farms
    let rotor_diameter = 126.0; // NREL 5MW
    let spacing = 7.0 * rotor_diameter; // 882m spacing
    
    let layout_x = Array1::from_vec(vec![
        0.0,
        spacing,
        2.0 * spacing,
        3.0 * spacing,
        4.0 * spacing,
    ]);
    let layout_y = Array1::from_vec(vec![0.0; 5]);
    let turbine_types = vec!["nrel_5MW".to_string(); 5];
    
    println!("Creating 5-turbine wind farm (linear layout):");
    for (i, x) in layout_x.iter().enumerate() {
        println!("  Turbine {}: position ({:.0}, {:.0}) m", i, x, layout_y[i]);
    }
    println!("  Spacing: {:.0} m ({:.1}D)", spacing, spacing / rotor_diameter);
    
    let farm = Farm::new(layout_x, layout_y, turbine_types)?;
    
    // Set up flow conditions
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]); // From West
    let turbulence_intensities = Array1::from_vec(vec![0.06]);
    
    println!("\nWind conditions:");
    println!("  Wind speed: 8.0 m/s");
    println!("  Wind direction: 270° (from West)");
    println!("  Turbulence intensity: 0.06");
    
    let flow_field = florus::core::FlowField::new(
        wind_speeds.clone(),
        wind_directions,
        0.0,    // wind_veer
        0.14,   // wind_shear
        1.225,  // air_density
        turbulence_intensities,
        90.0,   // reference_wind_height
    )?;
    
    // Create and run model
    let mut model = florus::FlorisModel {
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };
    
    model.initialize_grid()?;
    model.initialize_flow_field()?;
    model.run()?;
    
    println!("\nResults:");
    println!("========");
    
    let mut total_power = 0.0;
    for ti in 0..model.farm.n_turbines() {
        let vel = model.flow_field.u_sorted[[0, ti, 1, 1]];
        let rotor_d = model.farm.rotor_diameters[ti];
        let area = std::f64::consts::PI * (rotor_d / 2.0).powi(2);
        let turbine = &model.farm.turbine_map[ti];
        
        let power = if vel < turbine.cut_in_wind_speed || vel > turbine.cut_out_wind_speed {
            0.0
        } else {
            let cp = turbine.power_coefficient(vel);
            0.5 * model.flow_field.air_density * area * vel.powi(3) * cp
        }.min(turbine.rated_power);
        
        total_power += power;
        println!("  Turbine {}: velocity = {:.2} m/s, power = {:.0} kW", 
                 ti, vel, power / 1000.0);
    }
    
    // Calculate wake loss relative to first turbine
    let ref_power = model.farm.turbine_map[0].rated_power;
    println!("\nFarm Summary:");
    println!("  Total power: {:.1} MW", total_power / 1_000_000.0);
    println!("  Capacity factor: {:.1}%", total_power / (5.0 * ref_power) * 100.0);
    
    // Show wake cascade
    println!("\nWake Cascade Analysis:");
    let v0 = model.flow_field.u_sorted[[0, 0, 1, 1]];
    for ti in 1..model.farm.n_turbines() {
        let vel = model.flow_field.u_sorted[[0, ti, 1, 1]];
        let loss = (1.0 - vel / v0) * 100.0;
        println!("  Turbine {}: {:.1}% loss from upstream turbines", ti, loss);
    }
    
    println!("\nExample completed successfully!");
    
    Ok(())
}
