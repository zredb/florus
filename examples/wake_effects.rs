/// Example: Wake Effects Demonstration
///
/// Demonstrates wake modeling with 2-turbine setup where the second turbine
/// is in the wake of the first turbine

use florus::core::{Farm, FlowField};
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Wake Effects Demonstration\n");
    
    // Create a 2-turbine wind farm layout
    // Turbine 0 at x=0 (upstream), Turbine 1 at x=630m (downstream)
    let layout_x = Array1::from_vec(vec![0.0, 630.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 2];
    
    println!("Creating 2-turbine wind farm:");
    println!("  Turbine 0: position (0, 0) m - UPSTREAM");
    println!("  Turbine 1: position (630, 0) m - DOWNSTREAM");
    
    let farm = Farm::new(layout_x, layout_y, turbine_types)?;
    
    // Set up flow conditions (single wind direction from west)
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]); // From West
    let turbulence_intensities = Array1::from_vec(vec![0.06]);
    
    println!("\nWind conditions:");
    println!("  Wind speed: 8.0 m/s");
    println!("  Wind direction: 270° (from West)");
    println!("  Turbulence intensity: 0.06");
    
    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions,
        0.0,    // wind_veer
        0.14,   // wind_shear
        1.225,  // air_density
        turbulence_intensities,
        90.0,   // reference_wind_height
    )?;
    
    // Create FlorisModel and run simulation
    let mut model = florus::FlorisModel {
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };
    
    // Initialize grid and flow field
    model.initialize_grid()?;
    model.initialize_flow_field()?;
    
    println!("\nGrid initialized:");
    let grid = model.grid.as_ref().unwrap();
    println!("  Number of turbines: {}", grid.n_turbines());
    println!("  Number of flow conditions: {}", grid.n_findex());
    
    // Get initial velocities (before wake)
    let shape = model.flow_field.u_initial_sorted.shape();
    println!("\nFlow field shape: {:?}", shape);
    
    // Run the wake solver
    println!("\nRunning wake solver...");
    model.run()?;
    
    // Calculate and display power for each turbine
    println!("\nResults:");
    println!("========");
    
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
        
        println!("  Turbine {}: velocity = {:.2} m/s, power = {:.0} kW", 
                 ti, vel, power / 1000.0);
    }
    
    // Calculate wake loss
    let upstream_vel = model.flow_field.u_sorted[[0, 0, 1, 1]];
    let downstream_vel = model.flow_field.u_sorted[[0, 1, 1, 1]];
    let wake_loss = (1.0 - downstream_vel / upstream_vel) * 100.0;
    
    println!("\nWake Analysis:");
    println!("  Upstream velocity: {:.2} m/s", upstream_vel);
    println!("  Downstream velocity: {:.2} m/s", downstream_vel);
    println!("  Wake loss: {:.1}%", wake_loss);
    
    if downstream_vel < upstream_vel * 0.95 {
        println!("\n✓ Wake effect successfully modeled!");
    } else {
        println!("\nNote: Wake effect may need solver completion");
    }
    
    println!("\nExample completed successfully!");
    
    Ok(())
}
