/// Example: Basic FLORIS usage

use florus::core::{Farm, FlowField};
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    println!("FLORIS-RS Example: Basic Wind Farm Simulation\n");
    
    // Create a simple 3-turbine wind farm layout
    let layout_x = Array1::from_vec(vec![0.0, 630.0, 1260.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 3];
    
    println!("Creating wind farm with {} turbines", layout_x.len());
    let mut farm = Farm::new(layout_x, layout_y, turbine_types)?;
    
    // Set up flow conditions
    let wind_speeds = Array1::from_vec(vec![8.0, 10.0, 12.0, 14.0]);
    let wind_directions = Array1::from_vec(vec![270.0, 270.0, 270.0, 270.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06, 0.06, 0.06, 0.06]);
    
    println!("\nWind conditions:");
    println!("  Wind speeds: {:?} m/s", wind_speeds.as_slice().unwrap());
    println!("  Wind directions: {:?} degrees", wind_directions.as_slice().unwrap());
    println!("  Turbulence intensities: {:?}", turbulence_intensities.as_slice().unwrap());
    
    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions,
        0.0,    // wind_veer
        0.14,   // wind_shear
        1.225,  // air_density
        turbulence_intensities,
        90.0,   // reference_wind_height
    )?;
    
    // Initialize control arrays
    farm.initialize_control_arrays(flow_field.n_findex);
    
    // Demonstrate coordinate calculations
    println!("\nTurbine coordinates:");
    let coords = farm.coordinates();
    for i in 0..farm.n_turbines() {
        println!("  Turbine {}: ({:.1}, {:.1}) m", i, coords[[i, 0]], coords[[i, 1]]);
    }
    
    // Calculate turbine distances
    println!("\nInter-turbine distances:");
    let n = farm.n_turbines();
    for i in 0..n {
        for j in (i+1)..n {
            let dx = coords[[i, 0]] - coords[[j, 0]];
            let dy = coords[[i, 1]] - coords[[j, 1]];
            let dist = (dx * dx + dy * dy).sqrt();
            println!("  T{} to T{}: {:.1} m", i, j, dist);
        }
    }
    
    // Demonstrate wind speed at different heights
    println!("\nWind speed profile (at 8 m/s reference):");
    for height in [50.0, 90.0, 120.0, 150.0] {
        let ws = flow_field.wind_speed_at_height(height, 0);
        println!("  At {:.0} m: {:.2} m/s", height, ws);
    }
    
    println!("\nExample completed successfully!");
    
    Ok(())
}
