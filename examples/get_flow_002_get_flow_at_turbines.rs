//! Example: Extract Wind Speed at Points
//!
//! This example demonstrates the use of sample_flow_at_points to extract
//! wind speed information at user-specified locations in the flow.
//!
//! Specifically, this example returns the wind speed at a single x, y
//! location and four different heights over a sweep of wind directions.
//! This mimics the wind speed measurements of a met mast across all
//! wind directions (at a fixed free stream wind speed).
//!
//! Corresponds to: examples_get_flow/002_extract_wind_speed_at_points.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Extract Wind Speed at Points ===\n");

    // Instantiate FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set up a two-turbine farm
    let d = 126.0;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 3.0 * d]),
        &ndarray::arr1(&[0.0, 3.0 * d]),
    )?;

    println!("Turbine Layout:");
    println!("  T0: (0.0, 0.0)");
    println!("  T1: ({:.0}, {:.0})\n", 3.0 * d, 3.0 * d);

    // Set the wind direction to run 360 degrees
    let wd_array = ndarray::Array1::linspace(0.0, 359.0, 360);
    let n_wd = wd_array.len();
    let ws_array = ndarray::Array1::from_vec(vec![8.0; n_wd]);
    let ti_array = ndarray::Array1::from_vec(vec![0.06; n_wd]);

    fmodel.set_wind_conditions(ws_array.clone(), wd_array.clone(), ti_array.clone())?;

    // Run simulation first
    println!("Running simulation with 360 wind directions...\n");
    fmodel.run()?;

    // Simulate a met mast between the turbines
    let met_mast_option = 0; // Try 0, 1, 2, 3
    
    let (points_x, points_y) = match met_mast_option {
        0 => (vec![3.0 * d; 4], vec![0.0; 4]),
        1 => (vec![200.0; 4], vec![200.0; 4]),
        2 => (vec![20.0; 4], vec![20.0; 4]),
        3 => (vec![305.0; 4], vec![158.0; 4]),
        _ => panic!("Invalid met_mast_option"),
    };
    
    let points_z = vec![30.0, 90.0, 150.0, 250.0];
    let n_points = points_z.len();
    
    println!("Met Mast Location (option {}):", met_mast_option);
    println!("  Position: ({:.0}, {:.0})", points_x[0], points_y[0]);
    println!("  Heights: {:?} m\n", points_z);

    // Note: sample_flow_at_points samples at first findex only
    // To get all wind directions, we loop through them
    println!("Wind Speed at Met Mast (sample at every 30 degrees):");
    println!("WD (deg) | z=30m  | z=90m  | z=150m | z=250m");
    println!("---------|--------|--------|--------|--------");
    
    for i in (0..n_wd).step_by(30) {
        // Set single wind condition
        fmodel.set_wind_conditions(
            ndarray::arr1(&[ws_array[i]]),
            ndarray::arr1(&[wd_array[i]]),
            ndarray::arr1(&[ti_array[i]]),
        )?;
        fmodel.run()?;
        
        // Sample at this wind direction
        let u_at_points = fmodel.sample_flow_at_points(&points_x, &points_y, &points_z)?;
        
        print!("{:8.0} |", wd_array[i]);
        for z_idx in 0..n_points {
            print!(" {:6.2} |", u_at_points[z_idx]);
        }
        println!();
    }
    
    println!("\n=== Analysis ===");
    println!("The wind speed varies with height due to wind shear.");
    println!("As wind direction changes, the met mast experiences:");
    println!("  - Free stream flow when not in wake");
    println!("  - Reduced velocity when in turbine wake");
    println!("  - Different wake effects at different heights");

    println!("\n=== Example Complete ===");
    println!("\nNote: Full visualization would plot wind speed vs wind direction");
    println!("for each height, showing wake effects as wind direction changes.");

    Ok(())
}
