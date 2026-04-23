//! Example: Extract Turbulence Intensity at Points
//!
//! This example demonstrates the use of sample_ti_at_points to extract
//! turbulence intensity (TI) information at user-specified locations.
//!
//! Specifically, this example returns the TI at a single x, y location
//! and four different heights over a sweep of wind directions.
//! This mimics the TI measurements of a met mast across all wind directions.
//!
//! Corresponds to: examples_get_flow/003_extract_turbulence_intensity_at_points.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Extract Turbulence Intensity at Points ===\n");

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

    // Sample TI at different wind directions
    let wd_array = vec![0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0];
    let ws = 8.0;
    let ti_inflow = 0.06;
    
    println!("Turbulence Intensity at Met Mast (sample at every 30 degrees):");
    println!("WD (deg) | z=30m  | z=90m  | z=150m | z=250m");
    println!("---------|--------|--------|--------|--------");
    
    for wd in &wd_array {
        // Set single wind condition
        fmodel.set_wind_conditions(
            ndarray::arr1(&[ws]),
            ndarray::arr1(&[*wd]),
            ndarray::arr1(&[ti_inflow]),
        )?;
        fmodel.run()?;
        
        // Sample TI at this wind direction
        let ti_at_points = fmodel.sample_ti_at_points(&points_x, &points_y, &points_z)?;
        
        print!("{:8.0} |", wd);
        for z_idx in 0..n_points {
            print!(" {:6.4} |", ti_at_points[z_idx]);
        }
        println!();
    }
    
    println!("\n=== Analysis ===");
    println!("Turbulence intensity increases in turbine wakes.");
    println!("As wind direction changes, the met mast experiences:");
    println!("  - Ambient TI when not in wake (~0.06)");
    println!("  - Elevated TI when in turbine wake (>0.06)");
    println!("  - Different TI levels at different heights");

    println!("\n=== Example Complete ===");
    println!("\nNote: Full visualization would plot TI vs wind direction");
    println!("for each height, showing wake-induced turbulence.");

    Ok(())
}
