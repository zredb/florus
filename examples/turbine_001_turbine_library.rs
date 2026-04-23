//! Example: Reference Turbines
//!
//! This example demonstrates the reference turbines available in the turbine library.
//! For each turbine, it shows the power and thrust coefficient curves.
//!
//! Corresponds to: examples_turbine/001_reference_turbines.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Reference Turbines ===\n");

    // Initialize the FLORIS model with GCH
    let fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Get list of available turbine types from the library
    println!("Available turbine types in library:\n");
    
    // List common reference turbines
    let turbines = vec![
        "nrel_5MW",
        "iea_10MW",
        "iea_15MW",
        "iea_22MW",
    ];

    for (i, turbine_name) in turbines.iter().enumerate() {
        println!("{}. {}", i + 1, turbine_name);
    }

    println!("\n--- Demonstrating Default Turbine (NREL 5MW) ---\n");

    // Create a single turbine simulation
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0]),
        &ndarray::arr1(&[0.0]),
    )?;

    // Test at different wind speeds
    let wind_speeds = vec![3.0, 6.0, 9.0, 12.0, 15.0, 20.0, 25.0];
    let n_ws = wind_speeds.len();
    let wd_array = ndarray::Array1::from_vec(vec![270.0; n_ws]);
    let ws_array = ndarray::Array1::from_vec(wind_speeds.clone());
    let ti_array = ndarray::Array1::from_vec(vec![0.06; n_ws]);

    fmodel.set_wind_conditions(ws_array, wd_array, ti_array)?;

    println!("Wind Speed (m/s) | Power (kW)");
    println!("-----------------|------------");

    // Run simulation and get results
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers();
    
    for (i, ws) in wind_speeds.iter().enumerate() {
        let power_kw = powers[[i, 0]] / 1000.0;
        println!("{:16.1} | {:10.2}", ws, power_kw);
    }

    println!("\n--- Yaw Sensitivity Test ---\n");
    println!("Testing power loss due to yaw misalignment at 11 m/s:\n");

    // Set single wind speed for yaw test
    fmodel.set_wind_conditions(
        ndarray::arr1(&[11.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;

    let yaw_angles_test = vec![-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0];

    println!("Yaw Angle (deg) | Power (kW)    | Power Loss (%)");
    println!("----------------|---------------|----------------");

    // Get baseline power (0 deg yaw)
    fmodel.set_yaw_angles(ndarray::arr2(&[[0.0]]))?;
    fmodel.run()?;
    let baseline_power = fmodel.get_turbine_powers()[[0, 0]];

    for yaw in &yaw_angles_test {
        fmodel.set_yaw_angles(ndarray::arr2(&[[*yaw]]))?;
        fmodel.run()?;
        let power = fmodel.get_turbine_powers()[[0, 0]];
        let power_loss_pct = (1.0 - power / baseline_power) * 100.0;
        
        println!("{:15.1} | {:13.2} | {:14.2}", 
                 yaw, power / 1000.0, power_loss_pct);
    }

    println!("\n=== Example Complete ===");
    println!("\nNote: Full visualization with plots requires plotting library integration.");
    println!("This example demonstrates the data that would be plotted in Python version.");

    Ok(())
}
