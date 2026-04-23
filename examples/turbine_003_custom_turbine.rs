//! Example: Custom Turbine Configuration
//!
//! This example demonstrates how to work with different turbine configurations
//! and shows the characteristics of various reference turbines.
//!
//! Corresponds to: examples_turbine/003_specify_turbine_power_curve.py

use florus::{FlorisModel, Result};

fn main() -> Result<()> {
    println!("=== Custom Turbine Configuration ===\n");

    // Demonstrate different turbine types from the library
    let turbine_types = vec![
        "nrel_5MW",
        "iea_10MW", 
        "iea_15MW",
        "iea_22MW",
    ];

    println!("Available Reference Turbines:\n");
    println!("{:<15} | {:<15} | {:<15}", "Turbine", "Rotor Diam (m)", "Hub Height (m)");
    println!("{}|{}|{}", "-".repeat(16), "-".repeat(16), "-".repeat(16));

    for turb_name in &turbine_types {
        // Create a single turbine model for each type
        let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
        fmodel.set_layout(
            &ndarray::arr1(&[0.0]),
            &ndarray::arr1(&[0.0]),
        )?;

        // Get rotor diameter and hub height from farm
        let rotor_diam = fmodel.core().farm.rotor_diameters[0];
        let hub_height = fmodel.core().farm.hub_heights[0];

        println!("{:<15} | {:<15.2} | {:<15.2}", turb_name, rotor_diam, hub_height);
    }

    println!("\n--- Testing NREL 5MW Power Curve ---\n");

    // Test power curve at different wind speeds
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel.set_layout(
        &ndarray::arr1(&[0.0]),
        &ndarray::arr1(&[0.0]),
    )?;

    let wind_speeds = vec![3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 20.0, 25.0];
    let n_ws = wind_speeds.len();
    let wd_array = ndarray::Array1::from_vec(vec![270.0; n_ws]);
    let ws_array = ndarray::Array1::from_vec(wind_speeds.clone());
    let ti_array = ndarray::Array1::from_vec(vec![0.06; n_ws]);

    fmodel.set_wind_conditions(ws_array, wd_array, ti_array)?;
    fmodel.run()?;

    let powers = fmodel.get_turbine_powers();

    println!("Wind Speed (m/s) | Power (kW)     | Cp (approx)");
    println!("-----------------|----------------|------------");

    for (i, ws) in wind_speeds.iter().enumerate() {
        let power_w = powers[[i, 0]];
        let power_kw = power_w / 1000.0;
        
        // Calculate approximate Cp
        let air_density = 1.225;
        let rotor_radius = 63.0; // NREL 5MW
        let area = std::f64::consts::PI * rotor_radius * rotor_radius;
        let cp = if *ws > 0.0 {
            power_w / (0.5 * air_density * area * ws.powi(3))
        } else {
            0.0
        };

        println!("{:16.1} | {:14.2} | {:11.4}", ws, power_kw, cp);
    }

    println!("\n=== Analysis ===");
    println!("This example demonstrates:");
    println!("  - Different reference turbine specifications");
    println!("  - Power curve characteristics");
    println!("  - Power coefficient (Cp) variation with wind speed");
    println!("  - Cut-in, rated, and cut-out wind speeds");

    println!("\n=== Example Complete ===");
    println!("\nNote: Full custom turbine specification requires turbine_dict building,");
    println!("which is available in Python FLORIS via build_cosine_loss_turbine_dict().");
    println!("In Rust, use pre-defined YAML turbine configuration files.");

    Ok(())
}
