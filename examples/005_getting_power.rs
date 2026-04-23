/// Example 5: Getting Turbine and Farm Power
///
/// After setting the FlorisModel and running, the next step is typically to get the power output
/// of the turbines. FLORIS has several methods for getting power:
///
/// 1. `get_turbine_powers()`: Returns the power output of each turbine in the farm for each findex
///    (n_findex, n_turbines)
/// 2. `get_farm_power()`: Returns the total power output of the farm for each findex (n_findex)
///
/// This example demonstrates these methods using different wind conditions.

use florus::{Array1, FlorisModel};
use florus::core::turbines::TurbineLibrary;
use ndarray::Axis;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set to a 3-turbine layout
    let d = 126.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    println!("Layout: 3 turbines in a row (spacing: 5D, 10D)");
    println!("  n_turbines: {}", fmodel.n_turbines());

    //////////////////////////////////////////////////
    // Using TimeSeries-like conditions
    //////////////////////////////////////////////////

    // Set up conditions where wind direction sweeps from 250 to 290 degrees
    let wind_directions: Vec<f64> = (250..290).map(|x| x as f64).collect();
    let n_conditions = wind_directions.len();
    let wind_speeds = vec![9.9; n_conditions];
    let turbulence_intensities = vec![0.06; n_conditions];

    fmodel.set_wind_conditions(
        Array1::from_vec(wind_speeds),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(turbulence_intensities),
    )?;

    println!("\n========== TimeSeries-like Conditions ==========");
    println!("Number of conditions: {}", n_conditions);
    println!("Wind directions: {}° to {}°", wind_directions[0], wind_directions[n_conditions - 1]);

    // Run the model
    fmodel.run()?;

    // Get the turbine powers
    let turbine_powers = fmodel.get_turbine_powers();

    // Turbines powers will have shape (n_findex, n_turbines) where n_findex is the number of unique
    // wind conditions and n_turbines is the number of turbines in the farm
    println!("Turbine power has shape {:?}", turbine_powers.shape());

    // It is also possible to get the farm power directly
    let farm_power = fmodel.get_farm_power();

    // Farm power has length n_findex, and is the sum of the turbine powers
    println!("Farm power has shape {:?}", farm_power.shape());

    // Print some sample results
    println!("\nSample results (first 5 conditions):");
    for i in 0..5.min(n_conditions) {
        print!("  WD={:.0}°: ", wind_directions[i]);
        for ti in 0..turbine_powers.shape()[1] {
            print!("T{}={:.0}kW ", ti + 1, turbine_powers[[i, ti]] / 1000.0);
        }
        println!("| Farm={:.0}kW", farm_power[i] / 1000.0);
    }

    // It's possible to get these powers with wake losses disabled
    fmodel.run_no_wake()?;
    let farm_power_no_wake = fmodel.get_farm_power();

    // Calculate and print wake losses for a sample condition
    let sample_idx = 20; // Middle of the sweep
    if sample_idx < n_conditions {
        let wake_loss_pct = 100.0 * (farm_power_no_wake[sample_idx] - farm_power[sample_idx]) 
            / farm_power_no_wake[sample_idx];
        println!("\nWake loss at WD={}°: {:.2}%", 
                 wind_directions[sample_idx], wake_loss_pct);
    }

    //////////////////////////////////////////////////
    // Using WindRose-like conditions
    //////////////////////////////////////////////////

    println!("\n========== WindRose-like Conditions ==========");
    
    // Create a new model to demonstrate WindRose conditions
    let mut fmodel_rose = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set to a 3-turbine layout
    let d = 126.0;
    fmodel_rose.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;
    
    // Declare conditions for 2 wind directions and 3 wind speeds
    let wd_values = vec![270.0, 280.0];
    let ws_values = vec![8.0, 9.0, 10.0];
    
    let mut rose_wind_directions = Vec::new();
    let mut rose_wind_speeds = Vec::new();
    let mut rose_tis = Vec::new();
    
    for &wd in &wd_values {
        for &ws in &ws_values {
            rose_wind_directions.push(wd);
            rose_wind_speeds.push(ws);
            rose_tis.push(0.06);
        }
    }

    fmodel_rose.set_wind_conditions_with_rose(
        Array1::from_vec(rose_wind_speeds),
        Array1::from_vec(rose_wind_directions),
        Array1::from_vec(rose_tis),
        wd_values.clone(),
        ws_values.clone(),
    )?;

    println!("Wind directions: {:?}", wd_values);
    println!("Wind speeds: {:?}", ws_values);
    println!("Number of conditions (2 × 3): {}", fmodel_rose.n_findex());

    fmodel_rose.run()?;

    // Use reshaped methods to match Python output format
    let turbine_powers = fmodel_rose.get_turbine_powers_rose();

    println!("Shape of turbine powers: {:?}", turbine_powers.shape());

    let farm_power = fmodel_rose.get_farm_power_rose();

    println!("Shape of farm power: {:?}", farm_power.shape());

    // Print results organized by wind direction and speed
    println!("\nFarm power by condition:");
    for wd_idx in 0..wd_values.len() {
        for ws_idx in 0..ws_values.len() {
            print!("  WD={:.0}°: ", wd_values[wd_idx]);
            println!("WS={:.0}m/s → Farm Power={:.0}kW", 
                     ws_values[ws_idx], farm_power[[wd_idx, ws_idx]] / 1000.0);
        }
    }

    println!("\nExample 5 completed successfully!");

    Ok(())
}
