//! Example: Wind Data Comparisons
//!
//! This example demonstrates different wind data objects and their usage.
//! Shows TimeSeries, WindRose, and how to use them with FlorisModel.
//!
//! Corresponds to: examples_wind_data/001_wind_data_comparisons.py

use florus::{FlorisModel, Result};
use florus::wind_data::TimeSeries;

fn main() -> Result<()> {
    println!("=== Wind Data Comparisons ===\n");
    println!("This example compares different wind data representations:\n");

    // Initialize the FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // Set a simple 2-turbine layout
    fmodel.set_layout(
        &ndarray::arr1(&[0.0, 500.0]),
        &ndarray::arr1(&[0.0, 0.0]),
    )?;

    // 1. TimeSeries - Single time series of wind conditions
    println!("1. TimeSeries Wind Data");
    println!("   ---------------------");
    
    let wind_directions = ndarray::arr1(&[270.0, 280.0, 290.0]);
    let wind_speeds = ndarray::arr1(&[8.0, 9.0, 10.0]);
    let turbulence_intensities = ndarray::arr1(&[0.06, 0.07, 0.08]);
    
    let time_series = TimeSeries::new(wind_directions, wind_speeds, turbulence_intensities)?;
    
    println!("   Conditions:");
    let n_conditions = time_series.wind_directions.len();
    for i in 0..n_conditions {
        println!("     [{:2}] WD={:.0}°, WS={:.1} m/s, TI={:.2}", 
                 i,
                 time_series.wind_directions[i],
                 time_series.wind_speeds[i],
                 time_series.turbulence_intensities[i]);
    }
    
    fmodel.set_wind_data(&time_series)?;
    fmodel.run()?;
    
    let powers = fmodel.get_turbine_powers();
    println!("\n   Farm Power by Condition:");
    for i in 0..n_conditions {
        let farm_power: f64 = powers.row(i).sum();
        println!("     [{:2}] {:.2} kW", i, farm_power / 1000.0);
    }

    // 2. Uniform wind conditions (single condition)
    println!("\n\n2. Single Wind Condition");
    println!("   ----------------------");
    
    fmodel.set_wind_conditions(
        ndarray::arr1(&[8.0]),
        ndarray::arr1(&[270.0]),
        ndarray::arr1(&[0.06]),
    )?;
    
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers();
    
    println!("   Condition: WD=270°, WS=8.0 m/s, TI=0.06");
    println!("   Turbine Powers:");
    for ti in 0..powers.ncols() {
        println!("     T{:2}: {:.2} kW", ti, powers[[0, ti]] / 1000.0);
    }
    let farm_power: f64 = powers.row(0).sum();
    println!("   Farm Total: {:.2} kW", farm_power / 1000.0);

    // 3. Multiple wind speeds at same direction
    println!("\n\n3. Wind Speed Sweep");
    println!("   -----------------");
    
    let wind_speeds_sweep = ndarray::arr1(&[4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]);
    let n_ws = wind_speeds_sweep.len();
    let wd_array = ndarray::Array1::from_vec(vec![270.0; n_ws]);
    let ti_array = ndarray::Array1::from_vec(vec![0.06; n_ws]);
    
    fmodel.set_wind_conditions(wind_speeds_sweep.clone(), wd_array, ti_array)?;
    fmodel.run()?;
    
    let powers = fmodel.get_turbine_powers();
    
    println!("   Wind Direction: 270° (constant)");
    println!("   Wind Speed vs Farm Power:");
    println!("   WS (m/s) | Farm Power (kW)");
    println!("   ---------|----------------");
    for i in 0..n_ws {
        let farm_power: f64 = powers.row(i).sum();
        println!("   {:8.1} | {:14.2}", wind_speeds_sweep[i], farm_power / 1000.0);
    }

    println!("\n=== Example Complete ===");
    println!("\nNote: WindRose implementation requires additional WindData trait implementations.");
    println!("TimeSeries is the primary wind data representation currently available.");

    Ok(())
}
