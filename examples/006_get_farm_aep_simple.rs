/// Example 6: Getting Expected Power and AEP (Simplified)
///
/// This example demonstrates how to calculate expected farm power and AEP
/// using frequency-weighted averaging.

use florus::{Array1, FlorisModel};
use florus::core::turbines::TurbineLibrary;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Set to a 3-turbine layout
    let d = 126.0;
    fmodel.set_layout(
        &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
        &Array1::from_vec(vec![0.0, 0.0, 0.0]),
    )?;

    println!("Layout: 3 turbines in a row");
    println!("  Positions: [0, {:.0}, {:.0}] m", 5.0 * d, 10.0 * d);

    // Create a simple set of conditions with different wind speeds
    let wind_speeds = vec![8.0, 9.0, 10.0];
    let wind_directions = vec![270.0, 270.0, 270.0];
    let turbulence_intensities = vec![0.06, 0.06, 0.06];
    let n_conditions = wind_speeds.len();

    fmodel.set_wind_conditions(
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(turbulence_intensities),
    )?;

    println!("\nWind conditions: {} scenarios", n_conditions);
    for i in 0..n_conditions {
        println!("  Condition {}: {} m/s at {}°", i + 1, wind_speeds[i], wind_directions[i]);
    }

    // Run the model
    fmodel.run()?;

    println!("\n========== Calculating Expected Power and AEP ==========");

    // With uniform frequency, each condition has equal weight
    let uniform_freq = vec![1.0 / n_conditions as f64; n_conditions];
    
    let expected_farm_power = fmodel.get_expected_farm_power(
        Some(Array1::from_vec(uniform_freq.clone())),
        None,
    )?;
    
    // Calculate AEP: expected power (kW) * hours per year
    let hours_per_year = 365.0 * 24.0;
    let aep = expected_farm_power * hours_per_year;

    println!("\nUniform Frequency Case:");
    println!("  Expected farm power: {:.2} kW", expected_farm_power / 1000.0);
    println!("  AEP: {:.2} MWh/year", aep / 1e6);
    println!("  AEP: {:.3} GWh/year", aep / 1e9);

    // Now try with non-uniform frequencies
    println!("\nNon-Uniform Frequency Case:");
    // Give more weight to higher wind speeds
    let custom_freq = vec![0.2, 0.3, 0.5]; // 20%, 30%, 50%
    
    let expected_farm_power_custom = fmodel.get_expected_farm_power(
        Some(Array1::from_vec(custom_freq.clone())),
        None,
    )?;
    
    let aep_custom = expected_farm_power_custom * hours_per_year;

    println!("  Frequencies: {:?}", custom_freq);
    println!("  Expected farm power: {:.2} kW", expected_farm_power_custom / 1000.0);
    println!("  AEP: {:.2} MWh/year", aep_custom / 1e6);
    println!("  AEP: {:.3} GWh/year", aep_custom / 1e9);

    // Compare wake vs no-wake scenarios
    println!("\n========== Wake Loss Analysis ==========");
    
    // Run without wakes
    fmodel.run_no_wake()?;
    
    let expected_farm_power_no_wake = fmodel.get_expected_farm_power(
        Some(Array1::from_vec(uniform_freq)),
        None,
    )?;
    
    let aep_no_wake = expected_farm_power_no_wake * hours_per_year;
    let wake_losses = 100.0 * (aep_no_wake - aep) / aep_no_wake;

    println!("  AEP with wakes: {:.3} GWh/year", aep / 1e9);
    println!("  AEP without wakes: {:.3} GWh/year", aep_no_wake / 1e9);
    println!("  Wake losses: {:.2}%", wake_losses);

    println!("\nExample 6 completed successfully!");

    Ok(())
}
