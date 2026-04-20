/// Example 6: Getting Expected Power and AEP
///
/// The expected power of a farm is computed by multiplying the power output of the farm by the
/// frequency of each findex. This is done by the `get_expected_farm_power` method. The expected
/// AEP (Annual Energy Production) is computed by multiplying the expected power by the number of
/// hours in a year.
///
/// If wind data with frequencies is provided to the model, the expected power and AEP
/// can be computed directly. If not, a frequency table must be passed into these functions.

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
    println!("  Positions: [0, {:.*}, {:.*}] m", 1, 5.0 * d, 1, 10.0 * d);

    //////////////////////////////////////////////////
    // Using uniform frequency (TimeSeries-like)
    //////////////////////////////////////////////////

    // Create a simple set of conditions
    let wind_speeds = vec![8.0, 9.0, 10.0];
    let wind_directions = vec![270.0, 270.0, 270.0];
    let turbulence_intensities = vec![0.06, 0.06, 0.06];
    let n_conditions = wind_speeds.len();

    fmodel.set_wind_conditions(
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(turbulence_intensities),
    )?;

    println!("\n========== Uniform Frequency Case ==========");
    println!("Number of conditions: {}", n_conditions);

    // Run the model
    fmodel.run()?;

    // With uniform frequency, each condition has equal weight
    let uniform_freq = vec![1.0 / n_conditions as f64; n_conditions];
    
    let expected_farm_power = fmodel.get_expected_farm_power(
        Some(Array1::from_vec(uniform_freq.clone())),
        None,
    )?;
    
    // Calculate AEP: expected power (kW) * hours per year
    let hours_per_year = 365.0 * 24.0;
    let aep = expected_farm_power * hours_per_year;

    println!("Expected farm power: {:.2} kW", expected_farm_power / 1000.0);
    println!("AEP: {:.2} MWh/year", aep / 1e6);
    println!("AEP: {:.3} GWh/year", aep / 1e9);

    //////////////////////////////////////////////////
    // Using WindRose-like conditions with frequencies
    //////////////////////////////////////////////////

    // Create a simple wind rose with 2 directions and 3 speeds
    let wd_values = vec![270.0, 280.0];
    let ws_values = vec![8.0, 9.0, 10.0];
    
    let mut rose_wind_directions = Vec::new();
    let mut rose_wind_speeds = Vec::new();
    let mut rose_tis = Vec::new();
    let mut freq_table = Vec::new();
    
    // Create frequency table (uniform for this example)
    let n_wd = wd_values.len();
    let n_ws = ws_values.len();
    let total_bins = n_wd * n_ws;
    let uniform_bin_freq = 1.0 / total_bins as f64;
    
    for &wd in &wd_values {
        for &ws in &ws_values {
            rose_wind_directions.push(wd);
            rose_wind_speeds.push(ws);
            rose_tis.push(0.06);
            freq_table.push(uniform_bin_freq);
        }
    }

    println!("Debug: Arrays created - speeds={}, directions={}, tis={}, freq={}",
             rose_wind_speeds.len(), rose_wind_directions.len(), 
             rose_tis.len(), freq_table.len());

    fmodel.set_wind_conditions(
        Array1::from_vec(rose_wind_speeds),
        Array1::from_vec(rose_wind_directions),
        Array1::from_vec(rose_tis),
    )?;

    println!("\n========== WindRose-like Case ==========");
    println!("Wind directions: {:?}", wd_values);
    println!("Wind speeds: {:?}", ws_values);
    println!("Total bins: {} × {} = {}", n_wd, n_ws, total_bins);
    println!("Frequency per bin: {:.4}", uniform_bin_freq);

    // Run the model
    match fmodel.run() {
        Ok(_) => println!("Model run completed for {} conditions", fmodel.n_findex()),
        Err(e) => {
            eprintln!("Error running model: {}", e);
            return Err(e);
        }
    }

    // Get the expected farm power using the frequency table
    let expected_farm_power = fmodel.get_expected_farm_power(
        Some(Array1::from_vec(freq_table.clone())),
        None,
    )?;
    
    // Calculate AEP
    let hours_per_year = 365.0 * 24.0;
    let aep = expected_farm_power * hours_per_year;

    println!("\nAEP from wind rose: {:.3} GWh/year", aep / 1e9);

    // Run the model again without wakes to compute wake losses
    fmodel.run_no_wake()?;

    // Get the expected farm power without wake
    let expected_farm_power_no_wake = fmodel.get_expected_farm_power(
        Some(Array1::from_vec(freq_table)),
        None,
    )?;
    
    let aep_no_wake = expected_farm_power_no_wake * hours_per_year;

    // Compute the wake losses
    let wake_losses = 100.0 * (aep_no_wake - aep) / aep_no_wake;

    println!("AEP without wakes: {:.3} GWh/year", aep_no_wake / 1e9);
    println!("Wake losses: {:.2}%", wake_losses);

    println!("\nExample 6 completed successfully!");

    Ok(())
}
