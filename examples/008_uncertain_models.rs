/// Example 8: Uncertain Models (Simplified)
///
/// This example demonstrates the concept of uncertainty modeling in FLORIS.
/// Note: Full UncertainFlorisModel implementation is not yet available in Rust.
/// This example shows how to manually simulate wind direction uncertainty by
/// running multiple cases with perturbed wind directions.
///
/// In Python FLORIS, UncertainFlorisModel adds uncertainty to the inflow wind
/// direction by expanding the wind direction time series to include uncertainty
/// and then computing a Gaussian-weighted average of the results.

use florus::{Array1, FlorisModel};
use florus::core::turbines::TurbineLibrary;
use std::f64::consts::PI;

fn main() -> florus::Result<()> {
    TurbineLibrary::init_if_needed()?;

    // Instantiate FLORIS model
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // Define a two turbine farm
    let d = 126.0;
    let layout_x = vec![0.0, 6.0 * d];
    let layout_y = vec![0.0, 0.0];
    
    fmodel.set_layout(
        &Array1::from_vec(layout_x.clone()),
        &Array1::from_vec(layout_y.clone()),
    )?;

    println!("Layout: 2 turbines at [0, {:.0}]m (6D spacing)", 6.0 * d);
    println!("  n_turbines: {}", fmodel.n_turbines());

    // Define wind directions to sweep
    let wind_directions_base: Vec<f64> = (240..300).map(|x| x as f64).collect();
    let n_wd = wind_directions_base.len();
    let wind_speeds = vec![8.0; n_wd];
    let turbulence_intensities = vec![0.06; n_wd];

    //////////////////////////////////////////////////
    // Run nominal (no uncertainty) case
    //////////////////////////////////////////////////

    println!("\n========== Running Nominal Case ==========");
    fmodel.set_wind_conditions(
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(wind_directions_base.clone()),
        Array1::from_vec(turbulence_intensities.clone()),
    )?;

    fmodel.run()?;
    let turbine_powers_nom = fmodel.get_turbine_powers();
    let farm_powers_nom = fmodel.get_farm_power();

    println!("Nominal case completed.");
    println!("Sample results at WD=270°:");
    let idx_270 = 30; // 270 - 240 = 30
    if idx_270 < n_wd {
        println!(
            "  T1={:.0} kW, T2={:.0} kW, Farm={:.0} kW",
            turbine_powers_nom[[idx_270, 0]] / 1000.0,
            turbine_powers_nom[[idx_270, 1]] / 1000.0,
            farm_powers_nom[idx_270] / 1000.0
        );
    }

    //////////////////////////////////////////////////
    // Simulate uncertainty manually
    //////////////////////////////////////////////////

    println!("\n========== Simulating Wind Direction Uncertainty ==========");
    
    // For demonstration, we'll manually compute uncertain power at one wind direction
    // by averaging over a range of perturbed directions
    
    let wd_std = 3.0; // Standard deviation in degrees
    let wd_center = 270.0;
    let n_samples = 21; // Number of samples for uncertainty integration
    let wd_range = 4.0 * wd_std; // ±4 sigma covers ~99.99% of distribution
    
    let mut wd_samples = Vec::new();
    let mut weights = Vec::new();
    let mut total_weight = 0.0;
    
    // Generate samples with Gaussian weighting
    for i in 0..n_samples {
        let wd = wd_center - wd_range + (2.0 * wd_range * i as f64 / (n_samples - 1) as f64);
        let weight = gaussian_pdf(wd, wd_center, wd_std);
        wd_samples.push(wd);
        weights.push(weight);
        total_weight += weight;
    }
    
    // Normalize weights
    for w in weights.iter_mut() {
        *w /= total_weight;
    }
    
    // Create a new model for uncertainty simulation to avoid state issues
    let mut fmodel_unc = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    fmodel_unc.set_layout(
        &Array1::from_vec(layout_x),
        &Array1::from_vec(layout_y),
    )?;
    
    // Set wind conditions for all samples
    let ws_samples = vec![8.0; n_samples];
    let ti_samples = vec![0.06; n_samples];
    
    fmodel_unc.set_wind_conditions(
        Array1::from_vec(ws_samples),
        Array1::from_vec(wd_samples.clone()),
        Array1::from_vec(ti_samples),
    )?;
    
    fmodel_unc.run()?;
    let turbine_powers_unc = fmodel_unc.get_turbine_powers();
    
    // Compute weighted average
    let mut t1_power_unc = 0.0;
    let mut t2_power_unc = 0.0;
    let mut farm_power_unc = 0.0;
    
    for i in 0..n_samples {
        t1_power_unc += turbine_powers_unc[[i, 0]] * weights[i];
        t2_power_unc += turbine_powers_unc[[i, 1]] * weights[i];
        farm_power_unc += (turbine_powers_unc[[i, 0]] + turbine_powers_unc[[i, 1]]) * weights[i];
    }
    
    println!("Uncertainty simulation (wd_std = {:.1}°):", wd_std);
    println!("  Center WD: {:.0}°", wd_center);
    println!("  Samples: {} (range: {:.1}° to {:.1}°)", 
             n_samples, wd_samples[0], wd_samples[n_samples-1]);
    println!("\n  Weighted average powers:");
    println!("    T1={:.0} kW, T2={:.0} kW, Farm={:.0} kW",
             t1_power_unc / 1000.0, t2_power_unc / 1000.0, farm_power_unc / 1000.0);
    
    // Compare with nominal at center
    let t1_nom = turbine_powers_nom[[idx_270, 0]];
    let t2_nom = turbine_powers_nom[[idx_270, 1]];
    let farm_nom = farm_powers_nom[idx_270];
    
    println!("\n  Nominal powers at WD=270°:");
    println!("    T1={:.0} kW, T2={:.0} kW, Farm={:.0} kW",
             t1_nom / 1000.0, t2_nom / 1000.0, farm_nom / 1000.0);
    
    println!("\n  Difference due to uncertainty:");
    println!("    T1: {:+.1}%", 100.0 * (t1_power_unc - t1_nom) / t1_nom);
    println!("    T2: {:+.1}%", 100.0 * (t2_power_unc - t2_nom) / t2_nom);
    println!("    Farm: {:+.1}%", 100.0 * (farm_power_unc - farm_nom) / farm_nom);

    //////////////////////////////////////////////////
    // Calculate AEP with uniform frequencies
    //////////////////////////////////////////////////

    println!("\n========== AEP Calculation ==========");
    
    let aep_nom = fmodel.get_farm_aep_uniform(8760.0);
    println!("AEP without uncertainty: {:.2} GWh/year", aep_nom / 1e6);
    
    println!("\nExample 8 completed successfully!");
    println!("Note: Full UncertainFlorisModel will be implemented in future versions.");

    Ok(())
}

/// Gaussian probability density function
fn gaussian_pdf(x: f64, mean: f64, std: f64) -> f64 {
    let coeff = 1.0 / (std * (2.0 * PI).sqrt());
    let exponent = -0.5 * ((x - mean) / std).powi(2);
    coeff * exponent.exp()
}
