//! Example: The empirical Gauss velocity deficit parameters
//!
//! This example shows the effect of the velocity deficit parameters of the
//! empirical Gauss model on the wakes and the powers of a five-turbine row:
//! the wake expansion rates, the breakpoints and the mixing gain are varied.
//!
//! Scenario aligned 1:1 with
//! examples_emgauss/001_empirical_gauss_velocity_deficit_parameters.py.

use florus::floris_config::VelocityModelConfig;
use florus::floris_model::FlorisModel;
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    let num_in_row = 5;

    // The original model
    let mut fmodel = FlorisModel::from_file("examples/inputs/emgauss.yaml")?;
    let original_config = fmodel.wake_config().clone();
    let d = fmodel.core().farm.rotor_diameters()[0];
    fmodel.set_layout(
        &Array1::from_iter((0..num_in_row).map(|x| x as f64 * 5.0 * d)),
        &Array1::from_vec(vec![0.0; num_in_row]),
    )?;
    fmodel.set_wind_conditions(
        Array1::from_vec(vec![8.0]),
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![0.06]),
    )?;

    // The powers of the original model
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let original: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Original: {}", original.join(" "));

    // Increase the base recovery rate
    let mut config = original_config.clone();
    if let VelocityModelConfig::EmpiricalGauss {
        wake_expansion_rates,
        ..
    } = &mut config.model_strings.velocity_model
    {
        *wake_expansion_rates = vec![0.03, 0.015];
    }
    fmodel.set_wake_config(config)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let increased: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Increase base recovery: {}", increased.join(" "));

    // Add a rate and change the breakpoints
    let mut config = original_config.clone();
    if let VelocityModelConfig::EmpiricalGauss {
        wake_expansion_rates,
        breakpoints_d,
        ..
    } = &mut config.model_strings.velocity_model
    {
        wake_expansion_rates.push(0.0);
        *breakpoints_d = vec![5.0, 10.0];
    }
    fmodel.set_wake_config(config)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let added: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Add rate, change breakpoints: {}", added.join(" "));

    // Increase the wake-induced mixing gain
    let mut config = original_config.clone();
    if let VelocityModelConfig::EmpiricalGauss {
        mixing_gain_velocity,
        ..
    } = &mut config.model_strings.velocity_model
    {
        *mixing_gain_velocity = 3.0;
    }
    fmodel.set_wake_config(config)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let mixing: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Increase mixing gain: {}", mixing.join(" "));

    Ok(())
}
