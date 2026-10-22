//! Example: The empirical Gauss deflection parameters
//!
//! This example shows the effect of the deflection parameters of the
//! empirical Gauss model on the wakes and the powers of a five-turbine row
//! with yawed upstream turbines: the maximum deflection, the mixing gain and
//! the yaw-added mixing gain are varied.
//!
//! Scenario aligned 1:1 with
//! examples_emgauss/002_empirical_gauss_deflection_parameters.py.

use florus::floris_config::DeflectionModelConfig;
use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

fn main() -> anyhow::Result<()> {
    let num_in_row = 5;

    // The yaw angles of the first three turbines
    let yaw_angles = vec![20.0, 20.0, 10.0, 0.0, 0.0];
    println!(
        "Turbine yaw angles (degrees):  [{}]",
        yaw_angles
            .iter()
            .map(|y| format!("{:.0}", y))
            .collect::<Vec<String>>()
            .join(" ")
    );

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
    fmodel.set_yaw_angles(Array2::from_shape_vec((1, num_in_row), yaw_angles.clone())?)?;

    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let original: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Original: {}", original.join(" "));

    // Increase the maximum deflection attained
    let mut config = original_config.clone();
    if let DeflectionModelConfig::EmpiricalGauss {
        horizontal_deflection_gain_d,
        ..
    } = &mut config.model_strings.deflection_model
    {
        *horizontal_deflection_gain_d = 5.0;
    }
    fmodel.set_wake_config(config)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let deflection: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Increase max deflection: {}", deflection.join(" "));

    // Increase the influence of the wake added mixing
    let mut config = original_config.clone();
    if let DeflectionModelConfig::EmpiricalGauss {
        mixing_gain_deflection,
        ..
    } = &mut config.model_strings.deflection_model
    {
        *mixing_gain_deflection = 100.0;
    }
    fmodel.set_wake_config(config)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let mixing: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Increase mixing gain: {}", mixing.join(" "));

    // Increase the yaw-added mixing contribution (the yaw-added mixing gain
    // is included so that it is reflected in the deflection as well)
    let mut config = original_config.clone();
    if let DeflectionModelConfig::EmpiricalGauss {
        mixing_gain_deflection,
        yaw_added_mixing_gain,
        ..
    } = &mut config.model_strings.deflection_model
    {
        *mixing_gain_deflection = 100.0;
        *yaw_added_mixing_gain = 1.0;
    }
    fmodel.set_wake_config(config)?;
    fmodel.run()?;
    let powers = fmodel.get_turbine_powers()?;
    let yaw_added: Vec<String> = (0..num_in_row)
        .map(|t| format!("{:.3}", powers[[0, t]] / 1.0e6))
        .collect();
    println!("Increase yaw-added mixing gain: {}", yaw_added.join(" "));

    Ok(())
}
