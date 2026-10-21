//! Cross-validation of the empirical Gauss parameter examples against the
//! Python FLORIS reference values.
//!
//! Reference values come from `floris/examples/examples_emgauss/`:
//! - `001_empirical_gauss_velocity_deficit_parameters.py`
//! - `002_empirical_gauss_deflection_parameters.py`

use florus::floris_config::{DeflectionModelConfig, VelocityModelConfig};
use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

/// The turbine powers [MW] of the five-turbine row
fn powers(fmodel: &FlorisModel) -> Vec<f64> {
    let powers = fmodel.get_turbine_powers().unwrap();
    (0..powers.shape()[1])
        .map(|t| powers[[0, t]] / 1.0e6)
        .collect()
}

#[test]
fn test_velocity_deficit_parameters_match_python() {
    // Python: examples_emgauss/001_empirical_gauss_velocity_deficit_parameters.py
    let num_in_row = 5;
    let mut fmodel = FlorisModel::from_file("examples/inputs/emgauss.yaml").unwrap();
    let original_config = fmodel.wake_config().clone();
    let d = fmodel.core().farm.rotor_diameters()[0];
    fmodel
        .set_layout(
            &Array1::from_iter((0..num_in_row).map(|x| x as f64 * 5.0 * d)),
            &Array1::from_vec(vec![0.0; num_in_row]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .unwrap();

    // The original model
    fmodel.run().unwrap();
    let original = powers(&fmodel);
    let expected = [1.754, 0.678, 0.704, 0.722, 0.733];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (original[t] - value).abs() < 0.001,
            "original power of turbine {t}: {} MW (expected {value})",
            original[t]
        );
    }

    // Increase the base recovery rate
    let mut config = original_config.clone();
    if let VelocityModelConfig::EmpiricalGauss {
        wake_expansion_rates,
        ..
    } = &mut config.model_strings.velocity_model
    {
        *wake_expansion_rates = vec![0.03, 0.015];
    }
    fmodel.set_wake_config(config).unwrap();
    fmodel.run().unwrap();
    let increased = powers(&fmodel);
    let expected = [1.754, 0.791, 0.845, 0.874, 0.888];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (increased[t] - value).abs() < 0.001,
            "increased recovery power of turbine {t}: {} MW (expected {value})",
            increased[t]
        );
    }

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
    fmodel.set_wake_config(config).unwrap();
    fmodel.run().unwrap();
    let added = powers(&fmodel);
    let expected = [1.754, 0.671, 0.597, 0.559, 0.541];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (added[t] - value).abs() < 0.001,
            "added rate power of turbine {t}: {} MW (expected {value})",
            added[t]
        );
    }

    // Increase the wake-induced mixing gain
    let mut config = original_config.clone();
    if let VelocityModelConfig::EmpiricalGauss {
        mixing_gain_velocity,
        ..
    } = &mut config.model_strings.velocity_model
    {
        *mixing_gain_velocity = 3.0;
    }
    fmodel.set_wake_config(config).unwrap();
    fmodel.run().unwrap();
    let mixing = powers(&fmodel);
    let expected = [1.754, 0.678, 0.778, 0.835, 0.864];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (mixing[t] - value).abs() < 0.001,
            "mixing gain power of turbine {t}: {} MW (expected {value})",
            mixing[t]
        );
    }
}

#[test]
fn test_deflection_parameters_match_python() {
    // Python: examples_emgauss/002_empirical_gauss_deflection_parameters.py
    let num_in_row = 5;
    let yaw_angles = vec![20.0, 20.0, 10.0, 0.0, 0.0];

    let mut fmodel = FlorisModel::from_file("examples/inputs/emgauss.yaml").unwrap();
    let original_config = fmodel.wake_config().clone();
    let d = fmodel.core().farm.rotor_diameters()[0];
    fmodel
        .set_layout(
            &Array1::from_iter((0..num_in_row).map(|x| x as f64 * 5.0 * d)),
            &Array1::from_vec(vec![0.0; num_in_row]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .unwrap();
    fmodel
        .set_yaw_angles(Array2::from_shape_vec((1, num_in_row), yaw_angles.clone()).unwrap())
        .unwrap();

    // The original model
    fmodel.run().unwrap();
    let original = powers(&fmodel);
    let expected = [1.561, 0.829, 0.940, 0.872, 0.818];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (original[t] - value).abs() < 0.001,
            "original power of turbine {t}: {} MW (expected {value})",
            original[t]
        );
    }

    // Increase the maximum deflection attained
    let mut config = original_config.clone();
    if let DeflectionModelConfig::EmpiricalGauss {
        horizontal_deflection_gain_d,
        ..
    } = &mut config.model_strings.deflection_model
    {
        *horizontal_deflection_gain_d = 5.0;
    }
    fmodel.set_wake_config(config).unwrap();
    fmodel.run().unwrap();
    let deflection = powers(&fmodel);
    let expected = [1.561, 1.003, 1.117, 0.987, 0.879];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (deflection[t] - value).abs() < 0.001,
            "max deflection power of turbine {t}: {} MW (expected {value})",
            deflection[t]
        );
    }

    // Increase the influence of the wake added mixing
    let mut config = original_config.clone();
    if let DeflectionModelConfig::EmpiricalGauss {
        mixing_gain_deflection,
        ..
    } = &mut config.model_strings.deflection_model
    {
        *mixing_gain_deflection = 100.0;
    }
    fmodel.set_wake_config(config).unwrap();
    fmodel.run().unwrap();
    let mixing = powers(&fmodel);
    let expected = [1.561, 0.829, 0.891, 0.838, 0.805];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (mixing[t] - value).abs() < 0.001,
            "mixing gain power of turbine {t}: {} MW (expected {value})",
            mixing[t]
        );
    }

    // Increase the yaw-added mixing contribution
    let mut config = original_config.clone();
    if let DeflectionModelConfig::EmpiricalGauss {
        mixing_gain_deflection,
        yaw_added_mixing_gain,
        ..
    } = &mut config.model_strings.deflection_model
    {
        *mixing_gain_deflection = 100.0;
        *yaw_added_mixing_gain = 100.0;
    }
    fmodel.set_wake_config(config).unwrap();
    fmodel.run().unwrap();
    let yaw_added = powers(&fmodel);
    let expected = [1.561, 1.560, 1.703, 1.736, 0.678];
    for (t, value) in expected.iter().enumerate() {
        assert!(
            (yaw_added[t] - value).abs() < 0.001,
            "yaw-added mixing power of turbine {t}: {} MW (expected {value})",
            yaw_added[t]
        );
    }
}
