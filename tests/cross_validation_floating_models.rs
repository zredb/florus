//! Cross-validation of the floating-turbine examples against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/examples_floating/`:
//! - `001_floating_turbine_models.py` (the fixed, floating and defined
//!   floating turbine models)
//! - `003_tilt_driven_vertical_wake_deflection.py` (5 and 15 degrees of tilt)

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

#[test]
fn test_floating_turbine_models_match_python() {
    // Python: examples_floating/001_floating_turbine_models.py
    let wind_speeds: Vec<f64> = (3..25).map(|i| i as f64).collect();
    let n = wind_speeds.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();

    let mut models = vec![
        (
            "fixed",
            FlorisModel::from_file("examples/inputs/gch_fixed.yaml").unwrap(),
        ),
        (
            "floating",
            FlorisModel::from_file("examples/inputs/gch_floating.yaml").unwrap(),
        ),
        (
            "defined_floating",
            FlorisModel::from_file("examples/inputs/gch_floating_defined_floating.yaml").unwrap(),
        ),
    ];

    let mut powers: Vec<Vec<f64>> = Vec::new();
    let mut cts: Vec<Vec<f64>> = Vec::new();
    for (_, fmodel) in models.iter_mut() {
        fmodel.set_wind_data(&time_series).unwrap();
        fmodel.run().unwrap();
        let turbine_powers = fmodel.get_turbine_powers().unwrap();
        powers.push(
            (0..n)
                .map(|i| turbine_powers[[i, 0]] / 1000.0)
                .collect::<Vec<f64>>(),
        );
        let turbine_cts = fmodel.get_turbine_thrust_coefficients();
        cts.push((0..n).map(|i| turbine_cts[[i, 0]]).collect::<Vec<f64>>());
    }

    // Python reference (kW and the thrust coefficient)
    let cases = [
        (0usize, 36.520, 36.520, 36.520, 0.999900, 0.999900, 0.999900),
        (1, 175.867, 175.867, 175.867, 0.999900, 0.999900, 0.999900),
        (2, 400.180, 399.540, 400.180, 0.919042, 0.918209, 0.919042),
        (10, 5000.0, 5000.0, 5000.0, 0.405452, 0.402593, 0.405452),
        (21, 5000.0, 5000.0, 5000.0, 0.065005, 0.064973, 0.065005),
    ];
    for (i, p_fixed, p_float, p_defined, ct_fixed, ct_float, ct_defined) in cases {
        assert!(
            (powers[0][i] - p_fixed).abs() < 0.01,
            "fixed power at {} m/s",
            wind_speeds[i]
        );
        assert!(
            (powers[1][i] - p_float).abs() < 0.01,
            "floating power at {} m/s",
            wind_speeds[i]
        );
        assert!(
            (powers[2][i] - p_defined).abs() < 0.01,
            "defined floating power at {} m/s",
            wind_speeds[i]
        );
        assert!(
            (cts[0][i] - ct_fixed).abs() < 1e-6,
            "fixed Ct at {} m/s",
            wind_speeds[i]
        );
        assert!(
            (cts[1][i] - ct_float).abs() < 1e-6,
            "floating Ct at {} m/s",
            wind_speeds[i]
        );
        assert!(
            (cts[2][i] - ct_defined).abs() < 1e-6,
            "defined floating Ct at {} m/s",
            wind_speeds[i]
        );
    }
}

#[test]
fn test_tilt_driven_vertical_wake_deflection_matches_python() {
    // Python: examples_floating/003_tilt_driven_vertical_wake_deflection.py
    let num_in_row = 5;
    for (input_file, tilt, expected) in [
        (
            "examples/inputs/emgauss_floating_fixedtilt5.yaml",
            5.0,
            [1754.0, 678.5, 704.2, 722.2, 733.3],
        ),
        (
            "examples/inputs/emgauss_floating_fixedtilt15.yaml",
            15.0,
            [1655.5, 747.7, 778.7, 784.5, 783.9],
        ),
    ] {
        let mut fmodel = FlorisModel::from_file(input_file).unwrap();
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
        fmodel.run().unwrap();

        let powers = fmodel.get_turbine_powers().unwrap();
        for t in 0..num_in_row {
            assert!(
                (powers[[0, t]] / 1000.0 - expected[t]).abs() < 0.1,
                "tilt {tilt} deg, turbine {t}: {} kW (expected {})",
                powers[[0, t]] / 1000.0,
                expected[t]
            );
        }
    }
}
