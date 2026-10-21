//! Cross-validation of the turbine-library examples against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/examples_turbine/`:
//! - `001_reference_turbines.py` (the turbines of the library)
//! - `002_multiple_turbine_types.py` (a farm with two turbine types)

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

#[test]
fn test_reference_turbines_match_python() -> anyhow::Result<()> {
    // Python: examples_turbine/001_reference_turbines.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))
        .unwrap();

    // The wind speeds of the comparison
    let ws_array: Vec<f64> = (1..150).map(|i| 0.1 + 0.2 * i as f64).collect();
    let n = ws_array.len();

    // Python reference: (turbine, wind speed, power [kW], thrust coefficient)
    let cases = [
        ("nrel_5MW", 8.1, 1825.992, 0.787033),
        ("nrel_5MW", 11.1, 4639.270, 0.753053),
        ("iea_10MW", 8.1, 4556.441, 0.870173),
        ("iea_10MW", 11.1, 10000.000, 0.581471),
    ];

    for (turbine, ws, power, ct) in cases {
        fmodel.set_turbine_types(&[turbine.to_string()]).unwrap();
        fmodel.reset_operation();
        // The reference wind height follows the hub height of the turbine
        fmodel.assign_hub_height_to_ref_height().unwrap();
        fmodel
            .set_wind_conditions(
                Array1::from_vec(ws_array.clone()),
                Array1::from_vec(vec![270.0; n]),
                Array1::from_vec(vec![0.06; n]),
            )
            .unwrap();
        fmodel.run().unwrap();

        let powers = fmodel.get_turbine_powers().unwrap();
        let cts = fmodel.get_turbine_thrust_coefficients();
        let index = ws_array
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (*a - ws).abs().partial_cmp(&(*b - ws).abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap();

        assert!(
            (powers[[index, 0]] / 1000.0 - power).abs() < 0.01,
            "{turbine} at {ws} m/s: {} kW (expected {power})",
            powers[[index, 0]] / 1000.0
        );
        assert!(
            (cts[[index, 0]] - ct).abs() < 1e-6,
            "{turbine} at {ws} m/s: Ct {} (expected {ct})",
            cts[[index, 0]]
        );
    }

    Ok(())
}

#[test]
fn test_yaw_effect_matches_python() -> anyhow::Result<()> {
    // Python: examples_turbine/001_reference_turbines.py (the yaw sweep)
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))
        .unwrap();
    fmodel.set_turbine_types(&["nrel_5MW".to_string()]).unwrap();
    fmodel.reset_operation();
    fmodel.assign_hub_height_to_ref_height().unwrap();

    let n_yaw = 60;
    let yaw_angles: Vec<f64> = (0..n_yaw)
        .map(|i| -30.0 + 60.0 * i as f64 / (n_yaw - 1) as f64)
        .collect();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![11.0; n_yaw]),
            Array1::from_vec(vec![270.0; n_yaw]),
            Array1::from_vec(vec![0.06; n_yaw]),
        )
        .unwrap();
    fmodel
        .set_yaw_angles(Array2::from_shape_vec(
            (n_yaw, 1),
            yaw_angles.iter().copied().collect::<Vec<f64>>(),
        )?)
        .unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();

    // Python reference (kW)
    let cases = [
        (0usize, -30.0, 3467.856),
        (30, 0.508, 4519.082),
        (59, 30.0, 3467.856),
    ];
    for (i, yaw, power) in cases {
        assert!(
            (powers[[i, 0]] / 1000.0 - power).abs() < 0.01,
            "yaw {yaw} deg: {} kW (expected {power})",
            powers[[i, 0]] / 1000.0
        );
    }

    Ok(())
}

#[test]
fn test_multiple_turbine_types_match_python() -> anyhow::Result<()> {
    // Python: examples_turbine/002_multiple_turbine_types.py
    let mut fmodel =
        FlorisModel::from_file("examples/inputs/gch_multiple_turbine_types.yaml").unwrap();
    fmodel.run().unwrap();

    // Python reference (kW)
    let powers = fmodel.get_turbine_powers().unwrap();
    assert!((powers[[0, 0]] / 1000.0 - 1754.0).abs() < 0.1);
    assert!((powers[[0, 1]] / 1000.0 - 2548.6).abs() < 0.1);

    Ok(())
}
