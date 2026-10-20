//! Cross-validation of the multi-dimensional Cp/Ct tables against the Python
//! FLORIS reference values.
//!
//! Reference values come from `floris/examples/examples_multidim/`:
//! - `001_multi_dimensional_cp_ct.py` (wave period / wave height conditions)
//! - `002_multi_dimensional_cp_ct_2Hs.py` (two wave heights)
//! - `003_multi_dimensional_cp_ct_TI.py` (turbulence intensity condition)

use florus::floris_config::MultiDimConditionValue;
use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::MultidimConditions;
use florus::types::{Array1, Array2};
use florus::wind_data::TimeSeries;
use std::collections::HashMap;

#[test]
fn test_multidim_single_condition_matches_python() {
    // Python: examples_multidim/001_multi_dimensional_cp_ct.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .unwrap();
    fmodel.set_yaw_angles(Array2::zeros((1, 2))).unwrap();
    fmodel.run().unwrap();

    let powers = fmodel.get_turbine_powers().unwrap() / 1000.0;
    // Python reference (kW)
    assert!(
        (powers[[0, 0]] - 3115.68688793).abs() < 0.1,
        "turbine 0: {}",
        powers[[0, 0]]
    );
    assert!(
        (powers[[0, 1]] - 1234.25106658).abs() < 0.1,
        "turbine 1: {}",
        powers[[0, 1]]
    );
}

#[test]
fn test_multidim_multiple_conditions_match_python() {
    // Python: examples_multidim/001_multi_dimensional_cp_ct.py (9 conditions)
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0, 9.0, 10.0, 8.0, 9.0, 10.0, 8.0, 9.0, 10.0]),
            Array1::from_vec(vec![
                260.0, 260.0, 260.0, 270.0, 270.0, 270.0, 280.0, 280.0, 280.0,
            ]),
            Array1::from_vec(vec![0.06; 9]),
        )
        .unwrap();
    fmodel.set_yaw_angles(Array2::zeros((9, 2))).unwrap();
    fmodel.run().unwrap();

    let powers = fmodel.get_turbine_powers().unwrap() / 1000.0;
    // Python reference (kW), turbine 0
    let expected = [
        3115.68688793,
        4415.57562618,
        6007.59995332,
        3115.68688793,
        4415.57562618,
        6007.59995332,
        3115.68688793,
        4415.57562618,
        6007.59995332,
    ];
    for fi in 0..9 {
        assert!(
            (powers[[fi, 0]] - expected[fi]).abs() < 0.1,
            "condition {fi}: {} kW (expected {})",
            powers[[fi, 0]],
            expected[fi]
        );
    }
}

#[test]
fn test_multidim_two_wave_heights_match_python() {
    // Python: examples_multidim/002_multi_dimensional_cp_ct_2Hs.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();

    let n_wind_speeds = 16;
    let mut wind_speeds: Vec<f64> = (0..n_wind_speeds)
        .map(|i| 5.0 + (20.0 - 5.0) * i as f64 / (n_wind_speeds - 1) as f64)
        .collect();
    let first_sweep = wind_speeds.clone();
    wind_speeds.extend_from_slice(&first_sweep);
    let n_conditions = wind_speeds.len();

    let mut hs: Vec<f64> = vec![3.1; n_wind_speeds];
    hs.extend(std::iter::repeat(1.0).take(n_wind_speeds));

    let mut time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_conditions]),
        Array1::from_vec(wind_speeds),
        Array1::from_vec(vec![0.06; n_conditions]),
    )
    .unwrap();
    time_series.multidim_conditions = Some(MultidimConditions {
        tp: Array1::from_vec(vec![2.5; n_conditions]),
        hs: Some(Array1::from_vec(hs)),
        ti: None,
    });
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();

    let powers = fmodel.get_turbine_powers().unwrap() / 1000.0;
    // Python reference (kW): (turbine, condition, power)
    let expected = [
        (0usize, 0usize, 687.0),
        (0, 15, 7500.0),
        (0, 16, 1374.1),
        (0, 31, 15000.0),
        (1, 0, 184.5),
        (1, 16, 0.0),
        (2, 0, 137.3),
        (2, 16, 49.0),
    ];
    for (t, fi, value) in expected {
        assert!(
            (powers[[fi, t]] - value).abs() < 0.1,
            "condition {fi}, turbine {t}: {} kW (expected {value})",
            powers[[fi, t]]
        );
    }
}

#[test]
fn test_multidim_turbulence_intensity_matches_python() {
    // Python: examples_multidim/003_multi_dimensional_cp_ct_TI.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch_multi_dim_cp_ct_TI.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();

    let wind_speeds: Vec<f64> = (0..150).map(|i| 5.0 + 0.1 * i as f64).collect();
    let n_conditions = wind_speeds.len();
    let mut time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_conditions]),
        Array1::from_vec(wind_speeds),
        Array1::from_vec(vec![0.06; n_conditions]),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();

    // TI = 0.06
    let mut conditions: HashMap<String, MultiDimConditionValue> = HashMap::new();
    conditions.insert("TI".to_string(), MultiDimConditionValue::Scalar(0.06));
    fmodel.set_multidim_conditions(conditions).unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap() / 1000.0;
    // Python reference (kW)
    assert!((powers[[0, 0]] - 1375.2).abs() < 0.1, "{}", powers[[0, 0]]);
    assert!((powers[[0, 1]] - 0.0).abs() < 0.1, "{}", powers[[0, 1]]);
    assert!((powers[[0, 2]] - 43.9).abs() < 0.1, "{}", powers[[0, 2]]);

    // TI = 0.10 uses the nearest defined value (0.08)
    let mut conditions: HashMap<String, MultiDimConditionValue> = HashMap::new();
    conditions.insert("TI".to_string(), MultiDimConditionValue::Scalar(0.10));
    fmodel.set_multidim_conditions(conditions).unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap() / 1000.0;
    assert!((powers[[0, 0]] - 1512.8).abs() < 0.1, "{}", powers[[0, 0]]);
    assert!((powers[[0, 2]] - 48.3).abs() < 0.1, "{}", powers[[0, 2]]);
}
