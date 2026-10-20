//! Cross-validation of the heterogeneous inflow against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/examples_heterogeneous/`:
//! - `001_heterogeneous_inflow_single.py` (single condition)
//! - `002_heterogeneous_using_wind_data.py` (per-condition multipliers)
//! - `003_heterogeneous_speedup_by_wd_and_ws.py` (selection by direction/speed)
//! - `004_heterogeneous_2d_and_3d.py` (YAML configuration, 2D and 3D)

use florus::floris_model::FlorisModel;
use florus::heterogeneous_map::{HeterogeneousInflowConfig, HeterogeneousMap};
use florus::types::{Array1, Array2};

/// The heterogeneous box used by the FLORIS examples
fn example_box() -> (Array1, Array1) {
    (
        Array1::from_vec(vec![-500.0, -500.0, 1000.0, 1000.0]),
        Array1::from_vec(vec![-500.0, 1000.0, -500.0, 1000.0]),
    )
}

#[test]
fn test_heterogeneous_single_condition_matches_python() {
    // Python: examples_heterogeneous/001_heterogeneous_inflow_single.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 0.0, 500.0, 500.0]),
            &Array1::from_vec(vec![0.0, 500.0, 0.0, 500.0]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .unwrap();

    let (x, y) = example_box();
    fmodel
        .set_heterogeneous_inflow_config(HeterogeneousInflowConfig {
            x,
            y,
            z: None,
            wind_speeds: None,
            wind_directions: None,
            speed_multipliers: Array2::from_shape_vec((1, 4), vec![1.0, 1.25, 1.0, 1.25]).unwrap(),
        })
        .unwrap();
    fmodel.run().unwrap();

    let powers = fmodel.get_turbine_powers().unwrap() / 1000.0;
    let expected = [2248.2, 2800.1, 466.2, 601.5];
    for t in 0..4 {
        assert!(
            (powers[[0, t]] - expected[t]).abs() < 0.2,
            "turbine {t}: got {} kW, expected {} kW",
            powers[[0, t]],
            expected[t]
        );
    }
}

#[test]
fn test_heterogeneous_per_condition_multipliers_match_python() {
    // Python: examples_heterogeneous/002_heterogeneous_using_wind_data.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 0.0, 500.0, 500.0]),
            &Array1::from_vec(vec![0.0, 500.0, 0.0, 500.0]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; 4]),
            Array1::from_vec(vec![269.0, 270.0, 271.0, 282.0]),
            Array1::from_vec(vec![0.06; 4]),
        )
        .unwrap();

    let (x, y) = example_box();
    fmodel
        .set_heterogeneous_inflow_config(HeterogeneousInflowConfig {
            x,
            y,
            z: None,
            wind_speeds: None,
            wind_directions: None,
            speed_multipliers: Array2::from_shape_vec(
                (4, 4),
                vec![
                    1.0, 1.25, 1.0, 1.25, //
                    1.0, 1.25, 1.0, 1.25, //
                    1.0, 1.25, 1.0, 1.25, //
                    1.0, 1.35, 1.0, 1.35,
                ],
            )
            .unwrap(),
        })
        .unwrap();
    fmodel.run().unwrap();

    let powers = fmodel.get_turbine_powers().unwrap() / 1000.0;
    // Python reference (turbine, findex) in kW
    let expected = [
        [2248.2, 2248.2, 2248.2, 2446.9],
        [2800.1, 2800.1, 2800.1, 3294.4],
        [487.8, 466.2, 492.3, 2130.0],
        [623.5, 601.5, 627.7, 2856.0],
    ];
    for t in 0..4 {
        for fi in 0..4 {
            assert!(
                (powers[[fi, t]] - expected[t][fi]).abs() < 0.2,
                "turbine {t}, findex {fi}: got {} kW, expected {} kW",
                powers[[fi, t]],
                expected[t][fi]
            );
        }
    }
}

#[test]
fn test_heterogeneous_map_selection_by_wd_and_ws_matches_python() {
    // Python: examples_heterogeneous/003_heterogeneous_speedup_by_wd_and_ws.py
    let heterogeneous_map = HeterogeneousMap::new(
        Array1::from_vec(vec![0.0, 0.0, 250.0, 500.0, 500.0]),
        Array1::from_vec(vec![0.0, 500.0, 250.0, 0.0, 500.0]),
        Array2::from_shape_vec(
            (4, 5),
            vec![
                1.0, 1.0, 1.0, 1.0, 1.0, //
                1.0, 1.0, 1.0, 1.0, 1.0, //
                1.5, 1.0, 1.25, 1.5, 1.0, //
                1.0, 1.5, 1.25, 1.0, 1.5,
            ],
        )
        .unwrap(),
        None,
        Some(Array1::from_vec(vec![270.0, 270.0, 90.0, 90.0])),
        Some(Array1::from_vec(vec![5.0, 10.0, 5.0, 10.0])),
        florus::core::models::InterpMethod::Linear,
    )
    .unwrap();

    let config = heterogeneous_map
        .get_heterogeneous_inflow_config(
            Array1::from_vec(vec![275.0, 95.0, 75.0]),
            Array1::from_vec(vec![7.0, 6.2, 8.0]),
        )
        .unwrap();

    let expected = [
        [1.0, 1.0, 1.0, 1.0, 1.0],
        [1.5, 1.0, 1.25, 1.5, 1.0],
        [1.0, 1.5, 1.25, 1.0, 1.5],
    ];
    for fi in 0..3 {
        for p in 0..5 {
            assert!(
                (config.speed_multipliers[[fi, p]] - expected[fi][p]).abs() < 1e-9,
                "condition {fi}, point {p}: got {}, expected {}",
                config.speed_multipliers[[fi, p]],
                expected[fi][p]
            );
        }
    }
}

#[test]
fn test_heterogeneous_from_input_file_matches_python() {
    // Python: examples_heterogeneous/004_heterogeneous_2d_and_3d.py (2D case).
    // The heterogeneous inflow is defined in the input file.
    let mut fmodel =
        FlorisModel::from_file("examples/inputs/gch_heterogeneous_inflow.yaml").unwrap();
    fmodel.set_wind_shear(0.0).unwrap();
    fmodel.run().unwrap();

    let plane = fmodel.calculate_horizontal_plane(90.0, 200, 100).unwrap();
    let u = &plane.data.u;
    let min = u.iter().copied().fold(f64::INFINITY, f64::min);
    let max = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = u.iter().sum::<f64>() / u.len() as f64;

    // Python reference: u range -0.098 .. 15.357 m/s, mean 10.666 m/s
    assert!((min - (-0.098)).abs() < 0.05, "min u {min}");
    assert!((max - 15.357).abs() < 0.05, "max u {max}");
    assert!((mean - 10.666).abs() < 0.05, "mean u {mean}");
}

#[test]
fn test_heterogeneous_3d_matches_python() {
    // Python: examples_heterogeneous/004_heterogeneous_2d_and_3d.py (3D case)
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel.set_wind_shear(0.0).unwrap();
    fmodel
        .set_heterogeneous_inflow_config(HeterogeneousInflowConfig {
            x: Array1::from_vec(vec![
                -300.0, -300.0, -300.0, -300.0, 2600.0, 2600.0, 2600.0, 2600.0,
            ]),
            y: Array1::from_vec(vec![
                -300.0, 300.0, -300.0, 300.0, -300.0, 300.0, -300.0, 300.0,
            ]),
            z: Some(Array1::from_vec(vec![
                540.0, 540.0, 0.0, 0.0, 540.0, 540.0, 0.0, 0.0,
            ])),
            wind_speeds: None,
            wind_directions: None,
            speed_multipliers: Array2::from_shape_vec(
                (1, 8),
                vec![1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0],
            )
            .unwrap(),
        })
        .unwrap();
    fmodel.run().unwrap();

    let plane = fmodel.calculate_horizontal_plane(90.0, 200, 100).unwrap();
    let u = &plane.data.u;
    let min = u.iter().copied().fold(f64::INFINITY, f64::min);
    let max = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = u.iter().sum::<f64>() / u.len() as f64;

    // Python reference: u range 0.457 .. 14.667 m/s, mean 13.431 m/s
    assert!((min - 0.457).abs() < 0.05, "min u {min}");
    assert!((max - 14.667).abs() < 0.05, "max u {max}");
    assert!((mean - 13.431).abs() < 0.05, "mean u {mean}");
}
