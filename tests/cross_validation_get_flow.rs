//! Cross-validation of the flow-sampling APIs against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/examples_get_flow/`:
//! - `002_extract_wind_speed_at_points.py`
//! - `003_extract_turbulence_intensity_at_points.py`
//! - `004_plot_velocity_deficit_profiles.py`

use florus::floris_model::FlorisModel;
use florus::types::Array1;

/// A two-turbine farm with the met mast in between the turbines
fn met_mast_model() -> FlorisModel {
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    let d = 126.0;
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 3.0 * d]),
            &Array1::from_vec(vec![0.0, 3.0 * d]),
        )
        .unwrap();

    // The wind direction is swept over 360 degrees
    let n = 360;
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; n]),
            Array1::from_iter((0..n).map(|i| i as f64)),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();
    fmodel
}

#[test]
fn test_sample_flow_at_points_matches_python() {
    // Python: examples_get_flow/002_extract_wind_speed_at_points.py
    let fmodel = met_mast_model();

    let d = 126.0;
    let points_x = vec![3.0 * d; 4];
    let points_y = vec![0.0; 4];
    let points_z = [30.0, 90.0, 150.0, 250.0];

    let u_at_points = fmodel
        .sample_flow_at_points(&points_x, &points_y, &points_z)
        .unwrap();

    // Python reference (m/s) at a few wind directions
    let cases = [
        (0usize, [5.398580, 3.074048, 6.548713, 9.038363]),
        (1, [5.426605, 3.159507, 6.582709, 9.038452]),
        (2, [5.490698, 3.354895, 6.660456, 9.038656]),
        (90, [7.011896, 8.000000, 8.505735, 9.043441]),
        (180, [7.011896, 8.000000, 8.505735, 9.043441]),
        (270, [5.398580, 3.074048, 6.548713, 9.038363]),
    ];
    for (wd, expected) in cases {
        for (k, value) in expected.iter().enumerate() {
            assert!(
                (u_at_points[[wd, k]] - value).abs() < 1e-5,
                "wind direction {wd}, point {k}: {} m/s (expected {value})",
                u_at_points[[wd, k]]
            );
        }
    }
}

#[test]
fn test_sample_ti_at_points_matches_python() {
    // Python: examples_get_flow/003_extract_turbulence_intensity_at_points.py
    let fmodel = met_mast_model();

    let d = 126.0;
    let points_x = vec![3.0 * d; 4];
    let points_y = vec![0.0; 4];
    let points_z = [30.0, 90.0, 150.0, 250.0];

    let ti_at_points = fmodel
        .sample_ti_at_points(&points_x, &points_y, &points_z)
        .unwrap();

    let cases = [
        (0usize, [0.110623, 0.110623, 0.110623, 0.060000]),
        (1, [0.110627, 0.110627, 0.110627, 0.060000]),
        (2, [0.110638, 0.110638, 0.110638, 0.060000]),
        (90, [0.060000, 0.060000, 0.060000, 0.060000]),
        (270, [0.110623, 0.110623, 0.110623, 0.060000]),
    ];
    for (wd, expected) in cases {
        for (k, value) in expected.iter().enumerate() {
            assert!(
                (ti_at_points[[wd, k]] - value).abs() < 1e-6,
                "wind direction {wd}, point {k}: {} (expected {value})",
                ti_at_points[[wd, k]]
            );
        }
    }
}

#[test]
fn test_velocity_deficit_profiles_match_python() {
    // Python: examples_get_flow/004_plot_velocity_deficit_profiles.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))
        .unwrap();

    let d = 125.88;
    let downstream_dists = [3.0 * d, 5.0 * d, 7.0 * d];

    let profiles = fmodel
        .sample_velocity_deficit_profiles(
            "cross-stream",
            Some(&downstream_dists),
            None,
            100,
            None,
            Some(8.0),
            None,
            0.0,
            0.0,
            Some(90.0),
        )
        .unwrap();

    // Python reference: the deficit of every tenth point of the first profile
    let expected = [
        (-2.000, 0.000000),
        (-1.596, 0.000012),
        (-1.192, 0.001477),
        (-0.788, 0.045366),
        (-0.384, 0.338477),
        (0.020, 0.613539),
        (0.424, 0.270190),
        (0.828, 0.028908),
        (1.232, 0.000751),
        (1.636, 0.000005),
    ];
    let profile = &profiles[0];
    assert!(
        (profile.x1[0] - 3.0).abs() < 1e-9,
        "x1/D = {}",
        profile.x1[0]
    );
    for (j, (x2, deficit)) in expected.iter().enumerate() {
        let index = j * 10;
        assert!(
            (profile.x2[index] - x2).abs() < 1e-3,
            "x2/D at point {index}: {}",
            profile.x2[index]
        );
        assert!(
            (profile.velocity_deficit[index] - deficit).abs() < 1e-6,
            "deficit at point {index}: {}",
            profile.velocity_deficit[index]
        );
    }
}
