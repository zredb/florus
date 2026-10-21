//! Cross-validation of the operation-model example against the Python FLORIS
//! reference values.
//!
//! Reference values come from
//! `floris/examples/examples_operation_models/001_compare_yaw_loss.py`.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};

#[test]
fn test_operation_model_yaw_losses_match_python() -> anyhow::Result<()> {
    let n = 101;
    let yaw_max = 30.0;
    let yaw_angles: Vec<f64> = (0..n)
        .map(|i| -yaw_max + 2.0 * yaw_max * i as f64 / (n - 1) as f64)
        .collect();

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(&Array1::from_vec(vec![0.0]), &Array1::from_vec(vec![0.0]))
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![1.0; n]),
            Array1::from_vec(vec![270.0; n]),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();
    fmodel
        .set_yaw_angles(Array2::from_shape_vec(
            (n, 1),
            yaw_angles.iter().copied().collect::<Vec<f64>>(),
        )?)
        .unwrap();

    // Python reference (kW): the model, the wind speed and the powers at
    // -30 deg, 0 deg and +30 deg
    let cases = [
        ("simple", 8.0, 1753.954, 1753.954, 1753.954),
        ("simple", 11.5, 5000.0, 5000.0, 5000.0),
        ("cosine-loss", 8.0, 1339.856, 1753.954, 1339.856),
        ("cosine-loss", 11.5, 3955.922, 5000.0, 3955.922),
        ("unified-momentum", 8.0, 1371.941, 1753.954, 1371.941),
        ("unified-momentum", 11.5, 3908.773, 5000.0, 3908.773),
        // The controller-dependent model iterates, so the values agree to 0.01 kW
        ("controller-dependent", 8.0, 1349.717, 1728.652, 1356.395),
    ];

    for (model, wind_speed, p_minus, p_zero, p_plus) in cases {
        fmodel
            .set_wind_conditions(
                Array1::from_vec(vec![wind_speed; n]),
                Array1::from_vec(vec![270.0; n]),
                Array1::from_vec(vec![0.06; n]),
            )
            .unwrap();
        fmodel.set_operation_model(model).unwrap();
        fmodel.run().unwrap();
        let powers = fmodel.get_turbine_powers().unwrap();

        assert!(
            (powers[[0, 0]] / 1000.0 - p_minus).abs() < 0.01,
            "{model} at {wind_speed} m/s: {} kW (expected {p_minus})",
            powers[[0, 0]] / 1000.0
        );
        assert!(
            (powers[[n / 2, 0]] / 1000.0 - p_zero).abs() < 0.01,
            "{model} at {wind_speed} m/s, zero yaw: {} kW (expected {p_zero})",
            powers[[n / 2, 0]] / 1000.0
        );
        assert!(
            (powers[[n - 1, 0]] / 1000.0 - p_plus).abs() < 0.01,
            "{model} at {wind_speed} m/s: {} kW (expected {p_plus})",
            powers[[n - 1, 0]] / 1000.0
        );
    }

    Ok(())
}
