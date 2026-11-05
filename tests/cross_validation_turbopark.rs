//! Cross-validation of the TurbOPark example against the Python FLORIS
//! reference values.
//!
//! Reference values come from
//! `floris/examples/examples_turbopark/001_compare_turbopark_implementations.py`.

use florus::core::turbines::turbine_type::TurbineType;
use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

/// A constant thrust coefficient turbine (Python
/// `build_cosine_loss_turbine_dict`).
fn constant_ct_turbine() -> TurbineType {
    TurbineType::build(
        "ConstantCT",
        &[0.0, 30.0],
        Some(&[0.0, 1.0]),
        None,
        None,
        Some(&[0.75, 0.75]),
        None,
        100.0, // hub height
        120.0, // rotor diameter
        8.0,   // TSR
        1.225, // reference air density
        0.0,   // reference tilt
        1.88,  // cosine loss exponent yaw
        1.88,  // cosine loss exponent tilt
    )
    .expect("failed to build the turbine")
}

#[test]
fn test_turbopark_matches_python() {
    let turbine = constant_ct_turbine();

    // The original TurbOPark model
    let mut fmodel = FlorisModel::from_file("examples/inputs/turbopark_cubature.yaml").unwrap();
    fmodel.set_turbine_type(&turbine).unwrap();

    let wd_array: Vec<f64> = (0..900).map(|i| 225.0 + 0.1 * i as f64).collect();
    let n = wd_array.len();
    let wind_data = TimeSeries::new(
        Array1::from_vec(wd_array.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 600.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&wind_data).unwrap();
    fmodel.run().unwrap();
    let velocities = fmodel.turbine_average_velocities();

    // Python reference (the normalized waked speed at the downstream turbine)
    let cases = [(450usize, 0.779053), (451, 0.779082), (899, 1.000000)];
    for (index, expected) in cases {
        assert!(
            (velocities[[index, 1]] / 8.0 - expected).abs() < 1e-4,
            "wind direction {}: {} (expected {expected})",
            wd_array[index],
            velocities[[index, 1]] / 8.0
        );
    }
}

#[test]
fn test_turbopark_gauss_matches_python() {
    let turbine = constant_ct_turbine();

    // The TurbOPark-Gauss model
    let mut fmodel =
        FlorisModel::from_file("examples/inputs/turboparkgauss_cubature.yaml").unwrap();
    fmodel.set_turbine_type(&turbine).unwrap();

    let wd_array: Vec<f64> = (0..900).map(|i| 225.0 + 0.1 * i as f64).collect();
    let n = wd_array.len();
    let wind_data = TimeSeries::new(
        Array1::from_vec(wd_array.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 600.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&wind_data).unwrap();
    fmodel.run().unwrap();
    let velocities = fmodel.turbine_average_velocities();

    // Python reference
    let cases = [(450usize, 0.709920), (451, 0.709970), (899, 1.000000)];
    for (index, expected) in cases {
        assert!(
            (velocities[[index, 1]] / 8.0 - expected).abs() < 1e-4,
            "wind direction {}: {} (expected {expected})",
            wd_array[index],
            velocities[[index, 1]] / 8.0
        );
    }
}
