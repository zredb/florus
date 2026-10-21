//! Cross-validation of the load optimization examples against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/examples_load_optimization/`:
//! - `001_lti_and_voc.py` (the load turbulence intensity and the VOC)
//! - `002_row_opt_example.py` (the load optimization of a row of turbines)

use florus::core::turbines::{POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED};
use florus::floris_model::FlorisModel;
use florus::optimization::load_optimization::load_optimization::{
    compute_farm_revenue, compute_farm_voc, compute_lti, compute_net_revenue, compute_turbine_voc,
    optimize_power_setpoints,
};
use florus::types::{Array1, Array2};
use florus::utilities::linspace;
use florus::wind_data::TimeSeries;

#[test]
fn test_lti_and_voc_match_python() {
    // Python: examples_load_optimization/001_lti_and_voc.py
    let mut fmodel = FlorisModel::from_default().unwrap();
    fmodel.set_operation_model("simple-derating").unwrap();

    // The wind direction sweep
    let wind_directions: Vec<f64> = (0..360).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let ambient_lti = vec![0.1; n];
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();

    let d = fmodel.core().farm.rotor_diameters[0];
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, d * 7.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();

    let load_ti = compute_lti(&fmodel, &ambient_lti, 0.3, 10.0).unwrap();
    let voc = compute_turbine_voc(&fmodel, 2.0e-5, &ambient_lti, 0.3, 10.0, 1.0, 1.0).unwrap();

    // Python reference
    let cases = [
        (0usize, 0.100000, 0.100000, 6.144, 6.144),
        (269, 0.100000, 0.162440, 6.144, 10.980),
        (270, 0.100000, 0.162440, 6.144, 11.039),
    ];
    for (i, lti0, lti1, voc0, voc1) in cases {
        assert!(
            (load_ti[[i, 0]] - lti0).abs() < 1e-6,
            "LTI of turbine 0 at {} deg: {}",
            wind_directions[i],
            load_ti[[i, 0]]
        );
        assert!((load_ti[[i, 1]] - lti1).abs() < 1e-6, "LTI of turbine 1");
        assert!((voc[[i, 0]] - voc0).abs() < 1e-3, "VOC of turbine 0");
        assert!((voc[[i, 1]] - voc1).abs() < 1e-3, "VOC of turbine 1");
    }

    // The power setpoint sweep
    let n = 50;
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    let setpoints: Vec<f64> = (0..n)
        .map(|i| 5.0e6 + (0.001 - 5.0e6) * i as f64 / (n - 1) as f64)
        .collect();
    let power_setpoints = Array2::from_shape_vec(
        (n, 2),
        setpoints
            .iter()
            .flat_map(|s| vec![*s, 1.0e12])
            .collect::<Vec<f64>>(),
    )
    .unwrap();
    let ambient_lti = vec![0.1; n];

    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, d * 7.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel
        .set_operation(None, Some(power_setpoints), None, None, None, None)
        .unwrap();
    fmodel.run().unwrap();

    let load_ti = compute_lti(&fmodel, &ambient_lti, 0.3, 10.0).unwrap();
    let voc = compute_turbine_voc(&fmodel, 2.0e-5, &ambient_lti, 0.3, 10.0, 1.0, 1.0).unwrap();

    // Python reference
    let cases = [
        (0usize, 0.100000, 0.162440, 6.144, 11.039),
        (25, 0.100000, 0.162440, 6.144, 11.039),
        (49, 0.100000, 0.100000, 0.000, 6.144),
    ];
    for (i, lti0, lti1, voc0, voc1) in cases {
        assert!((load_ti[[i, 0]] - lti0).abs() < 1e-6, "LTI of turbine 0");
        assert!((load_ti[[i, 1]] - lti1).abs() < 1e-6, "LTI of turbine 1");
        assert!((voc[[i, 0]] - voc0).abs() < 1e-3, "VOC of turbine 0");
        assert!((voc[[i, 1]] - voc1).abs() < 1e-3, "VOC of turbine 1");
    }

    // The load ambient TI sweep
    let ambient_lti: Vec<f64> = (0..n)
        .map(|i| 0.05 + 0.2 * i as f64 / (n - 1) as f64)
        .collect();
    fmodel
        .set_operation(
            None,
            Some(Array2::from_elem((n, 2), 1.0e12)),
            None,
            None,
            None,
            None,
        )
        .unwrap();
    fmodel.run().unwrap();

    let load_ti = compute_lti(&fmodel, &ambient_lti, 0.3, 10.0).unwrap();
    let voc = compute_turbine_voc(&fmodel, 2.0e-5, &ambient_lti, 0.3, 10.0, 1.0, 1.0).unwrap();

    // Python reference
    let cases = [
        (0usize, 0.050000, 0.137428, 3.072, 9.339),
        (25, 0.152041, 0.198754, 9.342, 13.507),
        (49, 0.250000, 0.280868, 15.361, 19.087),
    ];
    for (i, lti0, lti1, voc0, voc1) in cases {
        assert!((load_ti[[i, 0]] - lti0).abs() < 1e-6, "LTI of turbine 0");
        assert!((load_ti[[i, 1]] - lti1).abs() < 1e-6, "LTI of turbine 1");
        assert!((voc[[i, 0]] - voc0).abs() < 1e-3, "VOC of turbine 0");
        assert!((voc[[i, 1]] - voc1).abs() < 1e-3, "VOC of turbine 1");
    }
}

#[test]
fn test_load_optimization_revenue_matches_python() {
    // Python: examples_load_optimization/002_row_opt_example.py
    let d = 126.0;
    let n_turbines = 3;
    let a = 4.0e-6;

    let mut fmodel = FlorisModel::from_default().unwrap();
    fmodel
        .set_layout(
            &Array1::from_iter((0..n_turbines).map(|i| i as f64 * d * 7.0)),
            &Array1::from_vec(vec![0.0; n_turbines]),
        )
        .unwrap();
    fmodel.set_operation_model("simple-derating").unwrap();

    let mut wind_directions: Vec<f64> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    let mut ambient_lti: Vec<f64> = Vec::new();
    for w_i in [240.0, 270.0] {
        for t_i in [0.1, 0.25] {
            for v_i in [1.0e-5, 1.0e-6] {
                wind_directions.push(w_i);
                values.push(v_i);
                ambient_lti.push(t_i);
            }
        }
    }
    let n = wind_directions.len();

    let mut time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    time_series.values = Array1::from_vec(values.clone());

    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();
    fmodel
        .set_operation(
            None,
            Some(Array2::from_elem((n, n_turbines), 5.0e6)),
            None,
            None,
            None,
            None,
        )
        .unwrap();
    fmodel.run().unwrap();

    let net_revenue = compute_net_revenue(&fmodel, a, &ambient_lti, 0.3, 10.0, 1.0, 1.0).unwrap();
    let farm_voc = compute_farm_voc(&fmodel, a, &ambient_lti, 0.3, 10.0, 1.0, 1.0).unwrap();
    let farm_revenue = compute_farm_revenue(&fmodel).unwrap();

    // Python reference (revenue, VOC, net revenue)
    let cases = [
        (0usize, 52.619, 3.687, 48.932),
        (2, 52.619, 9.216, 43.402),
        (6, 31.976, 10.684, 21.291),
        (7, 3.198, 10.684, -7.487),
    ];
    for (i, revenue, voc, net) in cases {
        assert!(
            (farm_revenue[i] - revenue).abs() < 1e-3,
            "farm revenue of condition {i}: {} (expected {revenue})",
            farm_revenue[i]
        );
        assert!(
            (farm_voc[i] - voc).abs() < 1e-3,
            "farm VOC of condition {i}"
        );
        assert!(
            (net_revenue[i] - net).abs() < 1e-3,
            "net revenue of condition {i}: {} (expected {net})",
            net_revenue[i]
        );
    }

    // The optimized power setpoints (Python reference: the derating search
    // picks the same levels)
    let power_setpoint_levels = linspace(POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED, 10);
    let (opt_power_setpoints, _) = optimize_power_setpoints(
        &mut fmodel,
        a,
        &ambient_lti,
        0.3,
        10.0,
        1.0,
        1.0,
        Some(&Array2::from_elem((n, n_turbines), 5.0e6)),
        &power_setpoint_levels,
    )
    .unwrap();

    // Python reference (the three power setpoints of every condition; the
    // disabled level is POWER_SETPOINT_DISABLED = 0.001 W)
    let expected_setpoints = [
        [5000000.0, 5000000.0, 5000000.0],
        [5000000.0, 5000000.0, 5000000.0],
        [5000000.0, 5000000.0, 5000000.0],
        [0.001, 0.001, 0.001],
        [5000000.0, 0.001, 5000000.0],
        [0.001, 0.001, 5000000.0],
        [5000000.0, 0.001, 5000000.0],
        [0.001, 0.001, 0.001],
    ];
    for (i, expected) in expected_setpoints.iter().enumerate() {
        for t in 0..n_turbines {
            assert!(
                (opt_power_setpoints[[i, t]] - expected[t]).abs() < 1e-6,
                "power setpoint of turbine {t} in condition {i}: {} W (expected {})",
                opt_power_setpoints[[i, t]],
                expected[t]
            );
        }
    }
}
