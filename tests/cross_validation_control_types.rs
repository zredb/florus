//! Cross-validation of the control-type examples against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/examples_control_types/`:
//! - `001_derating_control.py` (simple-derating operation model)
//! - `002_disable_turbines.py` (mixed operation model, disabled turbines)
//! - `004_helix_active_wake_mixing.py` (awc operation model with helix)
//! - `005_peak_shaving.py` (peak-shaving operation model)

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::YawOptimization;
use florus::types::{Array1, Array2};
use florus::wind_data::TimeSeries;

#[test]
fn test_yaw_optimization_single_ws_matches_python() {
    // Python: examples_control_optimization/001_opt_yaw_single_ws.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    let d = 126.0;
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();

    let wind_directions: Vec<f64> = (0..120).map(|i| 3.0 * i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();

    let mut yaw_opt = YawOptimizationSR::new();
    let result = yaw_opt.optimize(&mut fmodel, None).unwrap();

    // Python reference: the yaw angles stay at zero for most of the wind
    // directions (the turbines are in a row) and the farm power is
    // 5261863.3775 W; at 66 and 69 degrees the orientation of the farm
    // changes the wake losses slightly, and at 270 degrees the wind is
    // perpendicular to the row, so the optimizer yaws the turbines.
    let exact = [
        (0usize, 5261863.3775, [0.0, 0.0, 0.0]),
        (22, 5261858.6925, [0.0, 0.0, 0.0]),
        (23, 5261734.2066, [0.0, 0.0, 0.0]),
        (40, 5261863.3757, [0.0, 0.0, 0.0]),
        (90, 3275737.7312, [25.0, 22.65625, 0.0]),
    ];
    for i in 0..n {
        assert!(
            result.farm_powers_opt[i] <= 5261863.4 + 1.0,
            "optimized farm power at {} deg: {} W",
            wind_directions[i],
            result.farm_powers_opt[i]
        );
    }
    for (i, expected, yaw) in exact {
        assert!(
            (result.farm_powers_opt[i] - expected).abs() < 0.01,
            "optimized farm power at {} deg: {} W (expected {expected})",
            wind_directions[i],
            result.farm_powers_opt[i]
        );
        for t in 0..3 {
            assert!(
                (result.yaw_angles[[i, t]] - yaw[t]).abs() < 1e-6,
                "yaw angle of turbine {t} at {} deg: {} (expected {})",
                wind_directions[i],
                result.yaw_angles[[i, t]],
                yaw[t]
            );
        }
    }
}

#[test]
fn test_derating_control_matches_python() {
    // Python: examples_control_types/001_derating_control.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel.set_operation_model("simple-derating").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();

    let n = 50;
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![10.0; n]),
            Array1::from_vec(vec![270.0; n]),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();
    fmodel.run().unwrap();
    let powers_orig = fmodel.get_turbine_powers().unwrap();

    let setpoints: Vec<f64> = (0..n)
        .map(|i| 1.0 + (6.0e6 - 1.0) * i as f64 / (n - 1) as f64)
        .collect();
    let mut power_setpoints = Array2::zeros((n, 2));
    for (i, s) in setpoints.iter().enumerate() {
        power_setpoints[[i, 0]] = *s;
        power_setpoints[[i, 1]] = *s;
    }
    fmodel
        .set_operation(None, Some(power_setpoints), None, None, None, None)
        .unwrap();
    fmodel.run().unwrap();
    let powers_derated = fmodel.get_turbine_powers().unwrap();

    // Python reference (kW) at a few setpoints
    let cases = [
        (0usize, 0.0, 0.0, 0.0),
        (1, 122.4, 122.4, 122.4),
        (10, 1224.5, 1224.5, 1224.5),
        (25, 3061.2, 1604.3, 1604.3),
        (49, 3417.8, 1514.2, 1514.2),
    ];
    for (i, upstream, downstream, available_ds) in cases {
        assert!(
            (powers_orig[[i, 0]] / 1000.0 - 3417.8).abs() < 0.1,
            "original upstream power"
        );
        assert!(
            (powers_derated[[i, 0]] / 1000.0 - upstream).abs() < 0.1,
            "derated upstream power at setpoint {upstream}"
        );
        assert!(
            (powers_derated[[i, 1]] / 1000.0 - downstream).abs() < 0.1,
            "derated downstream power at setpoint {upstream}"
        );
        let _ = available_ds;
    }
}

#[test]
fn test_disable_turbines_matches_python() {
    // Python: examples_control_types/002_disable_turbines.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel.set_operation_model("mixed").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0, 8.0]),
            Array1::from_vec(vec![270.0, 270.0]),
            Array1::from_vec(vec![0.06, 0.06]),
        )
        .unwrap();

    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((2, 3), false);
    disable_turbines[[1, 0]] = true;
    disable_turbines[[1, 1]] = true;
    fmodel
        .set_operation(None, None, None, None, None, Some(disable_turbines))
        .unwrap();
    fmodel.run().unwrap();

    let powers = fmodel.get_turbine_powers().unwrap();
    // Python reference (kW)
    let expected = [[1753.95, 354.99, 339.56], [0.0, 0.0, 1753.95]];
    for fi in 0..2 {
        for t in 0..3 {
            assert!(
                (powers[[fi, t]] / 1000.0 - expected[fi][t]).abs() < 0.01,
                "condition {fi}, turbine {t}: {} kW",
                powers[[fi, t]] / 1000.0
            );
        }
    }

    // The effective wind speeds
    let velocities = fmodel.turbine_average_velocities();
    assert!((velocities[[0, 0]] - 7.97368581).abs() < 1e-6);
    assert!((velocities[[0, 1]] - 4.78380287).abs() < 1e-6);
    assert!((velocities[[0, 2]] - 4.71560113).abs() < 1e-6);
    assert!((velocities[[1, 2]] - 7.97368581).abs() < 1e-6);
}

#[test]
fn test_yaw_and_disabling_matches_python() {
    // Python: examples_control_types/003_setting_yaw_and_disabling.py
    let mut fmodel_2 = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    let mut fmodel_3 = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel_2.set_operation_model("mixed").unwrap();
    fmodel_3.set_operation_model("mixed").unwrap();
    fmodel_2
        .set_layout(
            &Array1::from_vec(vec![0.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel_3
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();

    let n = 50;
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    fmodel_2.set_wind_data(&time_series).unwrap();
    fmodel_3.set_wind_data(&time_series).unwrap();

    // The yaw angle of the most upstream turbine is swept from -20 to 20
    let upstream_yaw: Vec<f64> = (0..n)
        .map(|i| -20.0 + 40.0 * i as f64 / (n - 1) as f64)
        .collect();
    let mut yaw_2 = Array2::zeros((n, 2));
    let mut yaw_3 = Array2::zeros((n, 3));
    for (i, yaw) in upstream_yaw.iter().enumerate() {
        yaw_2[[i, 0]] = *yaw;
        yaw_3[[i, 0]] = *yaw;
    }

    // The middle turbine of the three-turbine case is disabled
    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((n, 3), false);
    for i in 0..n {
        disable_turbines[[i, 1]] = true;
    }

    fmodel_2
        .set_operation(Some(yaw_2), None, None, None, None, None)
        .unwrap();
    fmodel_3
        .set_operation(Some(yaw_3), None, None, None, None, Some(disable_turbines))
        .unwrap();

    fmodel_2.run().unwrap();
    fmodel_3.run().unwrap();
    let powers_2 = fmodel_2.get_turbine_powers().unwrap();
    let powers_3 = fmodel_3.get_turbine_powers().unwrap();

    // Python reference (kW): yaw, 2-turbine upstream/downstream, 3-turbine
    // upstream/downstream
    let cases = [
        (0usize, -20.0, 1561.3, 1090.7, 1561.3, 1090.7),
        (1, -19.18, 1576.1, 1073.6, 1576.1, 1073.6),
        (25, 0.41, 1753.9, 769.0, 1753.9, 769.0),
        (49, 20.0, 1561.3, 1106.9, 1561.3, 1106.9),
    ];
    for (i, yaw, p2u, p2d, p3u, p3d) in cases {
        assert!((upstream_yaw[i] - yaw).abs() < 0.01);
        assert!(
            (powers_2[[i, 0]] / 1000.0 - p2u).abs() < 0.1,
            "2-turbine upstream at {yaw} deg"
        );
        assert!(
            (powers_2[[i, 1]] / 1000.0 - p2d).abs() < 0.1,
            "2-turbine downstream at {yaw} deg"
        );
        assert!(
            (powers_3[[i, 0]] / 1000.0 - p3u).abs() < 0.1,
            "3-turbine upstream at {yaw} deg"
        );
        assert!(
            (powers_3[[i, 2]] / 1000.0 - p3d).abs() < 0.1,
            "3-turbine downstream at {yaw} deg"
        );
    }
}

#[test]
fn test_helix_awc_matches_python() {
    // Python: examples_control_types/004_helix_active_wake_mixing.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/emgauss_helix.yaml").unwrap();
    fmodel.set_operation_model("awc").unwrap();

    let n = 1;
    let awc_modes = ndarray::Array2::<String>::from_shape_vec(
        (n, 4),
        vec![
            "helix".to_string(),
            "baseline".to_string(),
            "baseline".to_string(),
            "baseline".to_string(),
        ],
    )
    .unwrap();
    let awc_amplitudes = Array2::from_shape_vec((n, 4), vec![2.5, 0.0, 0.0, 0.0]).unwrap();

    let d = 240.0;
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 4.0 * d, 0.0, 4.0 * d]),
            &Array1::from_vec(vec![0.0, 0.0, -3.0 * d, -3.0 * d]),
        )
        .unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; n]),
            Array1::from_vec(vec![270.0; n]),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();
    fmodel
        .set_operation(
            None,
            None,
            Some(awc_modes),
            Some(awc_amplitudes),
            None,
            None,
        )
        .unwrap();
    fmodel.run().unwrap();

    // Python reference (kW)
    let powers = fmodel.get_turbine_powers().unwrap();
    let expected = [6133.2, 2827.1, 6301.1, 2076.6];
    for t in 0..4 {
        assert!(
            (powers[[0, t]] / 1000.0 - expected[t]).abs() < 0.1,
            "turbine {t}: {} kW",
            powers[[0, t]] / 1000.0
        );
    }
}

#[test]
fn test_peak_shaving_matches_python() {
    // Python: examples_control_types/005_peak_shaving.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 1000.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();

    let wind_speeds: Vec<f64> = (0..100).map(|i| 30.0 * i as f64 / 99.0).collect();
    let n = wind_speeds.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(wind_speeds.clone()),
        Array1::from_vec(vec![0.10; n]),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();

    fmodel.set_operation_model("simple").unwrap();
    fmodel.run().unwrap();
    let powers_base = fmodel.get_turbine_powers().unwrap() / 1000.0;
    let thrust_base = fmodel.get_turbine_thrust_coefficients();

    fmodel.set_operation_model("peak-shaving").unwrap();
    fmodel.run().unwrap();
    let powers_peak = fmodel.get_turbine_powers().unwrap() / 1000.0;
    let thrust_peak = fmodel.get_turbine_thrust_coefficients();

    // Python reference (wind speed, Ct simple, Ct peak shaving, power simple,
    // power peak shaving)
    let cases = [
        (20usize, 6.06, 0.8590, 0.8590, 755.9, 755.9),
        (40, 12.12, 0.5312, 0.5215, 5000.0, 4928.9),
        (50, 15.15, 0.2440, 0.2440, 5000.0, 5000.0),
        (60, 18.18, 0.1410, 0.1410, 5000.0, 5000.0),
    ];
    for (i, ws, ct_base, ct_peak, p_base, p_peak) in cases {
        assert!(
            (wind_speeds[i] - ws).abs() < 0.01,
            "wind speed {}",
            wind_speeds[i]
        );
        assert!(
            (thrust_base[[i, 0]] - ct_base).abs() < 1e-4,
            "Ct without peak shaving at {ws} m/s: {}",
            thrust_base[[i, 0]]
        );
        assert!(
            (thrust_peak[[i, 0]] - ct_peak).abs() < 1e-4,
            "Ct with peak shaving at {ws} m/s: {}",
            thrust_peak[[i, 0]]
        );
        assert!(
            (powers_base[[i, 0]] - p_base).abs() < 0.1,
            "power without peak shaving at {ws} m/s: {}",
            powers_base[[i, 0]]
        );
        assert!(
            (powers_peak[[i, 0]] - p_peak).abs() < 0.1,
            "power with peak shaving at {ws} m/s: {}",
            powers_peak[[i, 0]]
        );
    }
}
