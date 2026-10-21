//! Cross-validation of the yaw optimization and wind data comparison examples
//! against the Python FLORIS reference values.
//!
//! Reference values come from:
//! - `floris/examples/examples_control_optimization/001_opt_yaw_single_ws.py`
//! - `floris/examples/examples_wind_data/001_wind_data_comparisons.py`

use florus::floris_model::FlorisModel;
use florus::optimization::yaw_optimization::yaw_optimizer_sr::YawOptimizationSR;
use florus::optimization::yaw_optimization::YawOptimization;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

#[test]
fn test_yaw_optimization_single_ws_matches_python() {
    // Python: examples_control_optimization/001_opt_yaw_single_ws.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();

    let wind_directions: Vec<f64> = (0..120).map(|i| 3.0 * i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();

    let d = 126.0;
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();

    let mut yaw_opt = YawOptimizationSR::new();
    let result = yaw_opt.optimize(&mut fmodel, None).unwrap();

    // Python reference: the yaw angles and the farm powers
    let cases = [
        (0usize, [0.0, 0.0, 0.0], 5261863.4, 5261863.4),
        (60, [0.0, 0.0, 0.0], 5261863.4, 5261863.4),
        (119, [0.0, 0.0, 0.0], 5261863.4, 5261863.4),
    ];
    for (i, yaw, farm_power_opt, farm_power_baseline) in cases {
        for t in 0..3 {
            assert!(
                (result.yaw_angles[[i, t]] - yaw[t]).abs() < 1e-6,
                "wind direction {}: yaw of turbine {t} = {}",
                wind_directions[i],
                result.yaw_angles[[i, t]]
            );
        }
        assert!(
            (result.farm_powers_opt[i] - farm_power_opt).abs() < 1.0,
            "optimized farm power at {} deg: {} W",
            wind_directions[i],
            result.farm_powers_opt[i]
        );
        assert!(
            (result.farm_powers_baseline[i] - farm_power_baseline).abs() < 1.0,
            "baseline farm power at {} deg: {} W",
            wind_directions[i],
            result.farm_powers_baseline[i]
        );
    }
}

#[test]
fn test_wind_data_comparisons_match_python() {
    // Python: examples_wind_data/001_wind_data_comparisons.py
    let n = 500;
    let wind_directions: Vec<f64> = (0..n)
        .map(|i| 270.0 + 0.04 * ((i % 37) as f64 - 18.0))
        .collect();
    let wind_speeds: Vec<f64> = (0..n)
        .map(|i| (8.0 + 0.16 * ((i % 23) as f64 - 11.0)).clamp(3.0, 50.0))
        .collect();
    let turbulence_intensities: Vec<f64> = (0..n)
        .map(|i| (0.1 + 0.001 * ((i % 17) as f64 - 8.0)).clamp(0.0, 0.25))
        .collect();
    let values: Vec<f64> = (0..n)
        .map(|i| (25.0 + 0.2 * ((i % 29) as f64 - 14.0)).clamp(0.0, 100.0))
        .collect();

    let mut time_series = TimeSeries::new(
        Array1::from_vec(wind_directions),
        Array1::from_vec(wind_speeds),
        Array1::from_vec(turbulence_intensities),
    )
    .unwrap();
    time_series.values = Array1::from_vec(values);
    let wind_rose = time_series.to_wind_rose(5.0, 1.0);

    // The AEP of the time series
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();
    let aep_time_series = fmodel.get_farm_aep(&time_series, 8760.0).unwrap();

    // The AEP of the wind rose
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&wind_rose).unwrap();
    fmodel.run().unwrap();
    let aep_wind_rose = fmodel.get_farm_aep(&wind_rose, 8760.0).unwrap();

    // Python reference (GWh)
    assert!(
        (aep_time_series / 1.0e9 - 21.043).abs() < 0.001,
        "time series AEP"
    );
    assert!(
        (aep_wind_rose / 1.0e9 - 21.031).abs() < 0.001,
        "wind rose AEP"
    );
}
