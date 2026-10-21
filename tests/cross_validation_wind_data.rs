//! Cross-validation of the wind-data examples against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/examples_wind_data/`:
//! - `002_generate_ti.py` (the turbulence intensity of a wind rose)
//! - `003_generate_value.py` (the value of a wind rose)

use florus::types::{Array1, Array2};
use florus::wind_data::{TIParams, TimeSeries, WindRose};

#[test]
fn test_iec_turbulence_intensity_matches_python() {
    // Python: examples_wind_data/002_generate_ti.py
    let wind_directions = Array1::from_vec(vec![250.0, 260.0, 270.0]);
    let wind_speeds = Array1::from_vec(vec![5.0, 6.0, 7.0, 8.0, 9.0, 10.0]);

    let mut wind_rose = WindRose::new(
        wind_directions,
        wind_speeds.clone(),
        Array2::from_elem((3, 6), 0.06),
        None,
        None,
        false,
        None,
        None,
    )
    .unwrap();

    // The custom function TI = 1 / wind_speed
    wind_rose.assign_ti_using_wd_ws_function(|_wd: f64, ws: f64| 1.0 / ws);
    let expected_custom = [0.200000, 0.166667, 0.142857, 0.125000, 0.111111, 0.100000];
    for (j, value) in expected_custom.iter().enumerate() {
        assert!(
            (wind_rose.ti_table[[0, j]] - value).abs() < 1e-6,
            "TI at {} m/s: {}",
            wind_speeds[j],
            wind_rose.ti_table[[0, j]]
        );
    }

    // The IEC 61400-1 method
    wind_rose.assign_ti_using_iec_method(Some(TIParams {
        iref: 0.07,
        ..Default::default()
    }));
    let expected_iec = [0.105700, 0.096833, 0.090500, 0.085750, 0.082056, 0.079100];
    for (j, value) in expected_iec.iter().enumerate() {
        assert!(
            (wind_rose.ti_table[[0, j]] - value).abs() < 1e-6,
            "TI at {} m/s: {}",
            wind_speeds[j],
            wind_rose.ti_table[[0, j]]
        );
    }
}

#[test]
fn test_value_assignment_matches_python() {
    // Python: examples_wind_data/003_generate_value.py
    let wind_directions = Array1::from_vec(vec![250.0, 260.0, 270.0]);
    let wind_speeds = Array1::from_iter((3..11).map(|i| i as f64));

    let mut wind_rose = WindRose::new(
        wind_directions,
        wind_speeds.clone(),
        Array2::from_elem((3, 8), 0.06),
        None,
        None,
        false,
        None,
        None,
    )
    .unwrap();

    // The custom function value = 100 / wind_speed
    wind_rose.assign_value_using_wd_ws_function(|_wd: f64, ws: f64| 100.0 / ws, false);
    let expected_custom = [
        33.333333, 25.000000, 20.000000, 16.666667, 14.285714, 12.500000, 11.111111, 10.000000,
    ];
    for (j, value) in expected_custom.iter().enumerate() {
        assert!(
            (wind_rose.value_table.as_ref().unwrap()[[0, j]] - value).abs() < 1e-6,
            "value at {} m/s: {}",
            wind_speeds[j],
            wind_rose.value_table.as_ref().unwrap()[[0, j]]
        );
    }

    // The piecewise linear function
    wind_rose.assign_value_piecewise_linear(1.425, 4.5, 0.0, -0.135, false, false);
    let expected_piecewise = [
        1.425000, 1.425000, 1.357500, 1.222500, 1.087500, 0.952500, 0.817500, 0.682500,
    ];
    for (j, value) in expected_piecewise.iter().enumerate() {
        assert!(
            (wind_rose.value_table.as_ref().unwrap()[[0, j]] - value).abs() < 1e-6,
            "value at {} m/s: {}",
            wind_speeds[j],
            wind_rose.value_table.as_ref().unwrap()[[0, j]]
        );
    }
}

#[test]
fn test_time_series_iec_turbulence_intensity_matches_python() {
    // Python: examples_wind_data/002_generate_ti.py (the time series part)
    let n = 100;
    let mut time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_iter((0..n).map(|i| 5.0 + 10.0 * i as f64 / (n - 1) as f64)),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    time_series.assign_ti_using_iec_method(Some(TIParams {
        iref: 0.07,
        ..Default::default()
    }));

    // The IEC values at the first and the last wind speed
    assert!((time_series.turbulence_intensities[0] - 0.105700).abs() < 1e-6);
    let last = n - 1;
    // TI at 15 m/s is Iref
    assert!((time_series.turbulence_intensities[last] - 0.070233).abs() < 1e-5);
}
