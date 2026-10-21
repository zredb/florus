//! Cross-validation of the wind data objects example against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/003_wind_data_objects.py`
//! (with deterministic inputs instead of the random signals of the example).

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::TimeSeries;

#[test]
fn test_wind_data_objects_match_python() {
    let n = 100;
    let wind_speeds: Vec<f64> = (0..n).map(|i| 8.0 + 0.1 * i as f64).collect();
    let wind_directions: Vec<f64> = (0..n).map(|i| 270.0 + 0.2 * i as f64).collect();
    let turbulence_intensities: Vec<f64> = (0..n).map(|i| 0.06 + 0.0001 * i as f64).collect();

    let mut time_series = TimeSeries::new(
        Array1::from_vec(wind_directions),
        Array1::from_vec(wind_speeds),
        Array1::from_vec(turbulence_intensities),
    )
    .unwrap();
    time_series.values = Array1::from_iter((0..n).map(|i| 25.0 + 0.05 * i as f64));

    // The time series converted into a wind rose
    let wind_rose = time_series.to_wind_rose(5.0, 1.0);
    assert_eq!(wind_rose.wind_directions.len(), 72);
    assert_eq!(wind_rose.wind_speeds.len(), 50);

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 630.0, 1260.0]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();

    // The powers with the time series
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();
    // Python reference (kW)
    let expected = [
        [1754.0, 436.4, 506.7],
        [1826.0, 459.8, 531.9],
        [1900.5, 485.3, 559.6],
    ];
    for (i, row) in expected.iter().enumerate() {
        for (t, power) in row.iter().enumerate() {
            assert!(
                (powers[[i, t]] / 1000.0 - power).abs() < 0.1,
                "time series, condition {i}, turbine {t}: {} kW (expected {power})",
                powers[[i, t]] / 1000.0
            );
        }
    }

    // The powers with the wind rose
    fmodel.set_wind_data(&wind_rose).unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();
    let expected = [[1754.0, 438.2, 508.1], [2496.4, 658.2, 745.9]];
    for (i, row) in expected.iter().enumerate() {
        for (t, power) in row.iter().enumerate() {
            assert!(
                (powers[[i, t]] / 1000.0 - power).abs() < 0.1,
                "wind rose, condition {i}, turbine {t}: {} kW (expected {power})",
                powers[[i, t]] / 1000.0
            );
        }
    }
}
