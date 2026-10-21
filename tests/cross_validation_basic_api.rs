//! Cross-validation of the basic API examples against the Python FLORIS
//! reference values.
//!
//! Reference values come from `floris/examples/`:
//! - `005_getting_power.py` (the turbine and farm powers)
//! - `007_sweeping_variables.py` (the variable sweeps)
//! - `010_compare_farm_power_with_neighbor.py` (the farm power with a
//!   neighboring farm)

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};
use florus::wind_data::TimeSeries;

#[test]
fn test_getting_power_matches_python() {
    // Python: examples/005_getting_power.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 126.0 * 5.0, 126.0 * 10.0]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .unwrap();

    let wind_directions: Vec<f64> = (250..290).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![9.9; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();

    let turbine_powers = fmodel.get_turbine_powers().unwrap();
    let farm_power = fmodel.get_farm_power().unwrap();
    fmodel.run_no_wake().unwrap();
    let farm_power_no_wake = fmodel.get_farm_power().unwrap();

    // Python reference (kW)
    let cases = [
        (0usize, [3325.1, 3324.8, 3324.8], 9974.7, 9975.4, 0.01),
        (20, [3325.1, 879.5, 1019.8], 5224.4, 9975.4, 47.63),
        (39, [3325.1, 3324.3, 3324.4], 9973.8, 9975.4, 0.02),
    ];
    for (i, powers, farm, farm_no_wake, wake_losses) in cases {
        for (t, power) in powers.iter().enumerate() {
            assert!(
                (turbine_powers[[i, t]] / 1000.0 - power).abs() < 0.1,
                "condition {i}, turbine {t}: {} kW (expected {power})",
                turbine_powers[[i, t]] / 1000.0
            );
        }
        assert!((farm_power[i] / 1000.0 - farm).abs() < 0.1, "farm power");
        assert!(
            (farm_power_no_wake[i] / 1000.0 - farm_no_wake).abs() < 0.1,
            "farm power without wakes"
        );
        let losses = 100.0 * (farm_power_no_wake[i] - farm_power[i]) / farm_power_no_wake[i];
        assert!((losses - wake_losses).abs() < 0.01, "wake losses");
    }
}

#[test]
fn test_sweeping_variables_matches_python() {
    // Python: examples/007_sweeping_variables.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 126.0 * 5.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();

    // The turbulence intensity sweep
    let turbulence_intensities: Vec<f64> = (0..17).map(|i| 0.03 + 0.01 * i as f64).collect();
    let n = turbulence_intensities.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n]),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(turbulence_intensities.clone()),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();

    // Python reference (kW)
    let ti_cases = [(0usize, 331.1), (8, 770.9), (16, 1128.8)];
    for (i, downstream) in ti_cases {
        assert!(
            (powers[[i, 0]] / 1000.0 - 1754.0).abs() < 0.1,
            "upstream power"
        );
        assert!(
            (powers[[i, 1]] / 1000.0 - downstream).abs() < 0.1,
            "downstream power at TI {}: {} kW",
            turbulence_intensities[i],
            powers[[i, 1]] / 1000.0
        );
    }

    // The upstream yaw angle sweep
    let n_yaw = 100;
    let yaw_angles_upstream: Vec<f64> = (0..n_yaw)
        .map(|i| -30.0 + 60.0 * i as f64 / (n_yaw - 1) as f64)
        .collect();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_yaw]),
        Array1::from_vec(vec![8.0; n_yaw]),
        Array1::from_vec(vec![0.06; n_yaw]),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    let mut yaw_angles = Array2::zeros((n_yaw, 2));
    for (i, yaw) in yaw_angles_upstream.iter().enumerate() {
        yaw_angles[[i, 0]] = *yaw;
    }
    fmodel.set_yaw_angles(yaw_angles).unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();

    // Python reference (kW)
    let yaw_cases = [
        (0usize, 1339.9, 1010.0),
        (50, 1753.9, 436.9),
        (99, 1339.9, 1023.4),
    ];
    for (i, upstream, downstream) in yaw_cases {
        assert!(
            (powers[[i, 0]] / 1000.0 - upstream).abs() < 0.1,
            "upstream power at yaw {}",
            yaw_angles_upstream[i]
        );
        assert!(
            (powers[[i, 1]] / 1000.0 - downstream).abs() < 0.1,
            "downstream power at yaw {}",
            yaw_angles_upstream[i]
        );
    }

    // The disabled turbine combinations
    fmodel.reset_operation();
    fmodel.set_operation_model("mixed").unwrap();
    let disable_combinations = [[false, false], [true, false], [false, true], [true, true]];
    let n_combinations = disable_combinations.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0; n_combinations]),
        Array1::from_vec(vec![8.0; n_combinations]),
        Array1::from_vec(vec![0.06; n_combinations]),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((n_combinations, 2), false);
    for (i, combination) in disable_combinations.iter().enumerate() {
        disable_turbines[[i, 0]] = combination[0];
        disable_turbines[[i, 1]] = combination[1];
    }
    fmodel
        .set_operation(None, None, None, None, None, Some(disable_turbines))
        .unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();

    // Python reference (kW)
    let disable_cases = [
        (0usize, 1754.0, 436.4),
        (1, 0.0, 1754.0),
        (2, 1754.0, 0.0),
        (3, 0.0, 0.0),
    ];
    for (i, upstream, downstream) in disable_cases {
        assert!(
            (powers[[i, 0]] / 1000.0 - upstream).abs() < 0.1,
            "upstream power"
        );
        assert!(
            (powers[[i, 1]] / 1000.0 - downstream).abs() < 0.1,
            "downstream power"
        );
    }
}

#[test]
fn test_neighbor_farm_power_matches_python() {
    // Python: examples/010_compare_farm_power_with_neighbor.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    let d = 126.0;
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, d * 6.0, 0.0, d * 6.0]),
            &Array1::from_vec(vec![0.0, 0.0, d * 3.0, d * 3.0]),
        )
        .unwrap();

    let wind_directions: Vec<f64> = (0..90).map(|i| 4.0 * i as f64).collect();
    let n = wind_directions.len();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; n]),
            Array1::from_vec(wind_directions.clone()),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();
    fmodel.run().unwrap();
    let farm_power_base = fmodel.get_farm_power().unwrap();

    // Add the neighboring farm to the east
    fmodel
        .set_layout(
            &Array1::from_vec(vec![
                0.0,
                d * 6.0,
                0.0,
                d * 6.0,
                d * 12.0,
                d * 15.0,
                d * 12.0,
                d * 15.0,
            ]),
            &Array1::from_vec(vec![0.0, 0.0, d * 3.0, d * 3.0, 0.0, 0.0, d * 3.0, d * 3.0]),
        )
        .unwrap();

    let mut turbine_weights = Array2::zeros((n, 8));
    for fi in 0..n {
        for t in 0..4 {
            turbine_weights[[fi, t]] = 1.0;
        }
    }
    fmodel.run().unwrap();
    let farm_power_neighbor = fmodel
        .get_farm_power_weighted(Some(&turbine_weights))
        .unwrap();

    // Python reference (kW); the base farm differs by up to 0.4 kW in a few
    // wind directions (a wake-interaction difference of 0.01 %)
    let cases = [
        (0usize, 4125.5, 4125.5),
        (1, 4418.9, 4420.7),
        (89, 4416.9, 4416.9),
    ];
    for (i, base, neighbor) in cases {
        assert!(
            (farm_power_base[i] / 1000.0 - base).abs() < 0.5,
            "base farm power at {} deg: {} kW (expected {base})",
            wind_directions[i],
            farm_power_base[i] / 1000.0
        );
        assert!(
            (farm_power_neighbor[i] / 1000.0 - neighbor).abs() < 0.5,
            "farm power with the neighbor at {} deg: {} kW (expected {neighbor})",
            wind_directions[i],
            farm_power_neighbor[i] / 1000.0
        );
    }
}
