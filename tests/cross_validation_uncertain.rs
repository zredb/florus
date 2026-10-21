//! Cross-validation of the uncertain model against the Python FLORIS reference
//! values.
//!
//! Reference values come from `floris/examples/008_uncertain_models.py`.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};
use florus::uncertain::UncertainFlorisModel;
use florus::wind_data::TimeSeries;

#[test]
fn test_parallel_uncertain_model_matches_python() {
    // Python: examples_uncertain/003_uncertain_model_with_parallelization.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    let mut ufmodel = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml").unwrap(),
        3.0,
        None,
    )
    .unwrap();
    let mut pufmodel = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml").unwrap(),
        3.0,
        None,
    )
    .unwrap();

    let wind_directions: Vec<f64> = (240..300).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();

    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, d * 6.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);

    fmodel.set_layout(&layout_x, &layout_y).unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    ufmodel.set_layout(&layout_x, &layout_y).unwrap();
    ufmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; n]),
            Array1::from_vec(wind_directions.clone()),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();
    pufmodel.set_layout(&layout_x, &layout_y).unwrap();
    pufmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; n]),
            Array1::from_vec(wind_directions.clone()),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();

    fmodel.run().unwrap();
    ufmodel.run().unwrap();
    pufmodel.run_parallel(2).unwrap();

    let nominal = fmodel.get_turbine_powers().unwrap();
    let serial = ufmodel.get_turbine_powers().unwrap();
    let parallel = pufmodel.get_turbine_powers().unwrap();

    // The nominal and the uncertain powers (Python reference)
    let cases = [
        (0usize, 1753.95, 1753.95, 1753.95, 1753.95),
        (30, 1753.95, 575.21, 1753.95, 754.06),
        (59, 1753.95, 1753.95, 1753.95, 1753.95),
    ];
    for (i, nom0, nom1, unc0, unc1) in cases {
        assert!((nominal[[i, 0]] / 1000.0 - nom0).abs() < 0.01, "nominal t0");
        assert!((nominal[[i, 1]] / 1000.0 - nom1).abs() < 0.01, "nominal t1");
        assert!((serial[[i, 0]] / 1000.0 - unc0).abs() < 0.01, "serial t0");
        assert!((serial[[i, 1]] / 1000.0 - unc1).abs() < 0.01, "serial t1");
    }

    // The parallel results are identical to the serial ones
    let difference = (&serial - &parallel).mapv(|v| v.abs());
    assert!(
        difference.iter().copied().fold(0.0, f64::max) < 1e-6,
        "the parallel and the serial results differ"
    );
}

#[test]
fn test_fix_yaw_to_nominal_direction() {
    // The yaw angles follow the wind direction when they are fixed to the
    // nominal direction (Python `fix_yaw_to_nominal_direction`).
    let mut model = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    model
        .set_layout(
            &Array1::from_vec(vec![0.0, 630.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();
    model
        .set_wind_conditions(
            Array1::from_vec(vec![8.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .unwrap();
    model
        .set_yaw_angles(Array2::from_shape_vec((1, 2), vec![10.0, 0.0]).unwrap())
        .unwrap();

    let mut uncertain = UncertainFlorisModel::new(model, 3.0, Some(vec![-6.0, 0.0, 6.0])).unwrap();
    uncertain.fix_yaw_to_nominal_direction = true;
    uncertain.expand_public().unwrap();

    // The yaw angles of the expanded conditions are shifted with the wind
    // direction
    let expanded = uncertain.expanded_model().unwrap();
    let yaw_angles = &expanded.core().farm.yaw_angles;
    let wind_directions = &expanded.core().flow_field.wind_directions;
    for i in 0..yaw_angles.shape()[0] {
        let offset = (wind_directions[i] - 270.0).rem_euclid(360.0);
        let expected = if offset > 180.0 {
            10.0 + offset - 360.0
        } else {
            10.0 + offset
        };
        assert!(
            (yaw_angles[[i, 0]] - expected).abs() < 1e-9,
            "yaw angle of condition {i}: {} (expected {expected})",
            yaw_angles[[i, 0]]
        );
    }
}

#[test]
fn test_uncertain_model_matches_python() {
    // Python: examples/008_uncertain_models.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    let mut ufmodel_3 = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml").unwrap(),
        3.0,
        None,
    )
    .unwrap();
    let mut ufmodel_5 = UncertainFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml").unwrap(),
        5.0,
        None,
    )
    .unwrap();

    // The wind direction is swept while the wind speed and the turbulence
    // intensity are held constant
    let wind_directions: Vec<f64> = (240..300).map(|i| i as f64).collect();
    let n = wind_directions.len();
    let time_series = TimeSeries::new(
        Array1::from_vec(wind_directions.clone()),
        Array1::from_vec(vec![8.0; n]),
        Array1::from_vec(vec![0.06; n]),
    )
    .unwrap();

    // A two turbine farm
    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, d * 6.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);

    fmodel.set_layout(&layout_x, &layout_y).unwrap();
    fmodel.set_wind_data(&time_series).unwrap();
    ufmodel_3.set_layout(&layout_x, &layout_y).unwrap();
    ufmodel_5.set_layout(&layout_x, &layout_y).unwrap();
    ufmodel_3
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; n]),
            Array1::from_vec(wind_directions.clone()),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();
    ufmodel_5
        .set_wind_conditions(
            Array1::from_vec(vec![8.0; n]),
            Array1::from_vec(wind_directions.clone()),
            Array1::from_vec(vec![0.06; n]),
        )
        .unwrap();

    fmodel.run().unwrap();
    ufmodel_3.run().unwrap();
    ufmodel_5.run().unwrap();

    // The sizes of the expanded and unique cases (Python prints them with
    // verbose=True)
    assert_eq!(ufmodel_3.n_unexpanded(), 60);
    assert_eq!(ufmodel_3.n_expanded(), 300);
    assert_eq!(ufmodel_3.n_unique(), 72);

    let nominal = fmodel.get_turbine_powers().unwrap();
    let uncertain_3 = ufmodel_3.get_turbine_powers().unwrap();
    let uncertain_5 = ufmodel_5.get_turbine_powers().unwrap();

    // Python reference (kW): the wind direction, the nominal and the
    // uncertain (3 and 5 degrees) turbine powers
    let cases = [
        (0usize, 1753.95, 1753.95, 1753.95, 1753.95, 1753.95, 1753.95),
        (10, 1753.95, 1753.95, 1753.95, 1753.58, 1753.95, 1748.30),
        (16, 1753.95, 1748.52, 1753.95, 1725.45, 1753.95, 1668.94),
        (30, 1753.95, 575.21, 1753.95, 754.06, 1753.95, 948.61),
        (59, 1753.95, 1753.95, 1753.95, 1753.95, 1753.95, 1753.95),
    ];
    for (i, nom0, nom1, u3_0, u3_1, u5_0, u5_1) in cases {
        let values = [
            (nominal[[i, 0]] / 1000.0, nom0),
            (nominal[[i, 1]] / 1000.0, nom1),
            (uncertain_3[[i, 0]] / 1000.0, u3_0),
            (uncertain_3[[i, 1]] / 1000.0, u3_1),
            (uncertain_5[[i, 0]] / 1000.0, u5_0),
            (uncertain_5[[i, 1]] / 1000.0, u5_1),
        ];
        for (value, expected) in values {
            assert!(
                (value - expected).abs() < 0.01,
                "wind direction {}: {} kW (expected {expected})",
                wind_directions[i],
                value
            );
        }
    }
}

/// The approximate model (Python's `ApproxFlorisModel`) is an uncertain
/// model with a single wind direction sample point and a reduced wind speed
/// resolution: it has to reproduce the expected power of the full model
/// closely, while solving far fewer conditions.
#[test]
fn test_approx_floris_model_reproduces_the_expected_power() {
    use florus::uncertain::ApproxFlorisModel;

    let d = 126.0;
    let layout_x = Array1::from_vec(vec![0.0, 5.0 * d, 10.0 * d]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);

    // A set of wind speeds spread over the operating range; the speeds are
    // 0.375 m/s apart, so that they are merged in pairs by the 0.5 m/s
    // resolution of the approximate model
    let wind_speeds = Array1::from_iter((0..40).map(|k| 5.0 + 0.375 * k as f64));
    let n = wind_speeds.len();
    let time_series = TimeSeries::new(
        Array1::from_elem(n, 270.0),
        wind_speeds,
        Array1::from_elem(n, 0.06),
    )
    .expect("time series");

    let mut model = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    model.set_layout(&layout_x, &layout_y).expect("layout");
    model.set_wind_data(&time_series).expect("wind data");
    model.run().expect("run");

    let mut approx = ApproxFlorisModel::new(
        FlorisModel::from_file("examples/inputs/gch.yaml").expect("model"),
        0.5,
    )
    .expect("approx model");
    approx.set_layout(&layout_x, &layout_y).expect("layout");
    approx.set_wind_data(&time_series).expect("wind data");
    approx.run().expect("run");

    // The approximate model solves far fewer conditions
    assert!(
        approx.model.n_unique() < n,
        "the approximate model solves {} conditions (expected fewer than {n})",
        approx.model.n_unique()
    );

    // The expected power is close to the power of the full model
    let exact = model
        .get_expected_farm_power(None, None)
        .expect("expected power");
    let approximate = approx
        .get_expected_farm_power(None)
        .expect("expected power");
    let difference = (approximate - exact).abs() / exact;
    assert!(
        difference < 0.01,
        "the expected power differs by {:.3} %: {approximate} W vs {exact} W",
        100.0 * difference
    );
}
