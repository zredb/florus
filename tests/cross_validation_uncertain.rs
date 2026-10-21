//! Cross-validation of the uncertain model against the Python FLORIS reference
//! values.
//!
//! Reference values come from `floris/examples/008_uncertain_models.py`.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::uncertain::UncertainFlorisModel;
use florus::wind_data::TimeSeries;

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
