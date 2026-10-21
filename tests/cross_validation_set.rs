//! Cross-validation of the "setting the model parameters" example against
//! the Python FLORIS reference values.
//!
//! Reference values come from `floris/examples/004_set.py`.

use florus::floris_model::FlorisModel;
use florus::types::{Array1, Array2};
use florus::wind_data::{TimeSeries, WindRose};

#[test]
fn test_set_example_matches_python() -> anyhow::Result<()> {
    // Python: examples/004_set.py
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;

    // The wind conditions, first with arrays, then with a TimeSeries and a
    // WindRose object
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![8.0, 9.0, 10.0]),
            Array1::from_vec(vec![270.0, 270.0, 270.0]),
            Array1::from_vec(vec![0.06, 0.06, 0.06]),
        )
        .unwrap();

    let time_series = TimeSeries::new(
        Array1::from_vec(vec![270.0, 270.0, 270.0]),
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array1::from_vec(vec![0.06, 0.06, 0.06]),
    )
    .unwrap();
    fmodel.set_wind_data(&time_series).unwrap();

    let wind_rose = WindRose::new(
        Array1::from_vec(vec![270.0]),
        Array1::from_vec(vec![8.0, 9.0, 10.0]),
        Array2::from_elem((1, 3), 0.06),
        None,
        None,
        false,
        None,
        None,
    )
    .unwrap();
    fmodel.set_wind_data(&wind_rose).unwrap();

    // The atmospheric conditions and the layout
    fmodel.set_wind_shear(0.2).unwrap();
    fmodel.set_air_density(1.1).unwrap();
    fmodel.set_reference_wind_height(92.0).unwrap();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 500.0]),
            &Array1::from_vec(vec![0.0, 0.0]),
        )
        .unwrap();

    // The control settings
    fmodel
        .set_yaw_angles(Array2::from_shape_vec(
            (3, 2),
            vec![0.0, 0.0, 25.0, 0.0, 0.0, 0.0],
        )?)
        .unwrap();
    fmodel.reset_operation();
    fmodel.set_operation_model("mixed").unwrap();

    // Shut down the front turbine for the first two wind conditions
    let mut disable_turbines = ndarray::Array2::<bool>::from_elem((3, 2), false);
    disable_turbines[[0, 0]] = true;
    disable_turbines[[1, 0]] = true;
    fmodel
        .set_operation(None, None, None, None, None, Some(disable_turbines))
        .unwrap();

    // ... and derate it afterwards (the setpoints take precedence)
    let rated_power = 5.0e6;
    let power_setpoints = Array2::from_shape_vec(
        (3, 2),
        vec![
            rated_power * 0.3,
            rated_power,
            rated_power * 0.3,
            rated_power,
            rated_power,
            rated_power,
        ],
    )
    .unwrap();
    fmodel
        .set_operation(None, Some(power_setpoints), None, None, None, None)
        .unwrap();

    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();

    // Python reference (kW)
    let expected = [[1500.0, 339.4], [1500.0, 850.6], [3051.1, 657.8]];
    for fi in 0..3 {
        for t in 0..2 {
            assert!(
                (powers[[fi, t]] / 1000.0 - expected[fi][t]).abs() < 0.1,
                "condition {fi}, turbine {t}: {} kW (expected {})",
                powers[[fi, t]] / 1000.0,
                expected[fi][t]
            );
        }
    }

    Ok(())
}
