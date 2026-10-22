//! Cross-validation of the wind resource grid (WRG) against the Python
//! FLORIS reference values.
//!
//! Reference values come from
//! `floris/examples/examples_wind_resource_grid/002_set_floris_with_wrg.py`.

use florus::floris_model::FlorisModel;
use florus::types::Array1;
use florus::wind_data::WindRoseWRG;

#[test]
fn test_wrg_expected_powers_match_python() {
    // Python: examples_wind_resource_grid/002_set_floris_with_wrg.py
    let mut wind_rose_wrg = WindRoseWRG::new("examples/inputs/wrg_example.wrg", None, None, None)
        .expect("failed to read the WRG file");

    // Use 2 m/s steps and a fixed TI of 7% for all wind directions
    wind_rose_wrg.set_wind_speeds(Array1::from_iter((0..13).map(|i| i as f64 * 2.0)));
    wind_rose_wrg.set_ti_table(0.07);

    // A line of turbines along the diagonals of the WRG grid
    let layout_x = Array1::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = Array1::from_vec(vec![0.0, 1000.0, 2000.0]);
    wind_rose_wrg
        .set_layout(layout_x.clone(), layout_y.clone())
        .unwrap();

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    fmodel.set_layout(&layout_x, &layout_y).unwrap();
    fmodel.set_wind_data(&wind_rose_wrg).unwrap();
    fmodel.run().unwrap();

    // Python reference (W)
    let expected = [2280015.29955901, 2272352.87711173, 2324102.91346193];
    let powers = fmodel.get_expected_turbine_powers(None).unwrap();
    for t in 0..3 {
        assert!(
            (powers[t] - expected[t]).abs() < 1.0,
            "turbine {t}: {} W (expected {})",
            powers[t],
            expected[t]
        );
    }

    // The farm power
    let farm_power = fmodel.get_expected_farm_power(None, None).unwrap();
    assert!(
        (farm_power - 6876471.090132672).abs() < 10.0,
        "farm power {farm_power} W"
    );
}

/// The wind rose of the wind resource grid with a 2 degree wind direction
/// step is upsampled from the 12 sectors of the WRG file (Python's
/// `get_wind_rose_at_point`).
#[test]
fn test_wrg_upsampled_wind_rose() {
    let mut wind_rose_wrg = WindRoseWRG::new(
        "examples/inputs/wrg_example.wrg",
        Some(2.0),
        Some(Array1::from_iter((0..21).map(|i| i as f64))),
        None,
    )
    .expect("failed to read the WRG file");
    wind_rose_wrg.set_ti_table(0.06);

    // The wind rose of the grid point at the origin
    let wind_rose = wind_rose_wrg.get_wind_rose_at_point(0.0, 0.0, None, None, None);

    // The wind directions are upsampled from the 12 sectors of the WRG file
    assert_eq!(wind_rose.wind_directions.len(), 180);
    assert_eq!(wind_rose.wind_speeds.len(), 21);
    assert!(
        wind_rose.freq_table.iter().all(|f| f.is_finite()),
        "the frequency table contains invalid values"
    );
    assert!(
        (wind_rose.freq_table.sum() - 1.0).abs() < 1e-9,
        "the frequencies sum to {}",
        wind_rose.freq_table.sum()
    );
}

/// The heterogeneous wind rose of the wind resource grid (Python's
/// `get_heterogeneous_wind_rose`): the speed multipliers of the grid are the
/// ratios of the expected power of the turbine to the expected power at the
/// reference location.
#[test]
fn test_wrg_heterogeneous_wind_rose() {
    let mut wind_rose_wrg = WindRoseWRG::new(
        "examples/inputs/wrg_example.wrg",
        Some(2.0),
        Some(Array1::from_iter((0..21).map(|i| i as f64))),
        None,
    )
    .expect("failed to read the WRG file");
    wind_rose_wrg.set_ti_table(0.06);

    let fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").unwrap();
    let wind_rose_het = wind_rose_wrg
        .get_heterogeneous_wind_rose(&fmodel, None, Some(0.0), Some(0.0), 9.0)
        .expect("heterogeneous wind rose");

    // The wind rose has the representative wind speed only
    assert_eq!(wind_rose_het.wind_speeds.len(), 1);
    assert!((wind_rose_het.wind_speeds[0] - 9.0).abs() < 1e-9);

    // The heterogeneous map describes the speed ups of the grid
    let map = wind_rose_het
        .heterogeneous_map
        .as_ref()
        .expect("the heterogeneous map");
    assert_eq!(map.x.len(), 6);
    assert_eq!(map.y.len(), 6);
    assert_eq!(map.speed_multipliers.shape(), &[12, 6]);
    for multiplier in map.speed_multipliers.iter() {
        assert!(
            (1.0..1.1).contains(multiplier),
            "the speed multiplier {multiplier} is out of range"
        );
    }

    // The wind rose with the heterogeneous map gives a non-zero power
    let mut fmodel = fmodel.clone();
    fmodel
        .set_layout(
            &Array1::from_vec(vec![500.0, 900.0]),
            &Array1::from_vec(vec![900.0, 1900.0]),
        )
        .unwrap();
    fmodel.set_wind_data(&wind_rose_het).unwrap();
    fmodel.run().unwrap();
    let farm_power = fmodel.get_expected_farm_power(None, None).unwrap();
    assert!(
        farm_power > 1.0e6,
        "the expected farm power is {farm_power} W"
    );
}
