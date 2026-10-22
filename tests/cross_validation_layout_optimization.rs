//! Cross-validation of the layout optimization example against the Python
//! FLORIS reference values.
//!
//! Reference values come from
//! `floris/examples/examples_layout_optimization/001_optimize_layout.py`.
//!
//! The wind resource of the example is drawn with numpy's random number
//! generator, so the Rust example uses `florus::numpy_random`, which
//! reproduces numpy's legacy `RandomState` exactly.

use florus::floris_model::FlorisModel;
use florus::numpy_random::NumpyRandomState;
use florus::optimization::layout_optimization::layout_optimization_base::{
    Boundary, LayoutOptimizationConfig,
};
use florus::optimization::layout_optimization::{LayoutOptimizationScipy, LayoutOptimizer};
use florus::types::{Array1, Array2};
use florus::wind_data::WindRose;

const D: f64 = 126.0;

/// The wind rose of the example: 72 wind directions, a frequency table
/// drawn with numpy's `np.random.seed(1)` and `np.random.randn(72)`, and a
/// value table that is high for the northerly and southerly winds.
fn wind_rose() -> (WindRose, Vec<f64>) {
    let wind_directions = Array1::from_iter((0..72).map(|i| i as f64 * 5.0));
    let wind_speeds = Array1::from_vec(vec![8.0]);

    let mut rng = NumpyRandomState::new(1);
    let mut frequencies = rng.standard_normal_vec(72);
    frequencies.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let frequencies: Vec<f64> = frequencies.iter().map(|v| v.abs()).collect();
    let total: f64 = frequencies.iter().sum();
    let freq_table = Array2::from_shape_fn((72, 1), |(i, _)| frequencies[i] / total);

    let raw_values: Vec<f64> = wind_directions
        .iter()
        .map(|&wd| (0.5 + 0.5 * (2.0 * wd.to_radians()).cos()).powi(10))
        .collect();
    let mean_value: f64 = raw_values.iter().sum::<f64>() / raw_values.len() as f64;
    let value_table = Array2::from_shape_fn((72, 1), |(i, _)| 25.0 * raw_values[i] / mean_value);

    let wind_rose = WindRose::new(
        wind_directions,
        wind_speeds,
        Array2::from_elem((72, 1), 0.06),
        Some(freq_table),
        Some(value_table),
        false,
        None,
        None,
    )
    .expect("wind rose");

    (wind_rose, frequencies)
}

/// A wind rose with every sixth wind direction (the optimization test uses
/// it, because it is six times cheaper to evaluate).
fn wind_rose_subset() -> (WindRose, Vec<f64>) {
    let (wind_rose, frequencies) = wind_rose();
    let index: Vec<usize> = (0..72).step_by(6).collect();

    let wind_directions = Array1::from_iter(index.iter().map(|&i| i as f64 * 5.0));
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let total: f64 = frequencies.iter().sum();
    let freq_table =
        Array2::from_shape_fn((index.len(), 1), |(i, _)| frequencies[index[i]] / total);
    let raw_values: Vec<f64> = wind_directions
        .iter()
        .map(|&wd| (0.5 + 0.5 * (2.0 * wd.to_radians()).cos()).powi(10))
        .collect();
    let mean_value: f64 = raw_values.iter().sum::<f64>() / raw_values.len() as f64;
    let value_table =
        Array2::from_shape_fn((index.len(), 1), |(i, _)| 25.0 * raw_values[i] / mean_value);

    let wind_rose = WindRose::new(
        wind_directions,
        wind_speeds,
        Array2::from_elem((index.len(), 1), 0.06),
        Some(freq_table),
        Some(value_table),
        false,
        None,
        None,
    )
    .expect("wind rose");

    (
        wind_rose,
        index.iter().map(|&i| frequencies[i] / total).collect(),
    )
}

/// The frequency table and the AEP/AVP of the initial layout must match
/// Python exactly.
#[test]
fn test_layout_optimization_wind_resource_matches_python() {
    let (wind_rose, frequencies) = wind_rose();

    // The frequencies of the wind rose (Python: the first entries of the
    // frequency table of the example)
    let expected = [
        0.0403453082,
        0.0361136713,
        0.0253148977,
        0.0244801772,
        0.0200280138,
    ];
    let total: f64 = frequencies.iter().sum();
    for (i, expected) in expected.iter().enumerate() {
        let value = frequencies[i] / total;
        assert!(
            (value - expected).abs() < 1e-9,
            "the frequency of the direction {}: {value} (expected {expected})",
            i * 5
        );
    }

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    fmodel.set_wind_data(&wind_rose).expect("wind data");
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 0.0, 6.0 * D, 6.0 * D]),
            &Array1::from_vec(vec![0.0, 4.0 * D, 0.0, 4.0 * D]),
        )
        .expect("layout");
    fmodel.run().expect("run");

    // Python: 56299.1812 MWh and 1527164.8088 dollars. The Rust values
    // differ by about 1e-6 relative (the waked turbines' powers differ by
    // 0.03 % for a few wind directions).
    let aep = fmodel.get_expected_farm_power(None, None).expect("aep") * 8760.0 / 1e6;
    let avp = fmodel
        .get_expected_farm_value(None, None, None)
        .expect("avp")
        * 8760.0
        / 1e6;
    assert!(
        (aep - 56299.1812).abs() < 0.5,
        "the AEP of the initial layout: {aep} MWh (expected 56299.1812 MWh)"
    );
    assert!(
        (avp - 1527164.8088).abs() < 100.0,
        "the AVP of the initial layout: {avp} dollars (expected 1527164.8088 dollars)"
    );
}

/// The genetic random search example
/// (`examples_layout_optimization/003_genetic_random_search.py`): the random
/// wind resource and the AEP of the initial layout must match Python.
#[test]
fn test_genetic_random_search_matches_python() {
    let (wind_rose, _) = wind_rose();

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    fmodel.set_wind_data(&wind_rose).expect("wind data");
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 0.0, 6.0 * D, 6.0 * D]),
            &Array1::from_vec(vec![0.0, 4.0 * D, 0.0, 4.0 * D]),
        )
        .expect("layout");
    fmodel.run().expect("run");

    // Python: Initial AEP = 56.3 GWh
    let aep = fmodel.get_expected_farm_power(None, None).expect("aep") * 8760.0 / 1e9;
    assert!(
        (aep - 56.3).abs() < 0.05,
        "the AEP of the initial layout: {aep} GWh (expected 56.3 GWh)"
    );
}

/// The layout optimization with heterogeneous inflow
/// (`examples_layout_optimization/002_optimize_layout_with_heterogeneity.py`):
/// the AEP of the initial layout must match Python, and the optimization
/// must improve it by about the same amount.
#[test]
fn test_layout_optimization_with_heterogeneity_matches_python() {
    use florus::heterogeneous_map::HeterogeneousInflowConfig;

    let d = 126.0;
    let size_d = 12.0;
    let wind_directions = Array1::from_vec(vec![90.0, 270.0]);
    let n_wds = wind_directions.len();

    let wind_rose = WindRose::new(
        wind_directions.clone(),
        Array1::from_vec(vec![8.0]),
        Array2::from_elem((n_wds, 1), 0.06),
        Some(Array2::from_elem((n_wds, 1), 1.0 / n_wds as f64)),
        None,
        false,
        None,
        None,
    )
    .expect("wind rose");

    let heterogeneous_config = HeterogeneousInflowConfig {
        x: Array1::from_vec(vec![0.0, size_d * d, 0.0, size_d * d]),
        y: Array1::from_vec(vec![-d, -d, d, d]),
        z: None,
        wind_speeds: None,
        wind_directions: Some(wind_directions.clone()),
        speed_multipliers: Array2::from_shape_fn((n_wds, 4), |(_, j)| [0.5, 1.0, 0.5, 1.0][j]),
    };

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    fmodel.set_wind_data(&wind_rose).expect("wind data");
    fmodel
        .set_heterogeneous_inflow_config(heterogeneous_config)
        .expect("heterogeneous inflow");
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.1, 0.3 * size_d * d, 0.6 * size_d * d]),
            &Array1::from_vec(vec![0.0, 0.0, 0.0]),
        )
        .expect("layout");
    fmodel.run().expect("run");

    // Python: 6198.3 MWh for the initial layout
    let base_aep = fmodel.get_expected_farm_power(None, None).expect("aep") * 8760.0 / 1e6;
    assert!(
        (base_aep - 6198.3).abs() < 0.5,
        "the AEP of the initial layout: {base_aep} MWh (expected 6198.3 MWh)"
    );

    let boundaries = Boundary::Polygon(vec![
        (0.0, 0.0),
        (size_d * d, 0.0),
        (size_d * d, 0.1),
        (0.0, 0.1),
        (0.0, 0.0),
    ]);
    let config = LayoutOptimizationConfig {
        boundaries: boundaries.clone(),
        min_dist: Some(2.0 * d),
        max_iterations: 40,
        ..Default::default()
    };
    let mut optimizer = LayoutOptimizationScipy::new(&fmodel, boundaries)
        .expect("optimizer")
        .with_config(config);
    let solution = optimizer.optimize().expect("optimize");

    fmodel.set_layout(&solution.x, &solution.y).expect("layout");
    fmodel.run().expect("run");
    let aep = fmodel.get_expected_farm_power(None, None).expect("aep") * 8760.0 / 1e6;

    // Python: 13930.0 MWh (an improvement of 124.7 %)
    assert!(
        (aep - 13930.0).abs() / 13930.0 < 0.05,
        "the AEP of the optimized layout: {aep} MWh (expected 13930.0 MWh)"
    );
}
/// Python example improves it by 4.5 % with scipy's SLSQP; the Rust
/// optimizer is a pattern search, which converges to a nearby optimum).
#[test]
fn test_layout_optimization_improves_the_aep() {
    let (wind_rose, _) = wind_rose_subset();

    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml").expect("model");
    fmodel.set_wind_data(&wind_rose).expect("wind data");
    fmodel
        .set_layout(
            &Array1::from_vec(vec![0.0, 0.0, 6.0 * D, 6.0 * D]),
            &Array1::from_vec(vec![0.0, 4.0 * D, 0.0, 4.0 * D]),
        )
        .expect("layout");

    // Python: 50764.446841 MWh for the initial layout
    fmodel.run().expect("run");
    let base_aep = fmodel.get_expected_farm_power(None, None).expect("aep") * 8760.0 / 1e6;
    assert!(
        (base_aep - 50764.446841).abs() < 0.5,
        "the AEP of the initial layout: {base_aep} MWh (expected 50764.446841 MWh)"
    );

    let boundaries = Boundary::Polygon(vec![
        (0.0, 0.0),
        (0.0, 1000.0),
        (1000.0, 1000.0),
        (1000.0, 0.0),
        (0.0, 0.0),
    ]);
    let config = LayoutOptimizationConfig {
        boundaries: boundaries.clone(),
        // The test runs the optimizer with a reduced number of iterations
        // (the example uses the default)
        max_iterations: 40,
        ..Default::default()
    };
    let mut optimizer = LayoutOptimizationScipy::new(&fmodel, boundaries)
        .expect("optimizer")
        .with_config(config);
    let solution = optimizer.optimize().expect("optimize");

    fmodel.set_layout(&solution.x, &solution.y).expect("layout");
    fmodel.run().expect("run");
    let aep = fmodel.get_expected_farm_power(None, None).expect("aep") * 8760.0 / 1e6;
    assert!(
        aep > base_aep * 1.02,
        "the optimized AEP is not higher than the initial one: {aep} MWh"
    );
}
