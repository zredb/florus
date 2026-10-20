//! Cross-validation of the floating-turbine wake behaviour against the
//! Python FLORIS reference values.
//!
//! Reference values come from `floris/examples/examples_floating/`:
//! `002_floating_vs_fixedbottom_farm.py` (empirical Gaussian wake model, a
//! 20-turbine farm with a fixed-bottom and a floating turbine).

use florus::floris_model::FlorisModel;
use florus::types::Array1;

/// The 5 x 4 grid of the example farm
fn farm_layout() -> (Array1, Array1) {
    let mut layout_x: Vec<f64> = Vec::new();
    let mut layout_y: Vec<f64> = Vec::new();
    for j in 0..4 {
        for i in 0..5 {
            layout_x.push(i as f64 * 630.0);
            layout_y.push(j as f64 * 630.0);
        }
    }
    (Array1::from_vec(layout_x), Array1::from_vec(layout_y))
}

fn run(input_file: &str) -> Vec<f64> {
    let mut fmodel = FlorisModel::from_file(input_file).unwrap();
    let (layout_x, layout_y) = farm_layout();
    fmodel.set_layout(&layout_x, &layout_y).unwrap();
    fmodel
        .set_wind_conditions(
            Array1::from_vec(vec![10.0]),
            Array1::from_vec(vec![270.0]),
            Array1::from_vec(vec![0.06]),
        )
        .unwrap();
    fmodel.run().unwrap();
    let powers = fmodel.get_turbine_powers().unwrap();
    (0..powers.shape()[1])
        .map(|t| powers[[0, t]] / 1000.0)
        .collect()
}

#[test]
fn test_floating_vs_fixed_bottom_matches_python() {
    // Python: examples_floating/002_floating_vs_fixedbottom_farm.py
    let fixed = run("examples/inputs/emgauss_fixed.yaml");
    let floating = run("examples/inputs/emgauss_floating.yaml");

    // Python reference (kW): the upstream turbine, the first waked turbine and
    // the last turbine of the first row
    let indices = [0usize, 1, 4];
    let expected_fixed = [3417.8, 1343.1, 1481.7];
    let expected_floating = [3377.2, 1385.7, 1501.3];
    for (k, i) in indices.iter().enumerate() {
        let i = *i;
        assert!(
            (fixed[i] - expected_fixed[k]).abs() < 0.1,
            "fixed bottom turbine {i}: {} kW (expected {})",
            fixed[i],
            expected_fixed[k]
        );
        assert!(
            (floating[i] - expected_floating[k]).abs() < 0.1,
            "floating turbine {i}: {} kW (expected {})",
            floating[i],
            expected_floating[k]
        );
    }

    // The farm powers
    let sum_fixed: f64 = fixed.iter().sum();
    let sum_floating: f64 = floating.iter().sum();
    assert!(
        (sum_fixed - 36674.6).abs() < 1.0,
        "fixed farm power {sum_fixed}"
    );
    assert!(
        (sum_floating - 36975.1).abs() < 1.0,
        "floating farm power {sum_floating}"
    );
}
