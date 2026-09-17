use serde::Deserialize;
use std::path::Path;

#[derive(Debug, Deserialize)]
pub struct WakeRegressionGolden {
    pub turbine_powers: Vec<Vec<f64>>,
    pub farm_power: Vec<f64>,
    pub n_findex: usize,
    pub n_turbines: usize,
    pub velocities_shape: Vec<usize>,
    pub config: WakeConfig,
}

#[derive(Debug, Deserialize)]
pub struct WakeConfig {
    pub velocity_model: String,
    pub deflection_model: String,
    pub combination_model: String,
    pub turbulence_model: String,
    pub yawed: bool,
}

#[derive(Debug, Deserialize)]
pub struct PowerGolden {
    pub turbine_powers: Vec<Vec<f64>>,
    pub farm_power: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct OperationModelGolden {
    pub turbine_powers: Vec<Vec<f64>>,
    pub farm_power: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct TimeSeriesGolden {
    pub timeseries: TimeSeriesData,
}

#[derive(Debug, Deserialize)]
pub struct TimeSeriesData {
    pub n_findex: usize,
    pub wind_speeds: Vec<f64>,
    pub wind_directions: Vec<f64>,
    pub turbulence_intensities: Vec<f64>,
}

#[derive(Debug, Deserialize)]
pub struct WindRoseGolden {
    pub windrose: WindRoseData,
}

#[derive(Debug, Deserialize)]
pub struct WindRoseData {
    pub n_wind_directions: usize,
    pub n_wind_speeds: usize,
    pub wind_directions: Vec<f64>,
    pub wind_speeds: Vec<f64>,
    pub freq_table_shape: Vec<usize>,
    pub freq_table_sum: f64,
}

pub fn load_wake_golden(name: &str) -> WakeRegressionGolden {
    let path = format!("tests/golden_values/wake_regression/{}.json", name);
    let data = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {}", path, e));
    serde_json::from_str(&data).unwrap_or_else(|e| panic!("Failed to parse {}: {}", path, e))
}

pub fn load_api_golden(name: &str) -> PowerGolden {
    let path = format!("tests/golden_values/floris_model_api/{}.json", name);
    let data = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {}", path, e));
    serde_json::from_str(&data).unwrap_or_else(|e| panic!("Failed to parse {}: {}", path, e))
}

pub fn load_operation_golden(name: &str) -> OperationModelGolden {
    let path = format!("tests/golden_values/operation_models/{}.json", name);
    let data = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {}: {}", path, e));
    serde_json::from_str(&data).unwrap_or_else(|e| panic!("Failed to parse {}: {}", path, e))
}

pub fn load_timeseries_golden() -> TimeSeriesGolden {
    let path = "tests/golden_values/wind_data/timeseries.json";
    let data = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {}: {}", path, e));
    serde_json::from_str(&data).unwrap_or_else(|e| panic!("Failed to parse {}: {}", path, e))
}

pub fn load_windrose_golden() -> WindRoseGolden {
    let path = "tests/golden_values/wind_data/windrose.json";
    let data = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {}: {}", path, e));
    serde_json::from_str(&data).unwrap_or_else(|e| panic!("Failed to parse {}: {}", path, e))
}

pub fn assert_relative_eq_slice(
    test: &[f64],
    baseline: &[f64],
    epsilon: f64,
    description: &str,
) {
    assert_eq!(test.len(), baseline.len(), "{}: length mismatch", description);
    for (i, (t, b)) in test.iter().zip(baseline.iter()).enumerate() {
        let diff = (t - b).abs();
        let rel_diff = if b.abs() > 1e-10 { diff / b.abs() } else { diff };
        assert!(
            diff < epsilon || rel_diff < epsilon,
            "{}[{}]: got {:.6}, expected {:.6}, diff {:.6}, rel_diff {:.6}",
            description, i, t, b, diff, rel_diff
        );
    }
}

pub fn assert_2d_close(
    test: &[Vec<f64>],
    baseline: &[Vec<f64>],
    epsilon: f64,
    description: &str,
) {
    assert_eq!(test.len(), baseline.len(), "{}: row count mismatch", description);
    for (i, (t_row, b_row)) in test.iter().zip(baseline.iter()).enumerate() {
        assert_relative_eq_slice(t_row, b_row, epsilon, &format!("{}[{}]", description, i));
    }
}
