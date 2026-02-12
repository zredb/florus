//! FLORIS Configuration Structures
//!
//! This module contains configuration structures that correspond to the
//! FLORIS YAML configuration file format (gch.yaml).

use crate::types::Float;
use serde::{Deserialize, Deserializer, Serialize};
use serde_yaml::Value;
use std::collections::HashMap;

/// Top-level FLORIS configuration structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlorisConfig {
    /// Name of the configuration (optional, for reference)
    #[serde(default)]
    pub name: String,

    /// Description of the configuration (optional, for reference)
    #[serde(default)]
    pub description: Option<String>,

    /// FLORIS version (optional, for reference)
    #[serde(default, rename = "floris_version")]
    pub floris_version: String,

    /// Logging configuration
    #[serde(default)]
    pub logging: LoggingConfig,

    /// Solver configuration
    pub solver: SolverConfig,

    /// Farm configuration
    pub farm: FarmConfig,

    /// Flow field configuration
    pub flow_field: FlowFieldConfig,

    /// Wake model configuration
    #[serde(default)]
    pub wake: WakeConfig,

    /// Turbine library path (optional)
    #[serde(default, rename = "turbine_library")]
    pub turbine_library: String,
}

/// Logging configuration
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct LoggingConfig {
    #[serde(default)]
    pub console: Option<ConsoleLoggingConfig>,
    #[serde(default)]
    pub file: Option<FileLoggingConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsoleLoggingConfig {
    #[serde(default = "default_true")]
    pub enable: bool,
    #[serde(default = "default_level")]
    pub level: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileLoggingConfig {
    #[serde(default)]
    pub enable: bool,
    #[serde(default = "default_level")]
    pub level: String,
}

fn default_true() -> bool {
    true
}

fn default_level() -> String {
    "WARNING".to_string()
}

/// Solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverConfig {
    /// Solver type: "turbine_grid" or "turbine_cubature_grid"
    #[serde(rename = "type")]
    pub solver_type: SolverType,

    /// Number of grid points per turbine
    #[serde(default, rename = "turbine_grid_points")]
    pub turbine_grid_points: usize,
}
impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            solver_type: SolverType::TurbineGrid,
            turbine_grid_points: 3,
        }
    }
}

impl SolverConfig {
    pub fn new(solver_type: SolverType, turbine_grid_points: usize) -> Self {
        Self {
            solver_type,
            turbine_grid_points,
        }
    }
    pub fn is_turbine_grid(&self) -> bool {
        matches!(self.solver_type, SolverType::TurbineGrid)
    }

    pub fn is_turbine_cubature_grid(&self) -> bool {
        matches!(self.solver_type, SolverType::TurbineCubatureGrid)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum SolverType {
    TurbineGrid,
    TurbineCubatureGrid,
}

/// Farm configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FarmConfig {
    /// X-coordinates for turbine locations
    #[serde(rename = "layout_x")]
    pub layout_x: Vec<Float>,

    /// Y-coordinates for turbine locations
    #[serde(rename = "layout_y")]
    pub layout_y: Vec<Float>,

    /// Turbine types (can be single type for all turbines or list matching layout)
    #[serde(rename = "turbine_type")]
    pub turbine_type: Vec<String>,
}

/// Flow field configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowFieldConfig {
    /// Wind speeds at reference height (m/s)
    #[serde(rename = "wind_speeds")]
    pub wind_speeds: Vec<Float>,

    /// Wind directions (degrees, 0=North, 90=East)
    #[serde(rename = "wind_directions")]
    pub wind_directions: Vec<Float>,

    /// Turbulence intensities
    #[serde(rename = "turbulence_intensities")]
    pub turbulence_intensities: Vec<Float>,

    /// Air density (kg/m^3)
    #[serde(rename = "air_density")]
    pub air_density: Float,

    /// Wind shear exponent
    #[serde(rename = "wind_shear")]
    pub wind_shear: Float,

    /// Wind veer (degrees)
    #[serde(rename = "wind_veer")]
    pub wind_veer: Float,

    /// Reference wind height (m). Use -1 to use hub height.
    #[serde(rename = "reference_wind_height")]
    pub reference_wind_height: Float,

    /// Multi-dimensional conditions for Cp/Ct interpolation
    #[serde(default, rename = "multidim_conditions")]
    pub multidim_conditions: Option<MultiDimConditions>,
}

/// Multi-dimensional conditions for advanced power/thrust calculations
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct MultiDimConditions {
    /// Wave period
    #[serde(rename = "Tp")]
    pub tp: Float,
    /// Significant wave height
    #[serde(rename = "Hs")]
    pub hs: Float,
}

/// Wake model configuration
#[derive(Debug, Clone, Serialize)]
pub struct WakeConfig {
    /// Wake model selection
    pub model_strings: WakeModelStringsConfig,

    /// Enable secondary steering effects
    pub enable_secondary_steering: bool,

    /// Enable yaw added recovery effects
    pub enable_yaw_added_recovery: bool,

    /// Enable active wake mixing (Empirical Gaussian model)
    pub enable_active_wake_mixing: bool,

    /// Enable transverse velocities across turbine rotors
    pub enable_transverse_velocities: bool,

    /// Enable wake mixing (general)
    pub enable_wake_mixing: bool,

    /// Use parallel calculation
    pub use_parallel_calc: bool,
}

// 辅助结构用于反序列化原始 YAML
#[derive(Debug, Deserialize)]
struct WakeConfigRaw {
    #[serde(default)]
    model_strings: ModelStringsRaw,
    #[serde(default)]
    wake_velocity_parameters: HashMap<String, Value>,
    #[serde(default)]
    wake_deflection_parameters: HashMap<String, Value>,
    #[serde(default)]
    wake_turbulence_parameters: HashMap<String, Value>,
    #[serde(default)]
    enable_secondary_steering: bool,
    #[serde(default)]
    enable_yaw_added_recovery: bool,
    #[serde(default)]
    enable_active_wake_mixing: bool,
    #[serde(default)]
    enable_transverse_velocities: bool,
    #[serde(default)]
    enable_wake_mixing: bool,
    #[serde(default)]
    use_parallel_calc: bool,
}

#[derive(Debug, Deserialize, Default)]
struct ModelStringsRaw {
    #[serde(default)]
    velocity_model: String,
    #[serde(default)]
    deflection_model: String,
    #[serde(default)]
    turbulence_model: String,
    #[serde(default)]
    combination_model: String,
}

impl<'de> Deserialize<'de> for WakeConfig {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let raw = WakeConfigRaw::deserialize(deserializer)?;

        // 根据 model_strings 中的名称查找参数并构建枚举
        let velocity_model = parse_velocity_model(
            &raw.model_strings.velocity_model,
            &raw.wake_velocity_parameters,
        )
        .map_err(serde::de::Error::custom)?;

        let deflection_model = parse_deflection_model(
            &raw.model_strings.deflection_model,
            &raw.wake_deflection_parameters,
        )
        .map_err(serde::de::Error::custom)?;

        let turbulence_model = parse_turbulence_model(
            &raw.model_strings.turbulence_model,
            &raw.wake_turbulence_parameters,
        )
        .map_err(serde::de::Error::custom)?;

        let combination_model = parse_combination_model(&raw.model_strings.combination_model)
            .map_err(serde::de::Error::custom)?;

        Ok(WakeConfig {
            model_strings: WakeModelStringsConfig {
                velocity_model,
                deflection_model,
                turbulence_model,
                combination_model,
            },
            enable_secondary_steering: raw.enable_secondary_steering,
            enable_yaw_added_recovery: raw.enable_yaw_added_recovery,
            enable_active_wake_mixing: raw.enable_active_wake_mixing,
            enable_transverse_velocities: raw.enable_transverse_velocities,
            enable_wake_mixing: raw.enable_wake_mixing,
            use_parallel_calc: raw.use_parallel_calc,
        })
    }
}

impl Default for WakeConfig {
    fn default() -> Self {
        Self {
            model_strings: WakeModelStringsConfig::default(),
            enable_secondary_steering: false,
            enable_yaw_added_recovery: false,
            enable_active_wake_mixing: false,
            enable_transverse_velocities: false,
            enable_wake_mixing: false,
            use_parallel_calc: false,
        }
    }
}

/// Wake model string identifiers
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WakeModelStringsConfig {
    /// Wake velocity deficit model: "gauss", "jensen", "turbopark", etc.
    #[serde(rename = "velocity_model")]
    pub velocity_model: VelocityModelConfig,

    /// Wake deflection model: "gauss", "jimenez", "empirical_gauss", etc.
    #[serde(rename = "deflection_model")]
    pub deflection_model: DeflectionModelConfig,

    /// Wake turbulence model: "crespo_hernandez", "none", etc.
    #[serde(rename = "turbulence_model")]
    pub turbulence_model: TurbulenceModelConfig,

    /// Wake combination model: "fls", "sosfs", "max", etc.
    #[serde(rename = "combination_model")]
    pub combination_model: CombinationModelConfig,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum VelocityModelConfig {
    Gauss {
        alpha: Float,
        beta: Float,
        ka: Float,
        kb: Float,
    },
    Jensen {
        we: Float,
    },
    Turbopark {
        kstar: Float,
        cstar: Float,
    },
    TurboparkGauss {
        kstar: Float,
        cstar: Float,
    },
    CC {
        a_s: Float,
        b_s: Float,
        c_s1: Float,
        c_s2: Float,
        a_f: Float,
        b_f: Float,
        c_f: Float,
        alpha_mod: Float,
    },
    CumulativeGaussCurl {
        a_s: Float,
        b_s: Float,
        c_s1: Float,
        c_s2: Float,
        a_f: Float,
        b_f: Float,
        c_f: Float,
        alpha_mod: Float,
    },
    EmpiricalGauss {
        ad: Float,
        bd: Float,
        alpha: Float,
        beta: Float,
    },
}

impl Default for VelocityModelConfig {
    fn default() -> Self {
        // 返回枚举名::变体名
        VelocityModelConfig::Gauss {
            alpha: 0.38,
            beta: 0.004,
            ka: 0.38,
            kb: 0.004,
        }
    }
}

impl VelocityModelConfig {
    pub fn model_name(&self) -> String {
        match self {
            VelocityModelConfig::Gauss { .. } => "gauss".to_string(),
            VelocityModelConfig::Jensen { .. } => "jensen".to_string(),
            VelocityModelConfig::Turbopark { .. } => "turbopark".to_string(),
            VelocityModelConfig::TurboparkGauss { .. } => "turbopark_gauss".to_string(),
            VelocityModelConfig::CC { .. } => "cc".to_string(),
            VelocityModelConfig::CumulativeGaussCurl { .. } => "cumulative_gauss_curl".to_string(),
            VelocityModelConfig::EmpiricalGauss { .. } => "empirical_gauss".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum DeflectionModelConfig {
    Gauss {
        alpha: Float,
        beta: Float,
        ad: Float,
        bd: Float,
        dm: Float,
        ka: Float,
        kb: Float,
    },
    Jimenez {
        ad: Float,
        bd: Float,
        kd: Float,
    },
    EmpiricalGauss {
        ad: Float,
        bd: Float,
        kd: Float,
    },
}

impl Default for DeflectionModelConfig {
    fn default() -> Self {
        Self::Gauss {
            alpha: 0.58,
            beta: 0.077,
            ad: 0.0,
            bd: 0.0,
            dm: 1.0,
            ka: 0.38,
            kb: 0.004,
        }
    }
}

impl DeflectionModelConfig {
    pub fn model_name(&self) -> String {
        match self {
            DeflectionModelConfig::Gauss { .. } => "gauss".to_string(),
            DeflectionModelConfig::Jimenez { .. } => "jimenez".to_string(),
            DeflectionModelConfig::EmpiricalGauss { .. } => "empirical_gauss".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum TurbulenceModelConfig {
    CrespoHernandez {
        initial: Float,
        constant: Float,
        ai: Float,
        downstream: Float,
    },
    WakeInducedMixing {
        constant: Float,
        ai: Float,
    },
    None,
}

impl Default for TurbulenceModelConfig {
    fn default() -> Self {
        Self::CrespoHernandez {
            initial: 0.0,
            constant: 0.0,
            ai: 0.0,
            downstream: 0.0,
        }
    }
}

impl TurbulenceModelConfig {
    pub fn model_name(&self) -> String {
        match self {
            TurbulenceModelConfig::CrespoHernandez { .. } => "crespo_hernandez".to_string(),
            TurbulenceModelConfig::WakeInducedMixing { .. } => "wake_induced_mixing".to_string(),
            TurbulenceModelConfig::None => "none".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum CombinationModelConfig {
    FLS,
    SOSFS,
    Max,
}

impl Default for CombinationModelConfig {
    fn default() -> Self {
        Self::SOSFS
    }
}

impl CombinationModelConfig {
    pub fn model_name(&self) -> String {
        match self {
            CombinationModelConfig::FLS => "fls".to_string(),
            CombinationModelConfig::SOSFS => "sosfs".to_string(),
            CombinationModelConfig::Max => "max".to_string(),
        }
    }
}

impl Default for WakeModelStringsConfig {
    fn default() -> Self {
        Self {
            velocity_model: VelocityModelConfig::default(),
            deflection_model: DeflectionModelConfig::default(),
            combination_model: CombinationModelConfig::default(),
            turbulence_model: TurbulenceModelConfig::default(),
        }
    }
}

// 辅助函数：根据模型名称和参数解析 VelocityModelConfig
fn parse_velocity_model(
    model_name: &str,
    params_map: &HashMap<String, Value>,
) -> Result<VelocityModelConfig, String> {
    if model_name.is_empty() {
        return Ok(VelocityModelConfig::default());
    }

    let model_key = model_name.to_lowercase();
    let params = params_map
        .get(&model_key)
        .ok_or_else(|| format!("Missing parameters for velocity model: {}", model_name))?;

    match model_key.as_str() {
        "gauss" => {
            let alpha = get_f64(params, "alpha").unwrap_or(0.38);
            let beta = get_f64(params, "beta").unwrap_or(0.004);
            let ka = get_f64(params, "ka").unwrap_or(0.38);
            let kb = get_f64(params, "kb").unwrap_or(0.004);
            Ok(VelocityModelConfig::Gauss {
                alpha,
                beta,
                ka,
                kb,
            })
        }
        "jensen" => {
            let we = get_f64(params, "we").unwrap_or(0.05);
            Ok(VelocityModelConfig::Jensen { we })
        }
        "turbopark" => {
            let kstar = get_f64(params, "kstar").unwrap_or(0.05);
            let cstar = get_f64(params, "cstar").unwrap_or(1.5);
            Ok(VelocityModelConfig::Turbopark { kstar, cstar })
        }
        "turbopark_gauss" => {
            let kstar = get_f64(params, "kstar").unwrap_or(0.05);
            let cstar = get_f64(params, "cstar").unwrap_or(1.5);
            Ok(VelocityModelConfig::TurboparkGauss { kstar, cstar })
        }
        "cc" | "cumulative_gauss_curl" => {
            let a_s = get_f64(params, "a_s").unwrap_or(0.179);
            let b_s = get_f64(params, "b_s").unwrap_or(0.012);
            let c_s1 = get_f64(params, "c_s1").unwrap_or(0.056);
            let c_s2 = get_f64(params, "c_s2").unwrap_or(0.133);
            let a_f = get_f64(params, "a_f").unwrap_or(3.11);
            let b_f = get_f64(params, "b_f").unwrap_or(-0.68);
            let c_f = get_f64(params, "c_f").unwrap_or(2.41);
            let alpha_mod = get_f64(params, "alpha_mod").unwrap_or(1.0);
            Ok(VelocityModelConfig::CC {
                a_s,
                b_s,
                c_s1,
                c_s2,
                a_f,
                b_f,
                c_f,
                alpha_mod,
            })
        }
        "empirical_gauss" => {
            let ad = get_f64(params, "ad").unwrap_or(0.0);
            let bd = get_f64(params, "bd").unwrap_or(0.0);
            let alpha = get_f64(params, "alpha").unwrap_or(0.38);
            let beta = get_f64(params, "beta").unwrap_or(0.004);
            Ok(VelocityModelConfig::EmpiricalGauss {
                ad,
                bd,
                alpha,
                beta,
            })
        }
        _ => Err(format!("Unknown velocity model: {}", model_name)),
    }
}

fn parse_deflection_model(
    model_name: &str,
    params_map: &HashMap<String, Value>,
) -> Result<DeflectionModelConfig, String> {
    if model_name.is_empty() {
        return Ok(DeflectionModelConfig::default());
    }

    let model_key = model_name.to_lowercase();
    let params = params_map
        .get(&model_key)
        .ok_or_else(|| format!("Missing parameters for deflection model: {}", model_name))?;

    match model_key.as_str() {
        "gauss" => {
            let alpha = get_f64(params, "alpha").unwrap_or(0.58);
            let beta = get_f64(params, "beta").unwrap_or(0.077);
            let ad = get_f64(params, "ad").unwrap_or(0.0);
            let bd = get_f64(params, "bd").unwrap_or(0.0);
            let dm = get_f64(params, "dm").unwrap_or(1.0);
            let ka = get_f64(params, "ka").unwrap_or(0.38);
            let kb = get_f64(params, "kb").unwrap_or(0.004);
            Ok(DeflectionModelConfig::Gauss {
                alpha,
                beta,
                ad,
                bd,
                dm,
                ka,
                kb,
            })
        }
        "jimenez" => {
            let ad = get_f64(params, "ad").unwrap_or(0.0);
            let bd = get_f64(params, "bd").unwrap_or(0.0);
            let kd = get_f64(params, "kd").unwrap_or(0.05);
            Ok(DeflectionModelConfig::Jimenez { ad, bd, kd })
        }
        "empirical_gauss" => {
            let ad = get_f64(params, "ad").unwrap_or(0.0);
            let bd = get_f64(params, "bd").unwrap_or(0.0);
            let kd = get_f64(params, "kd").unwrap_or(0.05);
            Ok(DeflectionModelConfig::EmpiricalGauss { ad, bd, kd })
        }
        _ => Err(format!("Unknown deflection model: {}", model_name)),
    }
}

fn parse_turbulence_model(
    model_name: &str,
    params_map: &HashMap<String, Value>,
) -> Result<TurbulenceModelConfig, String> {
    if model_name.is_empty() {
        // 当没有指定模型名称时，使用默认值
        return Ok(TurbulenceModelConfig::default());
    }

    if model_name.to_lowercase() == "none" {
        return Ok(TurbulenceModelConfig::None);
    }

    let model_key = model_name.to_lowercase();
    let params = params_map
        .get(&model_key)
        .ok_or_else(|| format!("Missing parameters for turbulence model: {}", model_name))?;

    match model_key.as_str() {
        "crespo_hernandez" => {
            let initial = get_f64(params, "initial").unwrap_or(0.0);
            let constant = get_f64(params, "constant").unwrap_or(0.0);
            let ai = get_f64(params, "ai").unwrap_or(0.0);
            let downstream = get_f64(params, "downstream").unwrap_or(0.0);
            Ok(TurbulenceModelConfig::CrespoHernandez {
                initial,
                constant,
                ai,
                downstream,
            })
        }
        "wake_induced_mixing" => {
            let constant = get_f64(params, "constant").unwrap_or(0.0);
            let ai = get_f64(params, "ai").unwrap_or(0.0);
            Ok(TurbulenceModelConfig::WakeInducedMixing { constant, ai })
        }
        "none" => Ok(TurbulenceModelConfig::None),
        _ => Err(format!("Unknown turbulence model: {}", model_name)),
    }
}

fn parse_combination_model(model_name: &str) -> Result<CombinationModelConfig, String> {
    if model_name.is_empty() {
        return Ok(CombinationModelConfig::default());
    }

    match model_name.to_lowercase().as_str() {
        "fls" => Ok(CombinationModelConfig::FLS),
        "sosfs" => Ok(CombinationModelConfig::SOSFS),
        "max" => Ok(CombinationModelConfig::Max),
        _ => Err(format!("Unknown combination model: {}", model_name)),
    }
}

// 辅助函数：从 Value 中提取 f64
fn get_f64(value: &Value, key: &str) -> Option<Float> {
    value.get(key).and_then(|v| match v {
        Value::Number(n) => n.as_f64().map(|f| f as Float),
        Value::String(s) => s.parse::<Float>().ok(),
        _ => None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gch_yaml() {
        // This tests that we can parse a minimal FLORIS config
        let yaml = r#"
name: Test Config
solver:
  type: turbine_grid
  turbine_grid_points: 3
farm:
  layout_x: [0.0, 630.0]
  layout_y: [0.0, 0.0]
  turbine_type: [nrel_5MW]
flow_field:
  wind_speeds: [8.0]
  wind_directions: [270.0]
  turbulence_intensities: [0.06]
  air_density: 1.225
  wind_shear: 0.12
  wind_veer: 0.0
  reference_wind_height: 90.0
wake:
  model_strings:
    velocity_model: gauss
    deflection_model: gauss
    combination_model: fls
    turbulence_model: crespo_hernandez
  wake_deflection_parameters:
    gauss:
      alpha: 0.58
      beta: 0.077
  wake_velocity_parameters:
    gauss:
      ka: 0.38
      kb: 0.004
  wake_turbulence_parameters:
    crespo_hernandez:
      initial: 0.1
      constant: 0.5
"#;
        let config: FlorisConfig = serde_yaml::from_str(yaml).unwrap();
        assert_eq!(config.name, "Test Config");
        assert_eq!(config.solver.solver_type, SolverType::TurbineGrid);
        assert_eq!(config.solver.turbine_grid_points, 3);
        assert_eq!(config.farm.layout_x.len(), 2);
        assert_eq!(config.flow_field.wind_speeds, vec![8.0]);

        // 现在 velocity_model 应该包含从 wake_velocity_parameters 中解析的参数
        match config.wake.model_strings.velocity_model {
            VelocityModelConfig::Gauss { ka, kb, .. } => {
                assert_eq!(ka, 0.38);
                assert_eq!(kb, 0.004);
            }
            _ => panic!("Expected Gauss velocity model"),
        }
    }

    #[test]
    fn test_parse_nested_wake_params() {
        let yaml = r#"
solver:
  type: turbine_grid
farm:
  layout_x: [0.0]
  layout_y: [0.0]
  turbine_type: [nrel_5MW]
flow_field:
  wind_speeds: [8.0]
  wind_directions: [270.0]
  turbulence_intensities: [0.06]
  air_density: 1.225
  wind_shear: 0.12
  wind_veer: 0.0
  reference_wind_height: 90.0
wake:
  model_strings:
    velocity_model: gauss
    deflection_model: gauss
    combination_model: fls
    turbulence_model: crespo_hernandez
  wake_deflection_parameters:
    gauss:
      ad: 0.0
      bd: 0.0
      alpha: 0.58
      beta: 0.077
      dm: 1.0
      ka: 0.38
      kb: 0.004
    jimenez:
      ad: 0.0
      bd: 0.0
      kd: 0.05
  wake_velocity_parameters:
    gauss:
      alpha: 0.58
      beta: 0.077
      ka: 0.38
      kb: 0.004
    cc:
      a_s: 0.179
      b_s: 0.012
      c_s1: 0.056
      c_s2: 0.133
      a_f: 3.11
      b_f: -0.68
      c_f: 2.41
      alpha_mod: 1.0
  wake_turbulence_parameters:
    crespo_hernandez:
      initial: 0.1
      constant: 0.5
      ai: 0.8
      downstream: -0.32
"#;
        let config: FlorisConfig = serde_yaml::from_str(yaml).unwrap();
        let wake = &config.wake;

        // Check gauss deflection params
        match &wake.model_strings.deflection_model {
            DeflectionModelConfig::Gauss { ad, alpha, .. } => {
                assert_eq!(*ad, 0.0);
                assert_eq!(*alpha, 0.58);
            }
            _ => panic!("Expected Gauss deflection model"),
        }

        // Check gauss velocity params
        match &wake.model_strings.velocity_model {
            VelocityModelConfig::Gauss { ka, kb, .. } => {
                assert_eq!(*ka, 0.38);
                assert_eq!(*kb, 0.004);
            }
            _ => panic!("Expected Gauss velocity model"),
        }

        // Check crespo hernandez turbulence params
        match &wake.model_strings.turbulence_model {
            TurbulenceModelConfig::CrespoHernandez { initial, ai, .. } => {
                assert_eq!(*initial, 0.1);
                assert_eq!(*ai, 0.8);
            }
            _ => panic!("Expected CrespoHernandez turbulence model"),
        }
    }

    #[test]
    fn test_default_wake_config() {
        let config: WakeConfig = serde_yaml::from_str("").unwrap();
        assert_eq!(
            config.model_strings.velocity_model,
            VelocityModelConfig::Gauss {
                alpha: 0.38,
                beta: 0.004,
                ka: 0.38,
                kb: 0.004,
            }
        );
        assert_eq!(
            config.model_strings.deflection_model,
            DeflectionModelConfig::Gauss {
                alpha: 0.58,
                beta: 0.077,
                ad: 0.0,
                bd: 0.0,
                dm: 1.0,
                ka: 0.38,
                kb: 0.004,
            }
        );
        assert_eq!(
            config.model_strings.combination_model,
            CombinationModelConfig::SOSFS
        );
        assert_eq!(
            config.model_strings.turbulence_model,
            TurbulenceModelConfig::CrespoHernandez {
                initial: 0.0,
                constant: 0.0,
                ai: 0.0,
                downstream: 0.0,
            }
        );
        assert!(!config.enable_secondary_steering);
    }
}
