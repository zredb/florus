use ndarray::{ArrayD, IxDyn};
use std::collections::HashMap;

// 风机类型的基本特征
pub trait TurbineModel {
    fn rotor_diameter(&self) -> f64;
    fn hub_height(&self) -> f64;
    fn power_coefficient(&self, wind_speed: f64, conditions: &OperatingConditions) -> f64;
    fn thrust_coefficient(&self, wind_speed: f64, conditions: &OperatingConditions) -> f64;
}

// 运行条件结构体
#[derive(Debug, Clone)]
pub struct OperatingConditions {
    pub wind_speed: f64,
    pub turbulence_intensity: Option<f64>,
    pub yaw_angle: Option<f64>,
    pub air_density: Option<f64>,
    pub blade_pitch: Option<f64>,
    pub rotor_speed: Option<f64>,
    pub wind_direction: Option<f64>,
    pub shear_exponent: Option<f64>,
    pub inflow_angle: Option<f64>,
}

// 基础风机实现
pub struct TurbineType {
    pub name: String,
    pub rotor_diameter: f64,
    pub hub_height: f64,
    pub tsr: Option<f64>,
    pub rated_power: Option<f64>,
    pub cut_in_wind_speed: Option<f64>,
    pub cut_out_wind_speed: Option<f64>,
    pub rated_wind_speed: Option<f64>,
    pub max_yaw_angle: Option<f64>,
    pub yaw_rate_limit: Option<f64>,
    pub generator_efficiency: Option<f64>,
    pub drive_train_efficiency: Option<f64>,
    pub additional_losses: Option<f64>,
    pub ti_ref: Option<f64>,
    pub power_thrust_table: PowerTrustParams,
    pub floating_tilt_table: Option<FloatingTiltTable>,
    pub correct_cp_ct_for_tilt: bool,
}
pub struct FloatingTiltTable {
    pub wind_speed: FloatArray,
    pub tilt_angle: FloatArray,
}

pub struct PowerTrustParams {
    pub cp_ct_table: CpCtTable,

    pub ref_air_density: Option<f64>,
    pub ref_tilt: Option<f64>,
    pub cosine_loss_exponent_yaw: Option<f64>,
    pub cosine_loss_exponent_tilt: Option<f64>,
    pub generator_efficiency: Float, // 发电机效率
    pub TI_ref: Float,               // 参考湍流强度
    // Peak shaving parameters
    pub peak_shaving_fraction: Option<f64>,
    pub peak_shaving_TI_threshold: Option<f64>,
    // Power thrust data file for multi-dimensional CP/CT
    pub power_thrust_data_file: Option<String>,
    pub helix: Option<Helix>,
    pub controller_dependent_turbine_parameters: Option<ControllerDependentTurbineParameters>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Helix {
    #[serde(rename = "helix_a")]
    pub a: f64,
    pub helix_power_b: f64,
    pub helix_power_c: f64,
    pub helix_thrust_b: f64,
    pub helix_thrust_c: f64,
}

/// Parameters used by controller-dependent operation models
#[derive(Debug, Clone, Deserialize, Serialize, Default)]
pub struct ControllerDependentTurbineParameters {
    #[serde(default)]
    pub rotor_solidity: Option<f64>,
    #[serde(default)]
    pub rated_rpm: Option<f64>,
    #[serde(default)]
    pub generator_efficiency: Option<f64>,
    #[serde(default)]
    pub rated_power: Option<f64>,
    #[serde(default)]
    pub rotor_diameter: Option<f64>,
    #[serde(default)]
    pub beta: Option<f64>,
    #[serde(default)]
    pub cd: Option<f64>,
    #[serde(default)]
    pub cl_alfa: Option<f64>,
    #[serde(default)]
    pub cp_ct_data_file: Option<String>,
}
pub fn load_turbine_type(path: &str) -> crate::Result<TurbineType> {
    let content = std::fs::read_to_string(path)?;
    load_turbine_type_from_str(&content)
}
fn load_turbine_type_from_str(content: &str) -> crate::Result<TurbineType> {
    let turbine_type: TurbineType = serde_yaml::from_str(content)?;
    Ok(turbine_type)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_turbine_type() {
        let yaml_content = r#"
name: TestTurbine
rotor_diameter: 100.0
hub_height: 100.0
"#;
        let turbine_type = load_turbine_type_from_str(yaml_content).unwrap();
        assert_eq!(turbine_type.name, "TestTurbine");
        assert_eq!(turbine_type.rotor_diameter, 100.0);
        assert_eq!(turbine_type.hub_height, 100.0);
    }
}
