use super::cp_ct_table::{CpCtTable, CsvDataLoader, InterpolationMethod, OneDimTable};
use crate::core::turbines::operation_models::{self, OperationModel};
use crate::core::turbines::TurbineTypeError;
use crate::types::{Array1, Float};
use log::warn;
use serde::de::{Deserializer, MapAccess, Visitor};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::result::Result;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum TurbineError {
    #[error("必需参数缺失: {0}")]
    MissingMandatoryParameter(String),

    #[error("无效参数值: {0} = {1}")]
    InvalidParameterValue(String, f64),

    #[error("Cp/Ct表验证失败: {0}")]
    TableValidationFailed(String),

    #[error("多维度表维度不匹配")]
    DimensionMismatch,

    #[error("插值失败: {0}")]
    InterpolationFailed(String),

    #[error("配置加载失败: {0}")]
    ConfigLoadError(#[from] serde_yaml::Error),

    #[error("文件操作失败: {0}")]
    IoError(#[from] std::io::Error),
}

// 运行条件结构体
#[derive(Debug, Clone, Serialize, Deserialize)]
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
#[derive(Debug, Clone)]
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
    pub operation_model: Box<dyn OperationModel>,
}

impl TurbineType {
    pub fn load_turbine_type(path: &str) -> Result<TurbineType, TurbineTypeError> {
        let content = std::fs::read_to_string(path)?;
        let mut turbine_type: TurbineType = serde_yaml::from_str(&content)?;
        // 检查并加载power_thrust_data_file（如果存在）
        if let Some(ref file_path) = turbine_type.power_thrust_table.power_thrust_data_file {
            // 构造相对于YAML文件的完整路径
            let yaml_dir = std::path::Path::new(path).parent().ok_or_else(|| {
                TurbineTypeError::IoError(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "Could not determine parent directory",
                ))
            })?;
            let full_path = yaml_dir.join(file_path);

            // 加载多维数据
            let multi_dim_table =
                CsvDataLoader::load_multidimensional_data(full_path).map_err(|e| {
                    TurbineTypeError::IoError(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("Failed to load power thrust data file: {}", e),
                    ))
                })?;

            // 更新turbine_type的cp_ct_table
            turbine_type.power_thrust_table.cp_ct_table =
                CpCtTable::MultiDimensional(multi_dim_table);
        }
        Ok(turbine_type)
    }
    fn load_turbine_type_from_str(content: &str) -> Result<TurbineType, TurbineError> {
        let turbine_type: TurbineType = serde_yaml::from_str(content)?;
        Ok(turbine_type)
    }
}

impl<'de> Deserialize<'de> for TurbineType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct Helper {
            rotor_diameter: f64,
            hub_height: f64,
            #[serde(rename = "TSR")]
            tsr: Option<f64>,
            #[serde(rename = "rated_power")]
            rated_power: Option<f64>,
            #[serde(rename = "cut_in_wind_speed")]
            cut_in_wind_speed: Option<f64>,
            #[serde(rename = "cut_out_wind_speed")]
            cut_out_wind_speed: Option<f64>,
            #[serde(rename = "rated_wind_speed")]
            rated_wind_speed: Option<f64>,
            #[serde(rename = "max_yaw_angle")]
            max_yaw_angle: Option<f64>,
            #[serde(rename = "yaw_rate_limit")]
            yaw_rate_limit: Option<f64>,
            #[serde(rename = "generator_efficiency")]
            generator_efficiency: Option<f64>,
            #[serde(rename = "drive_train_efficiency")]
            drive_train_efficiency: Option<f64>,
            #[serde(rename = "additional_losses")]
            additional_losses: Option<f64>,
            #[serde(rename = "ti_ref")]
            ti_ref: Option<f64>,
            #[serde(rename = "power_thrust_table")]
            power_thrust_table: PowerTrustParams,
            #[serde(rename = "floating_tilt_table")]
            floating_tilt_table: Option<FloatingTiltTable>,
            #[serde(rename = "correct_cp_ct_for_tilt")]
            correct_cp_ct_for_tilt: Option<bool>,
            #[serde(rename = "multi_dimensional_cp_ct")]
            multi_dimensional_cp_ct: Option<bool>,
            #[serde(rename = "turbine_type")]
            turbine_type: String,
            #[serde(rename = "operation_model")]
            operation_model: Option<String>,
            #[serde(rename = "power_thrust_data_file")]
            power_thrust_data_file: Option<String>,
        }
        let helper = Helper::deserialize(deserializer)?;
        let operation_model = helper
            .operation_model
            .unwrap_or_else(|| "cosine-loss".to_string());
        let om = match operation_model.as_str() {
            "simple" => {
                let simple_turbine: Box<dyn OperationModel + 'static> =
                    Box::new(operation_models::simple::SimpleTurbine);
                simple_turbine
            }
            "cosine-loss" | _ => Box::new(operation_models::cosine_loss::CosineLossTurbine)
                as Box<dyn OperationModel>,
            "simple-derating" => Box::new(operation_models::simple_derating::SimpleDeratingTurbine)
                as Box<dyn OperationModel>,
            "peak-shaving" => Box::new(operation_models::peak_shaving::PeakShavingTurbine)
                as Box<dyn OperationModel>,

            "mixed" => {
                Box::new(operation_models::mixed::MixedOperationTurbine) as Box<dyn OperationModel>
            }
            "awc" => Box::new(operation_models::awc::AWCTurbine) as Box<dyn OperationModel>,

            "unified-momentum" => {
                Box::new(operation_models::unified_momentum::UnifiedMomentumTurbine)
                    as Box<dyn OperationModel>
            }
            "controller-dependent" => {
                Box::new(operation_models::controller_dependent::ControllerDependentTurbine)
                    as Box<dyn OperationModel>
            }

            _ => {
                warn!(
                    "Unknown operation model specified: {}. Defaulting to CosineLossTurbine.",
                    operation_model
                );
                Box::new(operation_models::cosine_loss::CosineLossTurbine)
                    as Box<dyn OperationModel>
            }
        };

        Ok(TurbineType {
            name: helper.turbine_type,
            rotor_diameter: helper.rotor_diameter,
            hub_height: helper.hub_height,
            tsr: helper.tsr,
            rated_power: helper.rated_power,
            cut_in_wind_speed: helper.cut_in_wind_speed,
            cut_out_wind_speed: helper.cut_out_wind_speed,
            rated_wind_speed: helper.rated_wind_speed,
            max_yaw_angle: helper.max_yaw_angle,
            yaw_rate_limit: helper.yaw_rate_limit,
            generator_efficiency: helper.generator_efficiency,
            drive_train_efficiency: helper.drive_train_efficiency,
            additional_losses: helper.additional_losses,
            ti_ref: helper.ti_ref,
            power_thrust_table: helper.power_thrust_table,
            floating_tilt_table: helper.floating_tilt_table,
            correct_cp_ct_for_tilt: helper.correct_cp_ct_for_tilt.unwrap_or(false),
            operation_model: om,
        })
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct FloatingTiltTable {
    pub wind_speed: Array1,
    pub tilt_angle: Array1,
}

impl<'de> Deserialize<'de> for FloatingTiltTable {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct Helper {
            wind_speed: Vec<Float>,
            tilt: Vec<Float>, // 注意：YAML中是tilt，而不是tilt_angle
        }

        let helper = Helper::deserialize(deserializer)?;
        Ok(FloatingTiltTable {
            wind_speed: ndarray::Array1::from_vec(helper.wind_speed),
            tilt_angle: ndarray::Array1::from_vec(helper.tilt),
        })
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PowerTrustParams {
    #[serde(skip_deserializing)]
    pub cp_ct_table: CpCtTable,
    pub ref_air_density: Option<f64>,
    pub ref_tilt: Option<f64>,
    pub cosine_loss_exponent_yaw: Option<f64>,
    pub cosine_loss_exponent_tilt: Option<f64>,
    pub generator_efficiency: Float, // 发电机效率
    pub ti_ref: Float,               // 参考湍流强度
    // Peak shaving parameters
    pub peak_shaving_fraction: Option<f64>,
    pub peak_shaving_ti_threshold: Option<f64>,
    // Power thrust data file for multi-dimensional CP/CT
    pub power_thrust_data_file: Option<String>,
    pub helix: Option<Helix>,
    pub controller_dependent_turbine_parameters: Option<ControllerDependentTurbineParameters>,
}

impl<'de> Deserialize<'de> for PowerTrustParams {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct PowerTrustParamsVisitor;

        impl<'de> Visitor<'de> for PowerTrustParamsVisitor {
            type Value = PowerTrustParams;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a PowerTrustParams map")
            }

            fn visit_map<V>(self, mut map: V) -> Result<PowerTrustParams, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut ref_air_density = None;
                let mut ref_tilt = None;
                let mut cosine_loss_exponent_yaw = None;
                let mut cosine_loss_exponent_tilt = None;
                let mut generator_efficiency = None;
                let mut ti_ref = None;
                let mut peak_shaving_fraction = None;
                let mut peak_shaving_ti_threshold = None;
                let mut power_thrust_data_file = None;
                let mut helix = None;
                let mut controller_dependent_turbine_parameters = None;

                // 存储wind_speed, power, thrust_coefficient等数据以便后续处理
                let mut wind_speed: Option<Vec<Float>> = None;
                let mut power: Option<Vec<Float>> = None;
                let mut thrust_coefficient: Option<Vec<Float>> = None;
                let mut multi_dimensional_cp_ct: Option<bool> = None;
                let mut external_conditions: Option<Vec<Vec<Float>>> = None;

                // Helix 相关字段
                let mut helix_a: Option<f64> = None;
                let mut helix_power_b: Option<f64> = None;
                let mut helix_power_c: Option<f64> = None;
                let mut helix_thrust_b: Option<f64> = None;
                let mut helix_thrust_c: Option<f64> = None;

                while let Some(key) = map.next_key::<String>()? {
                    match key.as_str() {
                        "ref_air_density" => {
                            ref_air_density = Some(map.next_value()?);
                        }
                        "ref_tilt" => {
                            ref_tilt = Some(map.next_value()?);
                        }
                        "cosine_loss_exponent_yaw" => {
                            cosine_loss_exponent_yaw = Some(map.next_value()?);
                        }
                        "cosine_loss_exponent_tilt" => {
                            cosine_loss_exponent_tilt = Some(map.next_value()?);
                        }
                        "generator_efficiency" => {
                            generator_efficiency = map.next_value()?;
                        }
                        "ti_ref" => {
                            ti_ref = map.next_value()?;
                        }
                        "peak_shaving_fraction" => {
                            peak_shaving_fraction = Some(map.next_value()?);
                        }
                        "peak_shaving_ti_threshold" | "peak_shaving_TI_threshold" => {
                            peak_shaving_ti_threshold = Some(map.next_value()?);
                        }
                        "power_thrust_data_file" => {
                            power_thrust_data_file = Some(map.next_value()?);
                        }

                        "controller_dependent_turbine_parameters" => {
                            controller_dependent_turbine_parameters = Some(map.next_value()?);
                        }
                        "wind_speed" => {
                            wind_speed = Some(map.next_value()?);
                        }
                        "power" => {
                            power = Some(map.next_value()?);
                        }
                        "thrust_coefficient" => {
                            thrust_coefficient = Some(map.next_value()?);
                        }
                        // 处理 helix_ 开头的字段
                        "helix_a" => {
                            helix_a = Some(map.next_value()?);
                        }
                        "helix_power_b" => {
                            helix_power_b = Some(map.next_value()?);
                        }
                        "helix_power_c" => {
                            helix_power_c = Some(map.next_value()?);
                        }
                        "helix_thrust_b" => {
                            helix_thrust_b = Some(map.next_value()?);
                        }
                        "helix_thrust_c" => {
                            helix_thrust_c = Some(map.next_value()?);
                        }
                        // 如果multi_dimensional_cp_ct在PowerTrustParams内部，也处理它
                        _ => {
                            // 忽略未知字段或稍后处理
                            let _: serde::de::IgnoredAny = map.next_value()?;
                        }
                    }
                }
                // 构建 Helix 结构体，如果所有字段都存在的话
                if helix_a.is_some()
                    || helix_power_b.is_some()
                    || helix_power_c.is_some()
                    || helix_thrust_b.is_some()
                    || helix_thrust_c.is_some()
                {
                    helix = Some(Helix {
                        a: helix_a.unwrap_or(0.0),
                        helix_power_b: helix_power_b.unwrap_or(0.0),
                        helix_power_c: helix_power_c.unwrap_or(0.0),
                        helix_thrust_b: helix_thrust_b.unwrap_or(0.0),
                        helix_thrust_c: helix_thrust_c.unwrap_or(0.0),
                    });
                }
                // 根据wind_speed, power, thrust_coefficient创建CpCtTable
                let cp_ct_table = if let (Some(ws), Some(pow), Some(thrust)) =
                    (&wind_speed, &power, &thrust_coefficient)
                {
                    // 将功率转换为Cp值，推力系数就是Ct值
                    // 这里假设power是功率系数，如果实际是功率需要转换
                    CpCtTable::OneDimensional(OneDimTable {
                        wind_speeds: ndarray::Array1::from_vec(ws.clone()),
                        cp_values: ndarray::Array1::from_vec(pow.clone()),
                        ct_values: ndarray::Array1::from_vec(thrust.clone()),
                        interpolation: InterpolationMethod::Linear, // 默认插值方法
                    })
                } else if power_thrust_data_file.is_some() {
                    warn!(
                        "Power thrust data file specified but not loaded yet: {}",
                        power_thrust_data_file.as_ref().unwrap()
                    );
                    // 创建一个空的表格作为占位符，将在load_turbine_type函数中更新
                    CpCtTable::OneDimensional(OneDimTable {
                        wind_speeds: ndarray::Array1::zeros(0),
                        cp_values: ndarray::Array1::zeros(0),
                        ct_values: ndarray::Array1::zeros(0),
                        interpolation: InterpolationMethod::Linear,
                    })
                } else {
                    warn!(
                        "Power/Thrust data not provided in YAML, and no data file specified. Using empty Cp/Ct table as default."
                    );
                    // 创建一个空的表格作为默认值
                    CpCtTable::OneDimensional(OneDimTable {
                        wind_speeds: ndarray::Array1::zeros(0),
                        cp_values: ndarray::Array1::zeros(0),
                        ct_values: ndarray::Array1::zeros(0),
                        interpolation: InterpolationMethod::Linear,
                    })
                };

                // 返回一个临时的PowerTrustParams，其中cp_ct_table需要稍后设置
                Ok(PowerTrustParams {
                    cp_ct_table,
                    ref_air_density,
                    ref_tilt,
                    cosine_loss_exponent_yaw,
                    cosine_loss_exponent_tilt,
                    generator_efficiency: generator_efficiency.unwrap_or(0.96),
                    ti_ref: ti_ref.unwrap_or(0.06),
                    peak_shaving_fraction,
                    peak_shaving_ti_threshold,
                    power_thrust_data_file,
                    helix,
                    controller_dependent_turbine_parameters,
                })
            }
        }

        deserializer.deserialize_map(PowerTrustParamsVisitor)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Helix {
    #[serde(rename = "helix_a")]
    pub a: f64,
    pub helix_power_b: f64, // 字段名与YAML中的键名一致
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

#[cfg(test)]
mod tests {
    use crate::core::turbines::cp_ct_table::Dimension;
    use crate::core::turbines::cp_ct_table::TableType;

    use super::*;

    #[test]
    fn test_load_turbine_type_nrel_5mw() {
        let yaml_path = "turbine_library/nrel_5MW.yaml";
        let turbine_type = TurbineType::load_turbine_type(yaml_path);
        assert!(turbine_type.is_ok());
        let turbine_type = turbine_type.unwrap();

        // 验证基本参数
        assert_eq!(turbine_type.name, "nrel_5MW");
        assert_eq!(turbine_type.hub_height, 90.0);
        assert_eq!(turbine_type.rotor_diameter, 125.88);
        assert_eq!(turbine_type.tsr, Some(8.0));

        // 验证 power_thrust_table 参数
        let pt_table = &turbine_type.power_thrust_table;
        assert_eq!(pt_table.ref_air_density, Some(1.225));
        assert_eq!(pt_table.ref_tilt, Some(5.0));
        assert_eq!(pt_table.cosine_loss_exponent_tilt, Some(1.88));
        assert_eq!(pt_table.cosine_loss_exponent_yaw, Some(1.88));

        // 验证 helix 参数
        assert!(pt_table.helix.is_some());
        if let Some(helix) = &pt_table.helix {
            assert_eq!(helix.a, 1.802);
            assert_eq!(helix.helix_power_b, 4.568e-03);
            assert_eq!(helix.helix_power_c, 1.629e-10);
            assert_eq!(helix.helix_thrust_b, 1.027e-03);
            assert_eq!(helix.helix_thrust_c, 1.378e-06);
        }

        // 验证 peak shaving 参数
        assert_eq!(pt_table.peak_shaving_fraction, Some(0.2));
        assert_eq!(pt_table.peak_shaving_ti_threshold, Some(0.1));

        // 验证 controller_dependent_turbine_parameters
        if let Some(ctrl_params) = &pt_table.controller_dependent_turbine_parameters {
            assert_eq!(ctrl_params.rated_rpm, Some(12.1));
            assert_eq!(ctrl_params.rotor_solidity, Some(0.05132));
            assert_eq!(ctrl_params.generator_efficiency, Some(0.944));
            assert_eq!(ctrl_params.rated_power, Some(5000.0));
            assert_eq!(ctrl_params.rotor_diameter, Some(126.0));
            assert_eq!(ctrl_params.beta, Some(-0.45891));
            assert_eq!(ctrl_params.cd, Some(0.0040638));
            assert_eq!(ctrl_params.cl_alfa, Some(4.275049));
            assert_eq!(
                ctrl_params.cp_ct_data_file,
                Some("demo_cp_ct_surfaces/nrel_5MW_demo_cp_ct_surface.npz".to_string())
            );
        }

        // 验证其他参数
        assert_eq!(turbine_type.correct_cp_ct_for_tilt, false);

        // 验证 floating_tilt_table
        if let Some(floating_tilt) = &turbine_type.floating_tilt_table {
            assert_eq!(floating_tilt.wind_speed.len(), 7); // 7个风速值
            assert_eq!(floating_tilt.tilt_angle.len(), 7); // 7个倾斜角度值
            assert_eq!(floating_tilt.wind_speed[0], 4.0);
            assert_eq!(floating_tilt.tilt_angle[0], 5.0);
        }
    }
    #[test]
    fn test_load_turbine_type_iea_10mw() {
        let yaml_path = "turbine_library/iea_10MW.yaml";
        let turbine_type_result = TurbineType::load_turbine_type(yaml_path);
        if turbine_type_result.as_ref().is_err() {
            println!(
                "Error loading turbine type: {}",
                turbine_type_result.as_ref().err().unwrap()
            );
        }
        assert!(turbine_type_result.is_ok());
        let turbine_type = turbine_type_result.unwrap();

        // 验证基本参数
        assert_eq!(turbine_type.name, "iea_10MW");
        assert_eq!(turbine_type.hub_height, 119.0);
        assert_eq!(turbine_type.rotor_diameter, 198.0);
        assert_eq!(turbine_type.tsr, Some(8.0));

        // 验证 power_thrust_table 参数
        let pt_table = &turbine_type.power_thrust_table;
        assert_eq!(
            pt_table.cp_ct_table,
            CpCtTable::OneDimensional(OneDimTable {
                wind_speeds: ndarray::Array1::from_vec(vec![
                    0.0000, 2.9, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5,
                    12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0, 25.01, 50.0
                ]),
                cp_values: ndarray::Array1::from_vec(vec![
                    0.0, 0.0, 35.60156, 414.0606, 1009.90686, 1855.02326, 2963.01442, 4440.26484,
                    6330.82856, 7392.13274, 8514.32824, 9691.10578, 10000.00, 10000.00, 10000.00,
                    10000.00, 10000.00, 10000.00, 10000.00, 10000.00, 10000.00, 10000.00, 0.0, 0.0
                ]),
                ct_values: ndarray::Array1::from_vec(vec![
                    0.0, 0.0, 0.915, 0.926, 0.921, 0.895, 0.885, 0.873, 0.827, 0.789, 0.754, 0.721,
                    0.591, 0.49, 0.418, 0.318, 0.251, 0.203, 0.167, 0.119, 0.088, 0.049, 0.0, 0.0
                ]),
                interpolation: InterpolationMethod::Linear,
            })
        );
        assert_eq!(pt_table.ref_air_density, Some(1.225));
        assert_eq!(pt_table.ref_tilt, Some(6.0));
        assert_eq!(pt_table.cosine_loss_exponent_yaw, Some(1.88));
        assert_eq!(pt_table.cosine_loss_exponent_tilt, Some(1.88));

        // 验证 helix 参数
        assert!(pt_table.helix.is_some());
        if let Some(helix) = &pt_table.helix {
            assert_eq!(helix.a, 1.719);
            assert_eq!(helix.helix_power_b, 4.823e-03);
            assert_eq!(helix.helix_power_c, 2.314e-10);
            assert_eq!(helix.helix_thrust_b, 1.157e-03);
            assert_eq!(helix.helix_thrust_c, 1.167e-04);
        }

        // 验证 controller_dependent_turbine_parameters
        if let Some(ctrl_params) = &pt_table.controller_dependent_turbine_parameters {
            assert_eq!(ctrl_params.rotor_solidity, Some(0.03500415472147307));
            assert_eq!(ctrl_params.rated_rpm, Some(8.6));
            assert_eq!(ctrl_params.generator_efficiency, Some(0.944));
            assert_eq!(ctrl_params.rated_power, Some(10000.0));
            assert_eq!(ctrl_params.rotor_diameter, Some(198.0));
            assert_eq!(ctrl_params.beta, Some(-3.8233819218614817));
            assert_eq!(ctrl_params.cd, Some(0.004612981322772105));
            assert_eq!(ctrl_params.cl_alfa, Some(4.602140680380394));
            assert_eq!(
                ctrl_params.cp_ct_data_file,
                Some("demo_cp_ct_surfaces/iea_10MW_demo_cp_ct_surface.npz".to_string())
            );
        }

        // 验证其他参数
        assert_eq!(turbine_type.correct_cp_ct_for_tilt, false);
    }
    #[test]
    fn test_load_turbine_type_iea_15mw() {
        let yaml_path = "turbine_library/iea_15MW.yaml";
        let turbine_type_result = TurbineType::load_turbine_type(yaml_path);
        if turbine_type_result.as_ref().is_err() {
            println!(
                "Error loading turbine type: {}",
                turbine_type_result.as_ref().err().unwrap()
            );
        }
        assert!(turbine_type_result.is_ok());
        let turbine_type = turbine_type_result.unwrap();

        // 验证基本参数
        assert_eq!(turbine_type.name, "iea_15MW");
        assert_eq!(turbine_type.hub_height, 150.0);
        assert_eq!(turbine_type.rotor_diameter, 242.24);
        assert_eq!(turbine_type.tsr, Some(8.0));

        // 验证 power_thrust_table 参数
        let pt_table = &turbine_type.power_thrust_table;
        assert_eq!(pt_table.ref_air_density, Some(1.225));
        assert_eq!(pt_table.ref_tilt, Some(6.0));
        assert_eq!(pt_table.cosine_loss_exponent_yaw, Some(1.88));
        assert_eq!(pt_table.cosine_loss_exponent_tilt, Some(1.88));

        // 验证 helix 参数
        assert!(pt_table.helix.is_some());
        if let Some(helix) = &pt_table.helix {
            assert_eq!(helix.a, 1.809);
            assert_eq!(helix.helix_power_b, 4.828e-03);
            assert_eq!(helix.helix_power_c, 4.017e-11);
            assert_eq!(helix.helix_thrust_b, 1.390e-03);
            assert_eq!(helix.helix_thrust_c, 5.084e-04);
        }

        // 验证 controller_dependent_turbine_parameters
        if let Some(ctrl_params) = &pt_table.controller_dependent_turbine_parameters {
            assert_eq!(ctrl_params.rotor_solidity, Some(0.031018237027995298));
            assert_eq!(ctrl_params.rated_rpm, Some(7.55));
            assert_eq!(ctrl_params.generator_efficiency, Some(0.95756));
            assert_eq!(ctrl_params.rated_power, Some(15000.00));
            assert_eq!(ctrl_params.rotor_diameter, Some(242.24));
            assert_eq!(ctrl_params.beta, Some(-3.098605491003358));
            assert_eq!(ctrl_params.cd, Some(0.004426686198054057));
            assert_eq!(ctrl_params.cl_alfa, Some(4.546410770937916));
            assert_eq!(
                ctrl_params.cp_ct_data_file,
                Some("demo_cp_ct_surfaces/iea_15MW_demo_cp_ct_surface.npz".to_string())
            );
        }

        // 验证其他参数
        assert_eq!(turbine_type.correct_cp_ct_for_tilt, false);
    }

    #[test]
    fn test_load_turbine_type_iea_22mw() {
        let yaml_path = "turbine_library/iea_22MW.yaml";
        let turbine_type_result = TurbineType::load_turbine_type(yaml_path);
        if turbine_type_result.as_ref().is_err() {
            println!(
                "Error loading turbine type: {}",
                turbine_type_result.as_ref().err().unwrap()
            );
        }
        assert!(turbine_type_result.is_ok());
        let turbine_type = turbine_type_result.unwrap();

        // 验证基本参数
        assert_eq!(turbine_type.name, "iea_22MW");
        assert_eq!(turbine_type.hub_height, 170.0);
        assert_eq!(turbine_type.rotor_diameter, 284.0);
        assert_eq!(turbine_type.tsr, Some(9.15));

        // 验证 power_thrust_table 参数
        let pt_table = &turbine_type.power_thrust_table;
        assert_eq!(pt_table.ref_air_density, Some(1.225));
        assert_eq!(pt_table.ref_tilt, Some(6.0));
        assert_eq!(pt_table.cosine_loss_exponent_yaw, Some(1.88));
        assert_eq!(pt_table.cosine_loss_exponent_tilt, Some(1.88));

        // IEA 22MW doesn't seem to have helix parameters in the YAML file
        assert!(pt_table.helix.is_none());

        // 验证 controller_dependent_turbine_parameters (如果存在)
        if let Some(ctrl_params) = &pt_table.controller_dependent_turbine_parameters {
            // 检查一些预期的值
            assert_eq!(ctrl_params.rated_power, Some(22000.0));
        }

        // 验证其他参数
        assert_eq!(turbine_type.correct_cp_ct_for_tilt, false);
    }
    #[test]
    fn test_load_turbine_type_iea_15mw_multi_dim_cp_ct() {
        let yaml_path = "turbine_library/iea_15MW_multi_dim_cp_ct.yaml";
        let turbine_type_result = TurbineType::load_turbine_type(yaml_path);
        if turbine_type_result.as_ref().is_err() {
            println!(
                "Error loading turbine type: {}",
                turbine_type_result.as_ref().err().unwrap()
            );
        }
        assert!(turbine_type_result.is_ok());
        let turbine_type = turbine_type_result.unwrap();

        // 验证基本参数
        assert_eq!(turbine_type.name, "iea_15MW_multi_dim_cp_ct");
        assert_eq!(turbine_type.hub_height, 150.0);
        assert_eq!(turbine_type.rotor_diameter, 242.24);
        assert_eq!(turbine_type.tsr, Some(8.0));

        // 验证多维CP/CT标志
        // 注意：由于multi_dimensional_cp_ct在Helper结构中定义，但在最终的TurbineType中没有保存，
        // 所以我们不能直接验证这个值

        // 验证 power_thrust_table 参数
        let pt_table = &turbine_type.power_thrust_table;
        assert_eq!(pt_table.ref_air_density, Some(1.225));
        assert_eq!(pt_table.ref_tilt, Some(6.0));
        assert_eq!(pt_table.cosine_loss_exponent_yaw, Some(1.88));
        assert_eq!(pt_table.cosine_loss_exponent_tilt, Some(1.88));
        assert_eq!(
            pt_table.power_thrust_data_file,
            Some("iea_15MW_multi_dim_Tp_Hs.csv".to_string())
        );

        assert_eq!(
            pt_table.cp_ct_table.table_type(),
            TableType::MultiDimensional
        );
        let supported_dims = &pt_table.cp_ct_table.supported_dimensions();
        assert!(supported_dims.contains(&Dimension::WindSpeed));
        assert!(supported_dims.contains(&Dimension::WaveHeight));
        assert!(supported_dims.contains(&Dimension::WavePeriod));
        assert_eq!(supported_dims.len(), 3);

        // 验证其他参数
        assert_eq!(turbine_type.correct_cp_ct_for_tilt, false);
    }

    #[test]
    fn test_load_turbine_type_iea_15mw_floating_multi_dim_cp_ct() {
        let yaml_path = "turbine_library/iea_15MW_floating_multi_dim_cp_ct.yaml";
        let turbine_type_result = TurbineType::load_turbine_type(yaml_path);
        if turbine_type_result.as_ref().is_err() {
            println!(
                "Error loading turbine type: {}",
                turbine_type_result.as_ref().err().unwrap()
            );
        }
        assert!(turbine_type_result.is_ok());
        let turbine_type = turbine_type_result.unwrap();

        // 验证基本参数
        assert_eq!(turbine_type.name, "iea_15MW_floating");
        assert_eq!(turbine_type.hub_height, 150.0);
        assert_eq!(turbine_type.rotor_diameter, 242.24);
        assert_eq!(turbine_type.tsr, Some(8.0));

        // 验证 power_thrust_table 参数
        let pt_table = &turbine_type.power_thrust_table;
        assert_eq!(pt_table.ref_air_density, Some(1.225));
        assert_eq!(pt_table.ref_tilt, Some(6.0));
        assert_eq!(pt_table.cosine_loss_exponent_yaw, Some(1.88));
        assert_eq!(pt_table.cosine_loss_exponent_tilt, Some(1.88));
        assert_eq!(
            pt_table.power_thrust_data_file,
            Some("iea_15MW_multi_dim_Tp_Hs.csv".to_string())
        );

        // 验证 floating_tilt_table
        assert!(turbine_type.floating_tilt_table.is_some());
        if let Some(floating_tilt) = &turbine_type.floating_tilt_table {
            assert_eq!(floating_tilt.wind_speed.len(), 22); // 22个风速值
            assert_eq!(floating_tilt.tilt_angle.len(), 22); // 22个倾斜角度值
            assert_eq!(floating_tilt.wind_speed[0], 3.5);
            assert_eq!(floating_tilt.wind_speed[21], 24.5);
            assert_eq!(floating_tilt.tilt_angle[0], 5.406938261);
            assert_eq!(floating_tilt.tilt_angle[21], 7.45510389);
        }

        // 验证 correct_cp_ct_for_tilt
        assert_eq!(turbine_type.correct_cp_ct_for_tilt, true);
    }
}
