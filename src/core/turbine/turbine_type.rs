use serde::{Deserialize, Serialize};
use serde_yaml::Value;
use std::fs::File;
use std::fmt;

use crate::Array1;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct LookupTable {
    pub wind_speeds: Array1,
    pub values: Array1,
}

pub type PowerTable = LookupTable;
pub type ThrustTable = LookupTable;

impl LookupTable {
    pub fn interpolate(&self, wind_speed: f64) -> f64 {
        let ws = &self.wind_speeds;
        let n = ws.len();

        if n == 0 {
            return 0.0;
        }

        if wind_speed <= ws[0] {
            return self.values[0];
        }

        if wind_speed >= ws[n - 1] {
            return self.values[n - 1];
        }

        let mut lo = 0;
        let mut hi = n - 1;

        while lo < hi {
            let mid = (lo + hi) / 2;
            if wind_speed < ws[mid] {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }

        if lo >= n {
            return self.values[n - 1];
        }

        let lo_val = if lo > 0 {
            self.values[lo - 1]
        } else {
            self.values[0]
        };
        let _hi_val = self.values[lo];

        if hi == lo {
            return lo_val;
        }

        let x0 = ws[lo];
        let x1 = ws[hi];
        let y0 = self.values[lo];
        let y1 = self.values[hi];

        y0 + (y1 - y0) * (wind_speed - x0) / (x1 - x0)
    }
}

// 直接反序列化的TurbineType - 支持flat格式
#[derive(Debug, Clone, Deserialize)]
pub struct TurbineType {
    #[serde(rename = "turbine_type")]
    pub name: String,
    pub rotor_diameter: f64,
    pub hub_height: f64,
    #[serde(rename = "TSR")]
    pub tsr: f64,
    pub operation_model: String, // 使用String，通过方法转换为enum
    #[serde(default)]
    pub ref_tilt: Option<f64>,
    #[serde(default)]
    pub correct_cp_ct_for_tilt: Option<bool>,
    pub power_curve_wind_speeds: Vec<f64>,
    pub power_curve_powers: Vec<f64>,
    pub thrust_coefficient_wind_speeds: Vec<f64>,
    pub thrust_coefficient_values: Vec<f64>,
    #[serde(default)]
    pub controller_dependent_turbine_parameters: Option<Value>,
}

impl fmt::Display for TurbineType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}

impl TurbineType {
    /// 解析operation_model字符串为OperationModel枚举
    pub fn get_operation_model_enum(&self) -> OperationModel {
        match self.operation_model.as_str() {
            "cosine-loss" => OperationModel::CosineLoss,
            "null" => OperationModel::Null,
            "p8" => OperationModel::P8,
            "simple-yaw" => OperationModel::SimpleYaw,
            _ => OperationModel::Unknown,
        }
    }

    pub fn power_curve(&self) -> LookupTable {
        LookupTable {
            wind_speeds: ndarray::Array1::from(self.power_curve_wind_speeds.clone()),
            values: ndarray::Array1::from(self.power_curve_powers.clone()),
        }
    }

    pub fn thrust_curve(&self) -> LookupTable {
        LookupTable {
            wind_speeds: ndarray::Array1::from(self.thrust_coefficient_wind_speeds.clone()),
            values: ndarray::Array1::from(self.thrust_coefficient_values.clone()),
        }
    }

    pub fn get_ct(&self, wind_speed: f64) -> f64 {
        self.thrust_curve().interpolate(wind_speed)
    }

    pub fn get_power(&self, wind_speed: f64, yaw_angle_rad: f64) -> f64 {
        let base_power = self.power_curve().interpolate(wind_speed);
        let loss_factor = self
            .get_operation_model_enum()
            .power_loss_factor(yaw_angle_rad, 2.0);
        base_power * loss_factor
    }

    pub fn get_controller_param(&self, key: &str) -> Option<&Value> {
        self.controller_dependent_turbine_parameters
            .as_ref()?
            .get(key)
    }
}

pub fn load_turbine_type(path: &str) -> Result<TurbineType, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let turbine: TurbineType = serde_yaml::from_reader(file)?;
    Ok(turbine)
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "kebab-case")]
pub enum OperationModel {
    Null,
    CosineLoss,
    P8,
    SimpleYaw,
    #[serde(other)]
    Unknown,
}

impl OperationModel {
    pub fn power_loss_factor(&self, yaw_rad: f64, exponent: f64) -> f64 {
        match self {
            Self::CosineLoss => yaw_rad.cos().powf(exponent),
            Self::Null => 1.0,
            _ => unimplemented!("Model not implemented in Rust"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_nrel_5mw_turbine() {
        let result = load_turbine_type("turbine_library/nrel_5MW.yaml");
        assert!(result.is_ok(), "Failed to load: {:?}", result.err());

        let turbine = result.unwrap();
        assert_eq!(turbine.name, "nrel_5MW");
        assert_eq!(turbine.rotor_diameter, 126.0);
        assert_eq!(turbine.hub_height, 90.0);
        assert_eq!(turbine.tsr, 8.0);
        assert_eq!(turbine.operation_model, "cosine-loss");

        assert_eq!(turbine.power_curve_wind_speeds.len(), 54);
        assert_eq!(turbine.thrust_coefficient_wind_speeds.len(), 54);
    }

    #[test]
    fn test_load_iea_15mw_turbine() {
        let result = load_turbine_type("turbine_library/iea_15MW.yaml");
        assert!(result.is_ok(), "Failed to load: {:?}", result.err());

        let turbine = result.unwrap();
        assert_eq!(turbine.name, "iea_15MW");
        assert_eq!(turbine.rotor_diameter, 240.0);
        assert_eq!(turbine.hub_height, 150.0);
        assert_eq!(turbine.tsr, 8.0);
        assert_eq!(turbine.operation_model, "cosine-loss");
    }

    #[test]
    fn test_load_iea_10mw_turbine() {
        let result = load_turbine_type("turbine_library/iea_10MW.yaml");
        assert!(result.is_ok(), "Failed to load: {:?}", result.err());

        let turbine = result.unwrap();
        assert_eq!(turbine.name, "iea_10MW");
    }

    #[test]
    fn test_turbine_power_calculation() {
        let turbine = TurbineType {
            name: "test".to_string(),
            rotor_diameter: 100.0,
            hub_height: 80.0,
            tsr: 7.0,
            operation_model: "null".to_string(),
            power_curve_wind_speeds: vec![0.0, 5.0, 10.0, 15.0],
            power_curve_powers: vec![0.0, 500.0, 2000.0, 4500.0],
            thrust_coefficient_wind_speeds: vec![0.0, 5.0, 10.0, 15.0],
            thrust_coefficient_values: vec![0.0, 0.8, 0.6, 0.3],
            ref_tilt: None,
            correct_cp_ct_for_tilt: None,
            controller_dependent_turbine_parameters: None,
        };

        let power = turbine.get_power(10.0, 0.0);
        assert_eq!(power, 2000.0);
    }

    #[test]
    fn test_turbine_thrust_coefficient_calculation() {
        let turbine = TurbineType {
            name: "test".to_string(),
            rotor_diameter: 100.0,
            hub_height: 80.0,
            tsr: 7.0,
            operation_model: "null".to_string(),
            power_curve_wind_speeds: vec![0.0, 5.0, 10.0, 15.0],
            power_curve_powers: vec![0.0, 500.0, 2000.0, 4500.0],
            thrust_coefficient_wind_speeds: vec![0.0, 5.0, 10.0, 15.0],
            thrust_coefficient_values: vec![0.0, 0.8, 0.6, 0.3],
            ref_tilt: None,
            correct_cp_ct_for_tilt: None,
            controller_dependent_turbine_parameters: None,
        };

        let ct = turbine.get_ct(10.0);
        assert_eq!(ct, 0.6);
    }

    #[test]
    fn test_cosine_loss_power() {
        let turbine = TurbineType {
            name: "test".to_string(),
            rotor_diameter: 100.0,
            hub_height: 80.0,
            tsr: 7.0,
            operation_model: "cosine-loss".to_string(),
            power_curve_wind_speeds: vec![0.0, 5.0, 10.0, 15.0],
            power_curve_powers: vec![0.0, 500.0, 2000.0, 4500.0],
            thrust_coefficient_wind_speeds: vec![0.0, 5.0, 10.0, 15.0],
            thrust_coefficient_values: vec![0.0, 0.8, 0.6, 0.3],
            ref_tilt: None,
            correct_cp_ct_for_tilt: None,
            controller_dependent_turbine_parameters: None,
        };

        let power = turbine.get_power(10.0, 20.0f64.to_radians());
        assert!((power - 2000.0 * 0.883).abs() < 0.1);
    }
}
