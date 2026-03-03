use serde::{Deserialize, Serialize};
use serde_yaml::Value;
use std::fmt;
use std::fs::File;

use crate::Array1;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct LookupTable {
    pub keys: Array1,
    pub values: Array1,
}

pub type PowerTable = LookupTable;
pub type ThrustTable = LookupTable;
pub type TiltTable = LookupTable;

impl LookupTable {
    pub fn interpolate(&self, wind_speed: f64) -> f64 {
        let ws = &self.keys;
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

// 直接反序列化的TurbineType - 支持Python FLORIS格式
#[derive(Debug, Clone, Deserialize)]
pub struct TurbineType {
    #[serde(rename = "turbine_type")]
    pub name: String,
    pub rotor_diameter: f64,
    pub hub_height: f64,
    #[serde(rename = "TSR")]
    pub tsr: f64,
    pub operation_model: String,
    #[serde(default)]
    pub ref_tilt: Option<f64>,
    #[serde(default)]
    pub correct_cp_ct_for_tilt: Option<bool>,
    #[serde(default)]
    pub power_thrust_table: Option<PowerThrustTable>,
    // Legacy flat format support
    #[serde(default)]
    pub power_curve_wind_speeds: Vec<f64>,
    #[serde(default)]
    pub power_curve_powers: Vec<f64>,
    #[serde(default)]
    pub thrust_coefficient_wind_speeds: Vec<f64>,
    #[serde(default)]
    pub thrust_coefficient_values: Vec<f64>,
    #[serde(default)]
    pub controller_dependent_turbine_parameters: Option<Value>,

    #[serde(default)]
    pub floating_tilt_table: Option<TiltTable>,

    pub multi_dimensional_cp_ct: Option<bool>,
}

/// Power and thrust coefficient table matching Python FLORIS format
#[derive(Debug, Clone, Deserialize)]
pub struct PowerThrustTable {
    #[serde(default)]
    pub wind_speed: Vec<f64>,
    #[serde(default)]
    pub power: Vec<f64>,
    #[serde(default)]
    pub thrust_coefficient: Vec<f64>,
    #[serde(default)]
    pub ref_air_density: Option<f64>,
    #[serde(default)]
    pub ref_tilt: Option<f64>,
    #[serde(default)]
    pub cosine_loss_exponent_yaw: Option<f64>,
    #[serde(default)]
    pub cosine_loss_exponent_tilt: Option<f64>,
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

    /// Get wind speeds from either nested power_thrust_table or legacy flat format
    fn get_wind_speeds(&self) -> Vec<f64> {
        if let Some(ref table) = self.power_thrust_table {
            if !table.wind_speed.is_empty() {
                return table.wind_speed.clone();
            }
        }
        if !self.power_curve_wind_speeds.is_empty() {
            return self.power_curve_wind_speeds.clone();
        }
        if !self.thrust_coefficient_wind_speeds.is_empty() {
            return self.thrust_coefficient_wind_speeds.clone();
        }
        vec![]
    }

    /// Get power values from either nested power_thrust_table or legacy flat format
    fn get_powers(&self) -> Vec<f64> {
        if let Some(ref table) = self.power_thrust_table {
            if !table.power.is_empty() {
                return table.power.clone();
            }
        }
        self.power_curve_powers.clone()
    }

    /// Get thrust coefficient values from either nested power_thrust_table or legacy flat format
    fn get_thrust_coefficients(&self) -> Vec<f64> {
        if let Some(ref table) = self.power_thrust_table {
            if !table.thrust_coefficient.is_empty() {
                return table.thrust_coefficient.clone();
            }
        }
        self.thrust_coefficient_values.clone()
    }

    pub fn power_curve(&self) -> LookupTable {
        LookupTable {
            keys: ndarray::Array1::from(self.get_wind_speeds()),
            values: ndarray::Array1::from(self.get_powers()),
        }
    }

    pub fn thrust_curve(&self) -> LookupTable {
        LookupTable {
            keys: ndarray::Array1::from(self.get_wind_speeds()),
            values: ndarray::Array1::from(self.get_thrust_coefficients()),
        }
    }

    /// Get cosine loss exponent for yaw
    pub fn get_cosine_loss_exponent_yaw(&self) -> f64 {
        self.power_thrust_table
            .as_ref()
            .and_then(|t| t.cosine_loss_exponent_yaw)
            .unwrap_or(1.88)
    }

    /// Get cosine loss exponent for tilt
    pub fn get_cosine_loss_exponent_tilt(&self) -> f64 {
        self.power_thrust_table
            .as_ref()
            .and_then(|t| t.cosine_loss_exponent_tilt)
            .unwrap_or(1.88)
    }

    /// Get reference air density
    pub fn get_ref_air_density(&self) -> f64 {
        self.power_thrust_table
            .as_ref()
            .and_then(|t| t.ref_air_density)
            .unwrap_or(1.225)
    }

    /// Get reference tilt angle
    pub fn get_ref_tilt(&self) -> f64 {
        self.power_thrust_table
            .as_ref()
            .and_then(|t| t.ref_tilt)
            .or(self.ref_tilt)
            .unwrap_or(5.0)
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
    pub fn get_roter_radius(&self) -> f64 {
        self.rotor_diameter / 2.0
    }   

    pub fn get_rotor_area(&self) -> f64 {
        std::f64::consts::PI * (self.get_roter_radius()).powi(2)
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
            floating_tilt_table: None,
            multi_dimensional_cp_ct: None,
            power_thrust_table: None,
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
            floating_tilt_table: None,
            multi_dimensional_cp_ct: None,
            power_thrust_table: None,
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
            floating_tilt_table: None,
            multi_dimensional_cp_ct: None,
            power_thrust_table: None,
        };

        let power = turbine.get_power(10.0, 20.0f64.to_radians());
        assert!((power - 2000.0 * 0.883).abs() < 0.1);
    }
}
