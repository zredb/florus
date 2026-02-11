//! Wake Model Manager
//!
//! Manages all wake models (velocity, deflection, turbulence, combination)

use crate::types::NumericDict;
use crate::core::wake::{
    CombinationModel, DeflectionModel, TurbulenceModel, VelocityModel,
};
use std::collections::HashMap;

/// Wake model strings for identification
#[derive(Debug, Clone)]
pub struct WakeModelStrings {
    pub velocity_model: String,
    pub deflection_model: String,
    pub combination_model: String,
    pub turbulence_model: String,
}

/// Wake Model Manager
pub struct WakeModelManager {
    pub velocity_model: Box<dyn VelocityModel>,
    pub deflection_model: Box<dyn DeflectionModel>,
    pub turbulence_model: Box<dyn TurbulenceModel>,
    pub combination_model: Box<dyn CombinationModel>,
    pub model_params: HashMap<String, NumericDict>,
    pub turbine_type_params: HashMap<String, NumericDict>,
    pub enable_secondary_steering: bool,
    pub enable_yaw_added_recovery: bool,
    pub use_parallel_calc: bool,
}

impl std::fmt::Debug for WakeModelManager {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("WakeModelManager")
            .field("velocity_model", &"Box<dyn VelocityModel>")
            .field("deflection_model", &"Box<dyn DeflectionModel>")
            .field("turbulence_model", &"Box<dyn TurbulenceModel>")
            .field("combination_model", &"Box<dyn CombinationModel>")
            .field("model_params", &self.model_params.keys().collect::<Vec<_>>())
            .field("enable_secondary_steering", &self.enable_secondary_steering)
            .field("enable_yaw_added_recovery", &self.enable_yaw_added_recovery)
            .field("use_parallel_calc", &self.use_parallel_calc)
            .finish()
    }
}

impl Clone for WakeModelManager {
    fn clone(&self) -> Self {
        Self {
            velocity_model: Box::new(crate::core::wake::GaussVelocity::new(0.38, 0.004, 0.5)),
            deflection_model: Box::new(crate::core::wake::GaussVelocityDeflection::new(0.0, 0.0, 0.58, 0.077, 1.0)),
            turbulence_model: Box::new(crate::core::wake::CrespoHernandez::new(0.9, 0.9)),
            combination_model: Box::new(crate::core::wake::FLS),
            model_params: self.model_params.clone(),
            turbine_type_params: self.turbine_type_params.clone(),
            enable_secondary_steering: self.enable_secondary_steering,
            enable_yaw_added_recovery: self.enable_yaw_added_recovery,
            use_parallel_calc: self.use_parallel_calc,
        }
    }
}

impl WakeModelManager {
    /// Create a new WakeModelManager from model strings
    pub fn new(
        model_strings: WakeModelStrings,
        model_params: HashMap<String, NumericDict>,
        turbine_type_params: HashMap<String, NumericDict>,
        _turbine_type_templates: HashMap<String, NumericDict>,
        enable_secondary_steering: bool,
        enable_yaw_added_recovery: bool,
        use_parallel_calc: bool,
    ) -> anyhow::Result<Self> {
        let velocity_model = Self::create_velocity_model(
            &model_strings.velocity_model,
            &model_params,
            &turbine_type_params,
        )?;
        
        let deflection_model = Self::create_deflection_model(
            &model_strings.deflection_model,
            &model_params,
        )?;
        
        let turbulence_model = Self::create_turbulence_model(
            &model_strings.turbulence_model,
            &model_params,
        )?;
        
        let combination_model = Self::create_combination_model(
            &model_strings.combination_model,
        )?;

        Ok(Self {
            velocity_model,
            deflection_model,
            turbulence_model,
            combination_model,
            model_params,
            turbine_type_params,
            enable_secondary_steering,
            enable_yaw_added_recovery,
            use_parallel_calc,
        })
    }

    /// Create velocity model from string identifier
    fn create_velocity_model(
        model_name: &str,
        model_params: &HashMap<String, NumericDict>,
        _turbine_type_params: &HashMap<String, NumericDict>,
    ) -> anyhow::Result<Box<dyn VelocityModel>> {
        match model_name.to_lowercase().as_str() {
            "gauss" | "gaussian" => {
                // model_params contains nested dict like: {"gauss": {"ka": 0.38, "kb": 0.004, ...}}
                let gauss_params = model_params.get("gauss")
                    .cloned()
                    .unwrap_or_else(|| NumericDict {
                        data: std::collections::HashMap::new()
                    });

                let ka = gauss_params.data.get("ka")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.38);
                let kb = gauss_params.data.get("kb")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.004);
                let initial_wake_width = gauss_params.data.get("alpha")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.5);

                Ok(Box::new(
                    crate::core::wake::GaussVelocity::new(ka, kb, initial_wake_width)
                ))
            },
            "jensen" => {
                let jensen_params = model_params.get("jensen")
                    .cloned()
                    .unwrap_or_else(|| NumericDict {
                        data: std::collections::HashMap::new()
                    });

                let we = jensen_params.data.get("we")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.05);
                Ok(Box::new(
                    crate::core::wake::JensenVelocity::new(we, 0.0)
                ))
            },
            "turbopark" => Ok(Box::new(
                crate::core::wake_velocity::TurbOParkVelocityDeficit::default()
            )),
            "turboparkgauss" => Ok(Box::new(
                crate::core::wake_velocity::TurbOParkGaussVelocityDeficit::default()
            )),
            "gauss_legacy" | "gch" => Ok(Box::new(
                crate::core::wake_velocity::GaussLegacyVelocityDeficit::default()
            )),
            "empirical_gauss" | "empirical-gauss" => Ok(Box::new(
                crate::core::wake_velocity::EmpiricalGaussVelocityDeficit::default()
            )),
            "cumulative_gauss_curl" => Ok(Box::new(
                crate::core::wake_velocity::CumulativeCurlVelocityDeficit::default()
            )),
            "none" => Ok(Box::new(crate::core::wake::NoneVelocity::new())),
            _ => Err(anyhow::anyhow!(
                "Velocity model '{}' not implemented",
                model_name
            )),
        }
    }

    /// Create deflection model from string identifier

    /// Create deflection model from string identifier
    fn create_deflection_model(
        model_name: &str,
        model_params: &HashMap<String, NumericDict>,
    ) -> anyhow::Result<Box<dyn DeflectionModel>> {
        match model_name.to_lowercase().as_str() {
            "gauss" | "gaussian" => {
                // model_params contains nested dict like: {"gauss": {"ad": 0.0, "alpha": 0.58, ...}}
                let gauss_params = model_params.get("gauss")
                    .cloned()
                    .unwrap_or_else(|| NumericDict {
                        data: std::collections::HashMap::new()
                    });

                let ad = gauss_params.data.get("ad")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.0);
                let bd = gauss_params.data.get("bd")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.0);
                let alpha = gauss_params.data.get("alpha")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.58);
                let beta = gauss_params.data.get("beta")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.077);
                let dm = gauss_params.data.get("dm")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(1.0);

                Ok(Box::new(
                    crate::core::wake::GaussVelocityDeflection::new(ad, bd, alpha, beta, dm)
                ))
            },
            "jimenez" => {
                let jimenez_params = model_params.get("jimenez")
                    .cloned()
                    .unwrap_or_else(|| NumericDict {
                        data: std::collections::HashMap::new()
                    });

                let ad = jimenez_params.data.get("ad")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.0);
                let bd = jimenez_params.data.get("bd")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.0);
                let _kd = jimenez_params.data.get("kd")
                    .and_then(|v| match v { crate::types::ConfigValue::Float(f) => Some(*f), _ => None })
                    .unwrap_or(0.05);
                Ok(Box::new(
                    crate::core::wake::JimenezVelocityDeflection::new(ad, bd)
                ))
            },
            "empirical_gauss" | "empirical-gauss" => Ok(Box::new(
                crate::core::wake::EmpiricalGaussVelocityDeflection::new(0.01, 0.05, 0.5)
            )),
            "none" => Ok(Box::new(crate::core::wake::NoneVelocityDeflection::new())),
            _ => Err(anyhow::anyhow!(
                "Deflection model '{}' not implemented",
                model_name
            )),
        }
    }

    /// Create turbulence model from string identifier
    fn create_turbulence_model(
        model_name: &str,
        _model_params: &HashMap<String, NumericDict>,
    ) -> anyhow::Result<Box<dyn TurbulenceModel>> {
        match model_name.to_lowercase().as_str() {
            "crespo_hernandez" | "crespo-hernandez" => Ok(Box::new(
                crate::core::wake::CrespoHernandez::new(0.9, 0.9)
            )),
            "none" => Ok(Box::new(crate::core::wake::NoneTurbulence::new())),
            _ => Err(anyhow::anyhow!(
                "Turbulence model '{}' not implemented",
                model_name
            )),
        }
    }

    /// Create combination model from string identifier
    fn create_combination_model(
        model_name: &str,
    ) -> anyhow::Result<Box<dyn CombinationModel>> {
        match model_name.to_lowercase().as_str() {
            "fls" | "freestream_linear_superposition" => Ok(Box::new(
                crate::core::wake::FLS
            )),
            "max" | "maximum" => Ok(Box::new(crate::core::wake::MAX)),
            "sosfs" | "sos" => Ok(Box::new(crate::core::wake::SOSFS)),
            _ => Err(anyhow::anyhow!(
                "Combination model '{}' not implemented",
                model_name
            )),
        }
    }

    /// Get default model strings for common configurations
    pub fn default_gauss() -> WakeModelStrings {
        WakeModelStrings {
            velocity_model: "gauss".to_string(),
            deflection_model: "gauss".to_string(),
            combination_model: "fls".to_string(),
            turbulence_model: "crespo_hernandez".to_string(),
        }
    }

    pub fn default_jensen() -> WakeModelStrings {
        WakeModelStrings {
            velocity_model: "jensen".to_string(),
            deflection_model: "jimenez".to_string(),
            combination_model: "fls".to_string(),
            turbulence_model: "none".to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wake_model_manager_creation() {
        let model_strings = WakeModelManager::default_gauss();
        let manager = WakeModelManager::new(
            model_strings,
            HashMap::new(),
            HashMap::new(),
            HashMap::new(),
            false,
            false,
            false,
        );
        assert!(manager.is_ok());
    }

    #[test]
    fn test_default_gauss_models() {
        let strings = WakeModelManager::default_gauss();
        assert_eq!(strings.velocity_model, "gauss");
        assert_eq!(strings.deflection_model, "gauss");
        assert_eq!(strings.combination_model, "fls");
        assert_eq!(strings.turbulence_model, "crespo_hernandez");
    }

    #[test]
    fn test_default_jensen_models() {
        let strings = WakeModelManager::default_jensen();
        assert_eq!(strings.velocity_model, "jensen");
        assert_eq!(strings.deflection_model, "jimenez");
    }
}
