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
            velocity_model: Box::new(crate::core::wake::GaussVelocity::new(0.1, 0.05, 0.5)),
            deflection_model: Box::new(crate::core::wake::GaussVelocityDeflection::new(0.01, 0.05)),
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
        _model_params: &HashMap<String, NumericDict>,
        _turbine_type_params: &HashMap<String, NumericDict>,
    ) -> anyhow::Result<Box<dyn VelocityModel>> {
        match model_name.to_lowercase().as_str() {
            "gauss" | "gaussian" => Ok(Box::new(
                crate::core::wake::GaussVelocity::new(0.1, 0.05, 0.5)
            )),
            "jensen" => Ok(Box::new(
                crate::core::wake::JensenVelocity::new(0.1, 0.0)
            )),
            "none" => Ok(Box::new(crate::core::wake::NoneVelocity::new())),
            _ => Err(anyhow::anyhow!(
                "Velocity model '{}' not implemented",
                model_name
            )),
        }
    }

    /// Create deflection model from string identifier
    fn create_deflection_model(
        model_name: &str,
        _model_params: &HashMap<String, NumericDict>,
    ) -> anyhow::Result<Box<dyn DeflectionModel>> {
        match model_name.to_lowercase().as_str() {
            "gauss" | "gaussian" => Ok(Box::new(
                crate::core::wake::GaussVelocityDeflection::new(0.01, 0.05)
            )),
            "jimenez" => Ok(Box::new(
                crate::core::wake::JimenezVelocityDeflection::new(0.01, 0.05)
            )),
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
