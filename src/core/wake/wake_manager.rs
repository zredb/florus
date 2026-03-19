//! Wake Model Manager
//!
//! Manages all wake models (velocity, deflection, turbulence, combination)

use crate::core::wake::{CombinationModel, DeflectionModel, TurbulenceModel, VelocityModel};
use crate::floris_config::{
    CombinationModelConfig, DeflectionModelConfig, TurbulenceModelConfig, VelocityModelConfig,
    WakeConfig,
};

/// Wake model manager created from WakeConfig
///
/// This version is refactored to work directly with WakeConfig enums
/// instead of string-based model selection.
pub struct WakeModelManager {
    pub velocity_model: Box<dyn VelocityModel>,
    pub deflection_model: Box<dyn DeflectionModel>,
    pub turbulence_model: Box<dyn TurbulenceModel>,
    pub combination_model: Box<dyn CombinationModel>,
    pub enable_secondary_steering: bool,
    pub enable_yaw_added_recovery: bool,
    pub enable_active_wake_mixing: bool,
    pub enable_transverse_velocities: bool,
    pub enable_wake_mixing: bool,
    pub use_parallel_calc: bool,
}

impl std::fmt::Debug for WakeModelManager {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("WakeModelManager")
            .field("velocity_model", &"Box<dyn VelocityModel>")
            .field("deflection_model", &"Box<dyn DeflectionModel>")
            .field("turbulence_model", &"Box<dyn TurbulenceModel>")
            .field("combination_model", &"Box<dyn CombinationModel>")
            .field("enable_secondary_steering", &self.enable_secondary_steering)
            .field("enable_yaw_added_recovery", &self.enable_yaw_added_recovery)
            .field("enable_active_wake_mixing", &self.enable_active_wake_mixing)
            .field(
                "enable_transverse_velocities",
                &self.enable_transverse_velocities,
            )
            .field("enable_wake_mixing", &self.enable_wake_mixing)
            .field("use_parallel_calc", &self.use_parallel_calc)
            .finish()
    }
}

impl Clone for WakeModelManager {
    fn clone(&self) -> Self {
        Self {
            velocity_model: Box::new(crate::core::wake_velocity::gauss::GaussVelocity::new(
                0.58, 0.077, 0.38, 0.004,
            )),
            deflection_model: Box::new(
                crate::core::wake_deflection::gauss::GaussVelocityDeflection::new(
                    0.05, 0.0, 0.58, 0.077, 1.0,
                ),
            ),
            turbulence_model: Box::new(
                crate::core::wake_turbulence::crespo_hernandez::CrespoHernandez::new(0.9, 0.9),
            ),
            combination_model: Box::new(crate::core::wake_combination::FLS),
            enable_secondary_steering: self.enable_secondary_steering,
            enable_yaw_added_recovery: self.enable_yaw_added_recovery,
            enable_active_wake_mixing: self.enable_active_wake_mixing,
            enable_transverse_velocities: self.enable_transverse_velocities,
            enable_wake_mixing: self.enable_wake_mixing,
            use_parallel_calc: self.use_parallel_calc,
        }
    }
}

impl WakeModelManager {
    /// Create a new WakeModelManager from WakeConfig
    pub fn from_config(config: &WakeConfig) -> anyhow::Result<Self> {
        let velocity_model = Self::create_velocity_model(&config.model_strings.velocity_model)?;
        let deflection_model =
            Self::create_deflection_model(&config.model_strings.deflection_model)?;
        let turbulence_model =
            Self::create_turbulence_model(&config.model_strings.turbulence_model)?;
        let combination_model =
            Self::create_combination_model(&config.model_strings.combination_model)?;

        Ok(Self {
            velocity_model,
            deflection_model,
            turbulence_model,
            combination_model,
            enable_secondary_steering: config.enable_secondary_steering,
            enable_yaw_added_recovery: config.enable_yaw_added_recovery,
            enable_active_wake_mixing: config.enable_active_wake_mixing,
            enable_transverse_velocities: config.enable_transverse_velocities,
            enable_wake_mixing: config.enable_wake_mixing,
            use_parallel_calc: config.use_parallel_calc,
        })
    }

    /// Create velocity model from VelocityModelConfig enum
    fn create_velocity_model(
        config: &VelocityModelConfig,
    ) -> anyhow::Result<Box<dyn VelocityModel>> {
        match config {
            VelocityModelConfig::Gauss { alpha, beta, ka, kb } => {
                Ok(Box::new(
                    crate::core::wake_velocity::gauss::GaussVelocity::new(*alpha, *beta, *ka, *kb)
                ))
            }
            VelocityModelConfig::Jensen { we } => {
                // Jensen uses we (wake expansion rate) as kd
                Ok(Box::new(
                    crate::core::wake_velocity::jensen::JensenVelocity::new(*we, 0.0)
                ))
            }
            VelocityModelConfig::Turbopark { kstar, cstar } => {
                // TurbOPark parameters: a and sigma_max_rel
                // kstar controls wake width, cstar controls decay
                Ok(Box::new(
                    crate::core::wake_velocity::turbopark::TurbOParkVelocityDeficit::new(*kstar, *cstar)
                ))
            }
            VelocityModelConfig::TurboparkGauss { kstar, cstar: _ } => {
                Ok(Box::new(
                    crate::core::wake_velocity::turboparkgauss::TurbOParkGaussVelocityDeficit::new(*kstar, true)
                ))
            }
            VelocityModelConfig::CC {
                a_s,
                b_s,
                c_s1,
                c_s2,
                a_f,
                b_f,
                c_f,
                alpha_mod,
            } => {
                Ok(Box::new(
                    crate::core::wake_velocity::cumulative_gauss_curl::CumulativeCurlVelocityDeficit::new(
                        *alpha_mod, *a_s, *b_s, *c_s1, *c_s2,
                        *a_f, *b_f, *c_f,
                    )
                ))
            }
            VelocityModelConfig::CumulativeGaussCurl {
                a_s,
                b_s,
                c_s1,
                c_s2,
                a_f,
                b_f,
                c_f,
                alpha_mod,
            } => {
                Ok(Box::new(
                    crate::core::wake_velocity::cumulative_gauss_curl::CumulativeCurlVelocityDeficit::new(
                        *alpha_mod, *a_s, *b_s, *c_s1, *c_s2,
                        *a_f, *b_f, *c_f,
                    )
                ))
            }
            VelocityModelConfig::EmpiricalGauss { ad, bd: _, alpha, beta } => {
                Ok(Box::new(
                    crate::core::wake_velocity::empirical_gauss::EmpiricalGaussVelocityDeficit::new(*alpha, *beta, *ad)
                ))
            }
        }
    }

    /// Create deflection model from DeflectionModelConfig enum
    fn create_deflection_model(
        config: &DeflectionModelConfig,
    ) -> anyhow::Result<Box<dyn DeflectionModel>> {
        match config {
            DeflectionModelConfig::Gauss {
                alpha,
                beta,
                ad,
                bd: _,
                dm,
                ka,
                kb: _,
            } => {
                // Gauss deflection uses kd for wake expansion
                // Use ka as default kd if not specified separately
                let kd = *ka;
                Ok(Box::new(
                    crate::core::wake_deflection::gauss::GaussVelocityDeflection::new(kd, *ad, *alpha, *beta, *dm)
                ))
            }
            DeflectionModelConfig::Jimenez { ad, bd: _, kd } => {
                Ok(Box::new(
                    crate::core::wake_deflection::jimenez::JimenezVelocityDeflection::new(*kd, *ad)
                ))
            }
            DeflectionModelConfig::EmpiricalGauss { ad, bd, kd } => {
                Ok(Box::new(
                    crate::core::wake_deflection::empirical_gauss::EmpiricalGaussVelocityDeflection::new(*ad, *bd, *kd)
                ))
            }
        }
    }

    /// Create turbulence model from TurbulenceModelConfig enum
    fn create_turbulence_model(
        config: &TurbulenceModelConfig,
    ) -> anyhow::Result<Box<dyn TurbulenceModel>> {
        match config {
            TurbulenceModelConfig::CrespoHernandez {
                initial,
                constant,
                ai: _,
                downstream: _,
            } => Ok(Box::new(
                crate::core::wake_turbulence::crespo_hernandez::CrespoHernandez::new(
                    *initial, *constant,
                ),
            )),
            TurbulenceModelConfig::WakeInducedMixing { .. } => {
                // WakeInducedMixing not fully implemented yet, use CrespoHernandez as fallback
                Ok(Box::new(
                    crate::core::wake_turbulence::crespo_hernandez::CrespoHernandez::new(0.9, 0.9),
                ))
            }
            TurbulenceModelConfig::None => Ok(Box::new(
                crate::core::wake_turbulence::none::NoneTurbulence::new(),
            )),
        }
    }

    /// Create combination model from CombinationModelConfig enum
    fn create_combination_model(
        config: &CombinationModelConfig,
    ) -> anyhow::Result<Box<dyn CombinationModel>> {
        match config {
            CombinationModelConfig::FLS => Ok(Box::new(crate::core::wake_combination::FLS)),
            CombinationModelConfig::SOSFS => Ok(Box::new(crate::core::wake_combination::SOSFS)),
            CombinationModelConfig::Max => Ok(Box::new(crate::core::wake_combination::MAX)),
        }
    }

    /// Get default WakeModelManager (Gaussian models with FLS combination)
    pub fn default_gauss() -> anyhow::Result<Self> {
        let config = WakeConfig::default();
        Self::from_config(&config)
    }

    /// Get default WakeModelManager (Jensen models with FLS combination)
    pub fn default_jensen() -> anyhow::Result<Self> {
        let mut config = WakeConfig::default();
        config.model_strings.velocity_model = VelocityModelConfig::Jensen { we: 0.05 };
        config.model_strings.deflection_model = DeflectionModelConfig::Jimenez {
            ad: 0.0,
            bd: 0.0,
            kd: 0.05,
        };
        config.model_strings.turbulence_model = TurbulenceModelConfig::None;

        Self::from_config(&config)
    }
   

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wake_model_manager_from_config() {
        let config = WakeConfig::default();
        let manager = WakeModelManager::from_config(&config);
        assert!(manager.is_ok());
    }

    #[test]
    fn test_default_gauss() {
        let manager = WakeModelManager::default_gauss();
        assert!(manager.is_ok());
    }

    #[test]
    fn test_default_jensen() {
        let manager = WakeModelManager::default_jensen();
        assert!(manager.is_ok());
    }
}
