/// Wake model management
///
/// Corresponds to wake.py in Python implementation

pub mod wake_combination;
pub mod wake_deflection;
pub mod wake_turbulence;
pub mod wake_velocity;
pub mod wake_manager;

// Re-export base traits
pub use super::base::{
    BaseModel, CombinationModel, DeflectionModel,
    TurbulenceModel, VelocityModel,
};

// Re-export wake manager and its types
pub use wake_manager::WakeModelManager;
pub use wake_manager::WakeModelStrings;

// Re-export wake deflection models
pub use wake_deflection::{
    EmpiricalGaussVelocityDeflection, GaussVelocityDeflection, JimenezVelocityDeflection,
    NoneVelocityDeflection,
};

// Re-export wake turbulence models
pub use wake_turbulence::{CrespoHernandez, NoneTurbulence};

// Re-export wake velocity models
pub use wake_velocity::{GaussVelocity, JensenVelocity, NoneVelocity};

// Re-export combination models
pub use wake_combination::FLS;
pub use wake_combination::MAX;
pub use wake_combination::SOSFS;
