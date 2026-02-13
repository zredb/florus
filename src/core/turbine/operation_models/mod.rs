/// Turbine operation models
///
/// Modular turbine operation models organized by functionality:
/// - base.rs: Base types and trait definitions
/// - helpers.rs: Shared utility functions
/// - simple.rs: SimpleTurbine model
/// - cosine_loss.rs: CosineLossTurbine model
/// - simple_derating.rs: SimpleDeratingTurbine model
/// - awc.rs: AWCTurbine model (placeholder)
/// - peak_shaving.rs: PeakShavingTurbine model (placeholder)
/// - mixed.rs: MixedOperationTurbine model
/// - unified_momentum.rs: UnifiedMomentumTurbine model (full implementation)
/// - controller_dependent.rs: ControllerDependentTurbine model (framework for custom controls)

pub mod base;
pub mod helpers;
pub mod simple;
pub mod cosine_loss;
pub mod simple_derating;
pub mod awc;
pub mod peak_shaving;
pub mod mixed;
pub mod unified_momentum;
pub mod controller_dependent;

pub use base::*;
pub use helpers::*;

// Re-export main types
pub use base::{OperationModel, TurbineParameters, TurbineContext, POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED};

// Re-export operation models
pub use simple::SimpleTurbine;
pub use cosine_loss::CosineLossTurbine;
pub use simple_derating::SimpleDeratingTurbine;
pub use awc::AWCTurbine;
pub use peak_shaving::PeakShavingTurbine;
pub use mixed::MixedOperationTurbine;
pub use unified_momentum::UnifiedMomentumTurbine;
pub use controller_dependent::ControllerDependentTurbine;
