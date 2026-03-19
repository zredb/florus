pub mod cp_ct_table;
/// Turbine-related modules
///
/// Corresponds to turbine/ module in Python implementation
pub mod operation_models;
pub mod turbine;
pub mod turbine_calculations;
pub mod turbine_library;
pub mod turbine_type;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum TurbineTypeError {
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
    #[error("找不到文件: {0}")]
    InvalidConfiguration(String),
}

pub use turbine_library::TurbineLibrary;

pub use operation_models::{
    AWCTurbine, ControllerDependentTurbine, CosineLossTurbine, MixedOperationTurbine,
    PeakShavingTurbine, SimpleDeratingTurbine, SimpleTurbine, TurbineContext, TurbineParameters,
    UnifiedMomentumTurbine, POWER_SETPOINT_DEFAULT, POWER_SETPOINT_DISABLED,
};
pub use turbine::Turbine;
pub use turbine_calculations::*;
