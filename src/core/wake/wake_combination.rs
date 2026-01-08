/// Wake combination models
///
/// Corresponds to wake_combination/ module in Python implementation

pub mod fls;
pub mod max;
pub mod sosfs;

pub use fls::FLS;
pub use max::MAX;
pub use sosfs::SOSFS;

#[derive(Debug, Clone)]
pub struct WakeModelStrings {
    pub velocity_model: String,
    pub deflection_model: String,
    pub combination_model: String,
    pub turbulence_model: String,
}
