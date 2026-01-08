/// Wake velocity models
///
/// Corresponds to wake_velocity/ module in Python implementation

pub mod none;
pub mod gauss;
pub mod jensen;

pub use none::NoneVelocity;
pub use gauss::GaussVelocity;
pub use jensen::JensenVelocity;
