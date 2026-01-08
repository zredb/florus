/// Wake deflection models
///
/// Corresponds to wake_deflection/ module in Python implementation

pub mod none;
pub mod gauss;
pub mod jimenez;
pub mod empirical_gauss;

pub use none::NoneVelocityDeflection;
pub use gauss::GaussVelocityDeflection;
pub use jimenez::JimenezVelocityDeflection;
pub use empirical_gauss::EmpiricalGaussVelocityDeflection;
