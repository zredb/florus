pub mod empirical_gauss;
pub mod gauss;
pub mod jimenez;
/// Wake deflection models
///
/// Corresponds to wake_deflection/ module in Python implementation
pub mod none;

pub use empirical_gauss::EmpiricalGaussVelocityDeflection;
pub use gauss::GaussVelocityDeflection;
pub use jimenez::JimenezVelocityDeflection;
pub use none::NoneVelocityDeflection;
