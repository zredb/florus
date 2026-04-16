//! Turbine operation helper functions
//!
//! Shared utilities for turbine operation models

use crate::types::Array2;
use crate::types::Float;

/// Calculate rotor-normal axial induction factor from thrust coefficient
/// Based on classical actuator disk theory: a = (1 - sqrt(1 - Ct)) / 2
pub fn axial_induction_from_ct(ct: &Array2) -> Array2 {
    ct.mapv(|c| (1.0 - (1.0 - c).sqrt()) / 2.0)
}

/// Compute effective power by applying yaw cosine correction
///
/// # Arguments
/// * `power` - Base power values
/// * `yaw_degrees` - Yaw angles in degrees
/// * `yaw_exponent` - Cosine loss exponent for yaw
///
/// Returns power values scaled by cos(yaw)^exponent
pub fn yaw_cosine_correction(power: &Array2, yaw_degrees: &Array2, yaw_exponent: Float) -> Array2 {
    let yaw_rad = yaw_degrees.mapv(|y| y.to_radians());
    power * yaw_rad.mapv(|cosine| cosine.powf(yaw_exponent))
}

/// Convert degrees to radians
pub fn deg_to_rad(degrees: Float) -> Float {
    degrees.to_radians()
}

/// Convert radians to degrees
pub fn rad_to_deg(radians: Float) -> Float {
    radians.to_degrees()
}
