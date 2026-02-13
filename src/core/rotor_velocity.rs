/// Rotor velocity calculations
///
/// Corresponds to rotor_velocity.py in the Python implementation
use crate::types::{Float, Array2, Array3, Array4};
use crate::utilities::cosd;
use ndarray::{s, Axis};

/// Apply yaw cosine correction to rotor effective velocities
/// 
/// # Arguments
/// * `cosine_loss_exponent_yaw` - Yaw cosine loss exponent
/// * `yaw_angles` - Yaw angles in degrees (findex, n_turbines)
/// * `rotor_effective_velocities` - Effective velocities at rotor (findex, n_turbines)
pub fn rotor_velocity_yaw_cosine_correction(
    cosine_loss_exponent_yaw: Float,
    yaw_angles: &Array2,
    rotor_effective_velocities: &Array2,
) -> Array2 {
    let pw = cosine_loss_exponent_yaw / 3.0;
    
    let mut result = rotor_effective_velocities.clone();
    for ((i, j), vel) in result.indexed_iter_mut() {
        let yaw = yaw_angles[[i, j]];
        *vel *= cosd(yaw).powf(pw);
    }
    
    result
}

/// Apply tilt cosine correction to rotor effective velocities
/// 
/// # Arguments
/// * `tilt_angles` - Tilt angles in degrees
/// * `ref_tilt` - Reference tilt angles
/// * `cosine_loss_exponent_tilt` - Tilt cosine loss exponent
/// * `correct_cp_ct_for_tilt` - Whether to correct for tilt
/// * `rotor_effective_velocities` - Effective velocities at rotor
pub fn rotor_velocity_tilt_cosine_correction(
    tilt_angles: &Array2,
    ref_tilt: &Array2,
    cosine_loss_exponent_tilt: Float,
    correct_cp_ct_for_tilt: &ndarray::Array2<bool>,
    rotor_effective_velocities: &Array2,
) -> Array2 {
    let mut result = rotor_effective_velocities.clone();
    let exponent = cosine_loss_exponent_tilt / 3.0;
    
    for ((i, j), vel) in result.indexed_iter_mut() {
        if correct_cp_ct_for_tilt[[i, j]] {
            let tilt = tilt_angles[[i, j]];
            let ref_t = ref_tilt[[i, j]];
            *vel *= (cosd(tilt) / cosd(ref_t)).powf(exponent);
        }
    }
    
    result
}

/// Calculate simple arithmetic mean along specified axis
pub fn simple_mean(array: &Array3, axis: usize) -> Array2 {
    array.mean_axis(Axis(axis)).unwrap()
}

/// Calculate cubic mean along specified axis
/// 
/// Cubic mean = (mean(x³))^(1/3)
pub fn cubic_mean(array: &Array3, axis: usize) -> Array2 {
    let cubed = array.mapv(|x| x.powi(3));
    let mean = cubed.mean_axis(Axis(axis)).unwrap();
    mean.mapv(|x| x.cbrt())
}

/// Calculate simple mean with cubature weights
pub fn simple_cubature(array: &Array4, cubature_weights: &Array2, axis: usize) -> Array3 {
    // Flatten weights and normalize
    let weights_flat = cubature_weights.iter().cloned().collect::<Vec<_>>();
    let n = weights_flat.len() as Float;
    let sum: Float = weights_flat.iter().sum();
    let normalized_weights: Vec<Float> = weights_flat.iter().map(|w| w * n / sum).collect();
    
    // Apply weights and calculate mean
    let mut weighted = array.clone();
    for (idx, w) in normalized_weights.iter().enumerate() {
        weighted.slice_mut(s![.., .., idx, ..]).mapv_inplace(|x| x * w);
    }
    
    weighted.mean_axis(Axis(axis)).unwrap()
}

/// Calculate cubic mean with cubature weights
pub fn cubic_cubature(array: &Array4, cubature_weights: &Array2, axis: usize) -> Array3 {
    // Flatten weights and normalize
    let weights_flat = cubature_weights.iter().cloned().collect::<Vec<_>>();
    let n = weights_flat.len() as Float;
    let sum: Float = weights_flat.iter().sum();
    let normalized_weights: Vec<Float> = weights_flat.iter().map(|w| w * n / sum).collect();
    
    // Apply weights to cubed values
    let cubed = array.mapv(|x| x.powi(3));
    let mut weighted = cubed.clone();
    for (idx, w) in normalized_weights.iter().enumerate() {
        weighted.slice_mut(s![.., .., idx, ..]).mapv_inplace(|x| x * w);
    }
    
    let mean = weighted.mean_axis(Axis(axis)).unwrap();
    mean.mapv(|x| x.cbrt())
}

/// Averaging method enum
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AveragingMethod {
    SimpleMean,
    CubicMean,
    SimpleCubature,
    CubicCubature,
}

impl AveragingMethod {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "simple-mean" => Some(Self::SimpleMean),
            "cubic-mean" => Some(Self::CubicMean),
            "simple-cubature" => Some(Self::SimpleCubature),
            "cubic-cubature" => Some(Self::CubicCubature),
            _ => None,
        }
    }
}

/// Calculate average velocity using specified method
/// 
/// # Arguments
/// * `velocities` - Velocity array (findex, n_turbines, n_points, n_wd)
/// * `method` - Averaging method to use
/// * `cubature_weights` - Optional cubature weights for cubature methods
pub fn average_velocity(
    velocities: &Array4,
    method: AveragingMethod,
    cubature_weights: Option<&Array2>,
) -> crate::Result<Array3> {
    match method {
        AveragingMethod::SimpleMean => {
            // Average over axis 2 (n_points)
            Ok(velocities.mean_axis(Axis(2)).unwrap())
        }
        AveragingMethod::CubicMean => {
            let cubed = velocities.mapv(|x| x.powi(3));
            let mean = cubed.mean_axis(Axis(2)).unwrap();
            Ok(mean.mapv(|x| x.cbrt()))
        }
        AveragingMethod::SimpleCubature => {
            let weights = cubature_weights
                .ok_or_else(|| anyhow::anyhow!("Cubature weights required for simple-cubature"))?;
            Ok(simple_cubature(velocities, weights, 2))
        }
        AveragingMethod::CubicCubature => {
            let weights = cubature_weights
                .ok_or_else(|| anyhow::anyhow!("Cubature weights required for cubic-cubature"))?;
            Ok(cubic_cubature(velocities, weights, 2))
        }
    }
}

/// Calculate rotor effective velocity
/// 
/// This is the main function to compute the effective velocity at each turbine rotor
/// averaged over the rotor disk using the specified method.
pub fn rotor_effective_velocity(
    velocities: &Array4,
    method: AveragingMethod,
    cubature_weights: Option<&Array2>,
    yaw_angles: Option<&Array2>,
    tilt_angles: Option<&Array2>,
    ref_tilt: Option<&Array2>,
    cosine_loss_exponent_yaw: Option<Float>,
    cosine_loss_exponent_tilt: Option<Float>,
    correct_cp_ct_for_tilt: Option<&ndarray::Array2<bool>>,
) -> crate::Result<Array3> {
    // Calculate base average velocity
    let mut avg_vel = average_velocity(velocities, method, cubature_weights)?;
    
    // Apply yaw correction if parameters provided
    if let (Some(yaw), Some(exp_yaw)) = (yaw_angles, cosine_loss_exponent_yaw) {
        // Need to apply correction per findex
        for i in 0..avg_vel.shape()[0] {
            for j in 0..avg_vel.shape()[1] {
                let pw = exp_yaw / 3.0;
                let yaw_deg = yaw[[i, j]];
                avg_vel[[i, j, 0]] *= cosd(yaw_deg).powf(pw);
            }
        }
    }
    
    // Apply tilt correction if parameters provided
    if let (Some(tilt), Some(ref_t), Some(exp_tilt), Some(correct)) = 
        (tilt_angles, ref_tilt, cosine_loss_exponent_tilt, correct_cp_ct_for_tilt) {
        for i in 0..avg_vel.shape()[0] {
            for j in 0..avg_vel.shape()[1] {
                if correct[[i, j]] {
                    let exponent = exp_tilt / 3.0;
                    let tilt_deg = tilt[[i, j]];
                    let ref_tilt_deg = ref_t[[i, j]];
                    avg_vel[[i, j, 0]] *= (cosd(tilt_deg) / cosd(ref_tilt_deg)).powf(exponent);
                }
            }
        }
    }
    
    Ok(avg_vel)
}

/// Apply air density correction to velocities
/// 
/// Produces equivalent velocities at the reference air density
/// 
/// # Arguments
/// * `velocities` - Input velocities
/// * `air_density` - Current air density
/// * `ref_air_density` - Reference air density
pub fn rotor_velocity_air_density_correction(
    velocities: &Array2,
    air_density: Float,
    ref_air_density: Float,
) -> Array2 {
    let density_ratio = (air_density / ref_air_density).powf(1.0 / 3.0);
    velocities.mapv(|x| x * density_ratio)
}

/// Compute tilt angles for floating turbines
///
/// For floating turbines, the tilt angle may change based on wind speed.
/// This function interpolates the tilt angle based on rotor effective velocity.
pub fn compute_tilt_angles_for_floating_turbines(
    _tilt_angles: &Array2,
    _tilt_interp: Option<&Vec<Option<()>>>, // Placeholder for interpolator
    _rotor_effective_velocities: &Array2,
) -> Array2 {
    // Simplified implementation - would need proper interpolation in full version
    _tilt_angles.clone()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_simple_mean() {
        let data = Array::from_shape_vec(
            (2, 2, 3),
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0],
        ).unwrap();
        
        let result = simple_mean(&data, 2);
        assert_eq!(result.shape(), &[2, 2]);
        assert_relative_eq!(result[[0, 0]], 2.0);
        assert_relative_eq!(result[[0, 1]], 5.0);
    }
    
    #[test]
    fn test_cubic_mean() {
        let data = Array::from_shape_vec(
            (1, 1, 3),
            vec![1.0, 2.0, 3.0],
        ).unwrap();
        
        let result = cubic_mean(&data, 2);
        let expected = ((1.0_f64.powi(3) + 2.0_f64.powi(3) + 3.0_f64.powi(3)) / 3.0).cbrt();
        assert_relative_eq!(result[[0, 0]], expected, epsilon = 1e-10);
    }
    
    #[test]
    fn test_averaging_method_from_str() {
        assert_eq!(
            AveragingMethod::from_str("simple-mean"),
            Some(AveragingMethod::SimpleMean)
        );
        assert_eq!(
            AveragingMethod::from_str("cubic-mean"),
            Some(AveragingMethod::CubicMean)
        );
        assert_eq!(
            AveragingMethod::from_str("invalid"),
            None
        );
    }
}
