/// Rotor velocity calculations
///
/// Corresponds to rotor_velocity.py in the Python implementation
use crate::types::{Array1, Array2, Array3, Array4, ArrayD, Float};
use crate::utilities::cosd;
use ndarray::{s, Axis};
use std::fmt;

/// Helper function to convert ArrayD to Array2
fn convert_to_array2(arr: ArrayD) -> Array2 {
    // Ensure the array is 2D
    assert_eq!(arr.ndim(), 2, "Expected 2D array, got {}D", arr.ndim());
    let shape = arr.shape().to_vec();
    arr.into_shape((shape[0], shape[1])).unwrap()
}
/// Apply yaw cosine correction to rotor effective velocities
///
/// # Arguments
/// * `cosine_loss_exponent_yaw` - Yaw cosine loss exponent
/// * `yaw_angles` - Yaw angles in degrees (findex, n_turbines);
/// * `rotor_effective_velocities` - Effective velocities at rotor (findex, n_turbines)
/// * Array2第一维是风况索引 (n_findex), 第二维是风机索引 (n_turbines), 单位是度（degrees）
pub fn rotor_velocity_yaw_cosine_correction(
    cosine_loss_exponent_yaw: Float,
    yaw_angles: &Array2,
    rotor_effective_velocities: &Array2,
) -> Array2 {
    // 已验证 2026-05-13
    let pw = cosine_loss_exponent_yaw / 3.0;
    // Compute cosine correction factors for all elements at once
    let correction_factors = yaw_angles.mapv(|yaw| cosd(yaw).powf(pw));

    // Element-wise multiplication
    rotor_effective_velocities * &correction_factors
}

/// Apply tilt cosine correction to rotor effective velocities
///
/// # Arguments
/// * `tilt_angles` - Tilt angles in degrees
/// * `ref_tilt` - Reference tilt angles
/// * `cosine_loss_exponent_tilt` - Tilt cosine loss exponent
/// * `correct_cp_ct_for_tilt` - Whether to correct for tilt
/// * `rotor_effective_velocities` - Effective velocities at rotor
/// * Array2第一维是风况索引 (n_findex), 第二维是风机索引 (n_turbines), 单位是度（degrees）
pub fn rotor_velocity_tilt_cosine_correction(
    tilt_angles: &Array2,
    ref_tilt: &Array2,
    cosine_loss_exponent_tilt: Float,
    // tilt_interp: NDArrayObject, 暂时没有用到, 先省略
    correct_cp_ct_for_tilt: &ndarray::Array2<bool>,
    rotor_effective_velocities: &Array2,
) -> Array2 {
    // 已验证 2026-05-13
    let exponent = cosine_loss_exponent_tilt / 3.0;

    // Match Floris Python behavior:
    // When correct_cp_ct_for_tilt is False, use the original tilt_angles
    // When True, the tilt_angles have already been updated (e.g., for floating turbines)
    // In both cases, we apply the correction using the appropriate tilt angle

    // Precompute all correction factors using vectorized operations
    let tilt_cos = tilt_angles.mapv(|t| cosd(t));
    let ref_tilt_cos = ref_tilt.mapv(|t| cosd(t));
    let correction = (&tilt_cos / &ref_tilt_cos).mapv(|x| x.powf(exponent));

    // Apply corrections element-wise (matching Python's behavior)
    // Note: Python applies correction to ALL elements, but uses different tilt values
    // based on correct_cp_ct_for_tilt. Since tilt_angles already contains the correct
    // values (either original or computed for floating), we apply correction everywhere.
    rotor_effective_velocities * &correction
}

#[derive(Debug)]
pub enum AvgVelocityError {
    MissingWeights(&'static str),
    InvalidMethod(String),
    DimensionError(String),
}

impl fmt::Display for AvgVelocityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingWeights(m) => write!(f, "{}", m),
            Self::InvalidMethod(m) => write!(f, "Incorrect method given: {}", m),
            Self::DimensionError(e) => write!(f, "{}", e),
        }
    }
}

impl std::error::Error for AvgVelocityError {}

/// Averaging method enum
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AveragingMethod {
    SimpleMean,
    CubicMean,
    SimpleCubature,
    CubicCubature,
}

impl std::fmt::Display for AveragingMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AveragingMethod::SimpleMean => write!(f, "simple-mean"),
            AveragingMethod::CubicMean => write!(f, "cubic-mean"),
            AveragingMethod::SimpleCubature => write!(f, "simple-cubature"),
            AveragingMethod::CubicCubature => write!(f, "cubic-cubature"),
        }
    }
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
// ================= 核心均值函数 (支持多轴) =================

fn mean_axes_impl(array: &ArrayD, axes: &[usize]) -> ArrayD {
    if axes.is_empty() {
        return array.to_owned();
    }
    // 降序排序防止维度收缩后索引错位
    let mut sorted = axes.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));

    let mut res = array.to_owned();
    for &ax in &sorted {
        res = res.mean_axis(Axis(ax)).expect("Mean reduction failed");
    }
    res
}

pub fn simple_mean(array: &ArrayD, axes: &[usize]) -> ArrayD {
    mean_axes_impl(array, axes)
}

pub fn cubic_mean(array: &ArrayD, axes: &[usize]) -> ArrayD {
    mean_axes_impl(&array.mapv(|x| x.powi(3)), axes).mapv(|x| x.cbrt())
}

fn prepare_cubature_weights(array: &ArrayD, weights: &Array1) -> ArrayD {
    let ndim = array.ndim();
    assert!(
        ndim >= 3,
        "Array must have ≥3 dimensions for cubature integration"
    );

    let mut w = weights.to_owned();
    let w_len = w.len();
    let w_sum = w.sum();
    if w_sum == 0.0 {
        panic!("Cubature weights sum to zero");
    }
    // 归一化：使权重算术平均值为 1 (等价 Python: weights * len(weights) / sum(weights))
    let normalization_factor = w_len as f64 / w_sum;
    w.mapv_inplace(|v| v * normalization_factor);

    // 广播至 Axis(2)，对应 Python 的 weights[None, None, :, None]
    let mut bcast_shape = vec![1; ndim];
    bcast_shape[2] = w.len();
    w.into_shape(bcast_shape)
        .unwrap()
        .broadcast(array.dim())
        .unwrap()
        .to_owned()
}

pub fn simple_cubature(array: &ArrayD, cubature_weights: &Array1, axes: &[usize]) -> ArrayD {
    let w_broad = prepare_cubature_weights(array, cubature_weights);
    mean_axes_impl(&(array * w_broad), axes)
}

pub fn cubic_cubature(array: &ArrayD, cubature_weights: &Array1, axes: &[usize]) -> ArrayD {
    let w_broad = prepare_cubature_weights(array, cubature_weights);
    let prod = array.mapv(|x| x.powi(3)) * w_broad;
    mean_axes_impl(&prod, axes).mapv(|x| x.cbrt())
}

// ================= 主函数：average_velocity =================

/// 计算转子扫掠面内的平均有效风速
pub fn average_velocity(
    velocities: &ArrayD,
    ix_filter: Option<&[usize]>,
    method: &AveragingMethod,
    cubature_weights: Option<&Array1>,
) -> Result<Array2, AvgVelocityError> {
    // 已验证 2026-05-13
    // 1. 处理 ix_filter (等效 Python: velocities[:, ix_filter])
    let filtered = if let Some(indices) = ix_filter {
        // 边界检查
        let shape = velocities.shape();
        if shape.len() < 2 {
            return Err(AvgVelocityError::DimensionError(
                "Velocities must have at least 2 dimensions".into(),
            ));
        }
        let n_turbines = shape[1];
        for &idx in indices {
            if idx >= n_turbines {
                return Err(AvgVelocityError::DimensionError(format!(
                    "Turbine index {} out of bounds (max {})",
                    idx,
                    n_turbines - 1
                )));
            }
        }
        // 高级索引会创建新数组（与 Python 行为一致）
        velocities.select(Axis(1), indices)
    } else {
        velocities.to_owned()
    };

    // 2. 动态计算网格轴 (等效 Python: axis = tuple([2 + i for i in range(velocities.ndim - 2)]))
    let ndim = filtered.ndim();
    if ndim < 2 {
        return Err(AvgVelocityError::DimensionError(
            "Array must have ≥2 dimensions".into(),
        ));
    }
    let axes: Vec<usize> = (2..ndim).collect();

    // 3. 方法分发
    match method {
        AveragingMethod::SimpleMean => {
            let result = simple_mean(&filtered, &axes);
            Ok(convert_to_array2(result))
        }
        AveragingMethod::CubicMean => {
            let result = cubic_mean(&filtered, &axes);
            Ok(convert_to_array2(result))
        }
        AveragingMethod::SimpleCubature => {
            let w = cubature_weights.ok_or(AvgVelocityError::MissingWeights(
                "cubature_weights is required for 'simple-cubature' method.",
            ))?;
            let result = simple_cubature(&filtered, w, &axes);
            Ok(convert_to_array2(result))
        }
        AveragingMethod::CubicCubature => {
            let w = cubature_weights.ok_or(AvgVelocityError::MissingWeights(
                "cubature_weights is required for 'cubic-cubature' method.",
            ))?;
            let result = cubic_cubature(&filtered, w, &axes);
            Ok(convert_to_array2(result))
        }
        _ => Err(AvgVelocityError::InvalidMethod(method.to_string())),
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
    let velocities_dyn = velocities.view().into_dyn().to_owned();
    let cubature_weights_1d: Option<crate::types::Array1> =
        cubature_weights.map(|w| crate::types::Array1::from_iter(w.iter().copied()));
    let avg_vel_2d = average_velocity(&velocities_dyn, None, &method, cubature_weights_1d.as_ref())
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Convert Array2 to Array3 by adding a third dimension
    let mut avg_vel = Array3::zeros((avg_vel_2d.shape()[0], avg_vel_2d.shape()[1], 1));
    for i in 0..avg_vel_2d.shape()[0] {
        for j in 0..avg_vel_2d.shape()[1] {
            avg_vel[[i, j, 0]] = avg_vel_2d[[i, j]];
        }
    }

    // Apply yaw correction if parameters provided
    if let (Some(yaw), Some(exp_yaw)) = (yaw_angles, cosine_loss_exponent_yaw) {
        // Vectorized yaw correction
        let pw = exp_yaw / 3.0;
        let yaw_correction = yaw.mapv(|yaw_deg| cosd(yaw_deg).powf(pw));

        // Apply correction to all findex and turbine combinations
        for j in 0..avg_vel.shape()[1] {
            for i in 0..avg_vel.shape()[0] {
                avg_vel[[i, j, 0]] *= yaw_correction[[i, j]];
            }
        }
    }

    // Apply tilt correction if parameters provided
    if let (Some(tilt), Some(ref_t), Some(exp_tilt), Some(correct)) = (
        tilt_angles,
        ref_tilt,
        cosine_loss_exponent_tilt,
        correct_cp_ct_for_tilt,
    ) {
        // Precompute tilt corrections
        let exponent = exp_tilt / 3.0;
        let tilt_cos = tilt.mapv(|t| cosd(t));
        let ref_tilt_cos = ref_t.mapv(|t| cosd(t));
        let tilt_correction = (&tilt_cos / &ref_tilt_cos).mapv(|x| x.powf(exponent));

        // Apply corrections with masking
        for j in 0..avg_vel.shape()[1] {
            for i in 0..avg_vel.shape()[0] {
                if correct[[i, j]] {
                    avg_vel[[i, j, 0]] *= tilt_correction[[i, j]];
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
    // todo!("Implement tilt angle interpolation for floating turbines")
    // 目前floirs 还没有实现 2026-05-13
    _tilt_angles.clone()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array;

    #[test]
    fn test_simple_mean() {
        let data = Array::from_shape_vec(
            (2, 2, 3),
            vec![
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
            ],
        )
        .unwrap();

        let data_dyn = data.into_dyn();
        let result = simple_mean(&data_dyn, &[2]);
        let result_2d = convert_to_array2(result);
        assert_eq!(result_2d.shape(), &[2, 2]);
        assert_relative_eq!(result_2d[[0, 0]], 2.0);
        assert_relative_eq!(result_2d[[0, 1]], 5.0);
    }

    #[test]
    fn test_cubic_mean() {
        let data = Array::from_shape_vec((1, 1, 3), vec![1.0, 2.0, 3.0]).unwrap();

        let data_dyn = data.into_dyn();
        let result = cubic_mean(&data_dyn, &[2]);
        let result_2d = convert_to_array2(result);
        let expected = ((1.0_f64.powi(3) + 2.0_f64.powi(3) + 3.0_f64.powi(3)) / 3.0).cbrt();
        assert_relative_eq!(result_2d[[0, 0]], expected, epsilon = 1e-10);
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
        assert_eq!(AveragingMethod::from_str("invalid"), None);
    }

    #[test]
    fn test_rotor_velocity_tilt_cosine_correction() {
        // Test case matching Python Floris behavior
        let n_turbines = 4;
        let wind_speed = Array2::from_elem((1, n_turbines), 10.0);

        let tilt_angles = Array2::from_elem((1, n_turbines), 5.0);
        let ref_tilt = Array2::from_elem((1, n_turbines), 5.0);
        let cosine_loss_exponent_tilt = 1.88;

        // Test 1: When correct_cp_ct_for_tilt is all False, velocities should remain unchanged
        // (because tilt_angles == ref_tilt, so correction factor = 1.0)
        let correct_cp_ct_for_tilt = ndarray::Array2::<bool>::from_elem((1, n_turbines), false);
        let result = rotor_velocity_tilt_cosine_correction(
            &tilt_angles,
            &ref_tilt,
            cosine_loss_exponent_tilt,
            &correct_cp_ct_for_tilt,
            &wind_speed,
        );

        // Since tilt_angles == ref_tilt, cosd(5)/cosd(5) = 1.0, so result should equal input
        for i in 0..result.nrows() {
            for j in 0..result.ncols() {
                assert_relative_eq!(result[[i, j]], wind_speed[[i, j]], epsilon = 1e-10);
            }
        }

        // Test 2: When tilt angle > ref_tilt, velocity should decrease
        let tilt_angles_back = Array2::from_elem((1, n_turbines), 10.0); // Greater than ref_tilt (5.0)
        let correct_all_true = ndarray::Array2::<bool>::from_elem((1, n_turbines), true);
        let result_back = rotor_velocity_tilt_cosine_correction(
            &tilt_angles_back,
            &ref_tilt,
            cosine_loss_exponent_tilt,
            &correct_all_true,
            &wind_speed,
        );

        for i in 0..result_back.nrows() {
            for j in 0..result_back.ncols() {
                assert!(result_back[[i, j]] < wind_speed[[i, j]]);
            }
        }

        // Test 3: When tilt angle < ref_tilt, velocity should increase
        let tilt_angles_forward = Array2::from_elem((1, n_turbines), 0.0); // Less than ref_tilt (5.0)
        let result_forward = rotor_velocity_tilt_cosine_correction(
            &tilt_angles_forward,
            &ref_tilt,
            cosine_loss_exponent_tilt,
            &correct_all_true,
            &wind_speed,
        );

        for i in 0..result_forward.nrows() {
            for j in 0..result_forward.ncols() {
                assert!(result_forward[[i, j]] > wind_speed[[i, j]]);
            }
        }

        // Test 4: Symmetry around zero tilt
        let tilt_angle_sym = 3.0;
        let tilt_negative = Array2::from_elem((1, n_turbines), -tilt_angle_sym);
        let tilt_positive = Array2::from_elem((1, n_turbines), tilt_angle_sym);

        let result_negative = rotor_velocity_tilt_cosine_correction(
            &tilt_negative,
            &ref_tilt,
            cosine_loss_exponent_tilt,
            &correct_all_true,
            &wind_speed,
        );

        let result_positive = rotor_velocity_tilt_cosine_correction(
            &tilt_positive,
            &ref_tilt,
            cosine_loss_exponent_tilt,
            &correct_all_true,
            &wind_speed,
        );

        // cos(-x) == cos(x), so results should be equal
        for i in 0..result_negative.nrows() {
            for j in 0..result_negative.ncols() {
                assert_relative_eq!(
                    result_negative[[i, j]],
                    result_positive[[i, j]],
                    epsilon = 1e-10
                );
            }
        }
    }
}
