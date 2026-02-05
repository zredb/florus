//! Wind data structures for FLORIS-RS
//!
//! Provides wind data objects to hold ambient wind conditions including:
//! - TimeSeries: Time series wind data
//! - WindRose: Aggregated wind statistics by direction/speed bins
//! - WindTIRose: Wind rose with TI as an additional dimension
//! - WindRoseWRG: Wind rose from WAsP WRG file
//! - WindRoseByTurbine: Wind rose with separate wind rose for each turbine
//!
//! Corresponds to wind_data.py in Python FLORIS v4.6

pub mod timeseries;
pub mod traits;
pub mod wind_rose;
pub mod wind_rose_wrg;
pub mod wind_ti_rose;

use crate::types::ArrayView1;
use ordered_float::OrderedFloat;
pub use timeseries::TimeSeries;
pub use traits::{TIParams, WindData};
pub use wind_rose::WindRose;
pub use wind_rose_wrg::{RegularGridInterpolant, WRGData, WindRoseByTurbine, WindRoseWRG};
pub use wind_ti_rose::WindTIRose;

use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum ValidationError {
    #[error("Array must contain at least 2 elements")]
    InsufficientElements,

    #[error("{0} must be monotonically increasing")]
    NotMonotonicallyIncreasing(String),

    #[error("{0} must be evenly spaced")]
    NotEvenlySpaced(String),

    #[error("Wind directions have invalid step pattern")]
    InvalidWindDirectionSteps,

    #[error("Invalid shape for array with dimensions ({0}, {1})")]
    InvalidShape2(usize, usize),
}

// 结果类型别名
type ValidationResult<T> = Result<T, ValidationError>;

pub fn validate_wind_directions(wind_directions: &ArrayView1) -> ValidationResult<f64> {
    let n = wind_directions.len();

    if n < 2 {
        return Err(ValidationError::InsufficientElements);
    }

    if n == 1 {
        return Ok(0.0); // 单个元素，步长为0
    }

    // 检查单调递增
    for i in 1..n {
        if wind_directions[i] <= wind_directions[i - 1] {
            return Err(ValidationError::NotMonotonicallyIncreasing(
                "wind_directions".to_string(),
            ));
        }
    }

    // 检查等间距（考虑循环）
    check_and_identify_wind_direction_step(wind_directions)
}

/// 检查风速数组的有效性
pub fn validate_wind_speeds(wind_speeds: &ArrayView1) -> ValidationResult<f64> {
    let n = wind_speeds.len();

    if n < 2 {
        return Err(ValidationError::InsufficientElements);
    }

    if n == 1 {
        return Ok(0.0); // 单个元素，步长为0
    }

    // 检查单调递增
    for i in 1..n {
        if wind_speeds[i] <= wind_speeds[i - 1] {
            return Err(ValidationError::NotMonotonicallyIncreasing(
                "wind_speeds".to_string(),
            ));
        }
    }

    // 检查等间距
    check_and_identify_step_size(wind_speeds, false)
}

/// 主验证函数（对应Python原逻辑）
pub fn validate_wind_arrays(
    wind_directions: &ArrayView1,
    wind_speeds: &ArrayView1,
) -> ValidationResult<(f64, f64)> {
    let n_wd = wind_directions.len();
    let n_ws = wind_speeds.len();

    let mut wd_step = 0.0;
    let mut ws_step = 0.0;

    if n_wd > 1 {
        // 检查风向单调递增
        if !is_monotonically_increasing(wind_directions) {
            return Err(ValidationError::NotMonotonicallyIncreasing(
                "wind_directions".to_string(),
            ));
        }

        // 检查等间距（考虑循环）
        wd_step = check_and_identify_wind_direction_step(wind_directions)?;
    }

    if n_ws > 1 {
        // 检查风速单调递增
        if !is_monotonically_increasing(wind_speeds) {
            return Err(ValidationError::NotMonotonicallyIncreasing(
                "wind_speeds".to_string(),
            ));
        }

        // 检查等间距
        ws_step = check_and_identify_step_size(wind_speeds, true)?;
    }

    Ok((wd_step, ws_step))
}

// ========== 核心辅助函数 ==========

/// 检查数组是否严格单调递增
fn is_monotonically_increasing(arr: &ArrayView1) -> bool {
    arr.windows(2).into_iter().all(|w| w[1] > w[0])
}

/// 检查并识别风向步长（特殊处理360°循环）
fn check_and_identify_wind_direction_step(wind_directions: &ArrayView1) -> ValidationResult<f64> {
    let n = wind_directions.len();

    if n < 2 {
        return Err(ValidationError::InsufficientElements);
    }

    // 计算内部步长
    let steps: Vec<f64> = wind_directions
        .windows(2)
        .into_iter()
        .map(|w| w[1] - w[0])
        .collect();

    // 确认所有内部步长为正（已由单调性检查保证）
    if !steps.iter().all(|&s| s > 0.0) {
        return Err(ValidationError::NotMonotonicallyIncreasing(
            "wind_directions".to_string(),
        ));
    }

    // 计算最后一个元素到第一个元素的循环步长（加上360°）
    let last_step = wind_directions[0] - wind_directions[n - 1] + 360.0;

    match n {
        2 => {
            // 对于两个元素，返回较小的步长
            Ok(steps[0].min(last_step))
        }
        3 => {
            // 对于三个元素的特殊处理
            if (steps[0] - steps[1]).abs() < f64::EPSILON {
                // 如果两个内部步长相等
                Ok(steps[0])
            } else if (steps[0] - last_step).abs() < f64::EPSILON {
                Ok(steps[0])
            } else if (steps[1] - last_step).abs() < f64::EPSILON {
                Ok(steps[1])
            } else {
                Err(ValidationError::InvalidWindDirectionSteps)
            }
        }
        _ => {
            // 对于更多元素
            if are_all_close(&steps, steps[0]) {
                // 所有内部步长相等
                return Ok(steps[0]);
            }

            // 统计步长出现频率
            let mut step_counts = std::collections::HashMap::new();
            for &step in &steps {
                let key = OrderedFloat::from(round_step(step));
                *step_counts.entry(key).or_insert(0) += 1;
            }

            // 如果只有一种步长（已由上面处理）
            // 检查是否有两种步长，且其中一种只出现一次
            if step_counts.len() == 2 {
                let counts: Vec<&i32> = step_counts.values().collect();
                if *counts[0] == 1 || *counts[1] == 1 {
                    // 找到最常见的步长
                    let most_common_step = step_counts
                        .iter()
                        .max_by_key(|(_, &count)| count)
                        .map(|(&step, _)| step)
                        .unwrap();

                    // 检查循环步长是否等于最常见的步长
                    let rounded_last_step = OrderedFloat::from(round_step(last_step));
                    if (rounded_last_step - most_common_step).abs() < f64::EPSILON {
                        return Ok(most_common_step.into_inner());
                    }
                }
            }

            Err(ValidationError::NotEvenlySpaced(
                "wind_directions".to_string(),
            ))
        }
    }
}

/// 检查并识别普通数组的步长
fn check_and_identify_step_size(arr: &ArrayView1, is_wind_speed: bool) -> ValidationResult<f64> {
    let n = arr.len();

    if n < 2 {
        return Err(ValidationError::InsufficientElements);
    }

    // 计算步长
    let steps: Vec<f64> = arr.windows(2).into_iter().map(|w| w[1] - w[0]).collect();

    // 检查是否所有步长都近似相等
    if are_all_close(&steps, steps[0]) {
        Ok(steps[0])
    } else {
        let name = if is_wind_speed {
            "wind_speeds"
        } else {
            "generic array"
        };
        Err(ValidationError::NotEvenlySpaced(name.to_string()))
    }
}

/// 检查数组中所有值是否接近参考值
fn are_all_close(values: &[f64], reference: f64) -> bool {
    const EPSILON: f64 = 1e-10;
    values.iter().all(|&v| (v - reference).abs() < EPSILON)
}

/// 四舍五入步长到合理精度，用于比较
fn round_step(step: f64) -> f64 {
    (step * 1e9).round() / 1e9
}

// ========== 单元测试 ==========

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_wind_speeds_valid() {
        let speeds = vec![3.0, 4.0, 5.0, 6.0];
        let speeds = ArrayView1::from(&speeds);
        let result = validate_wind_speeds(&speeds);
        assert!(result.is_ok());
        assert!((result.unwrap() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_validate_wind_speeds_not_monotonic() {
        let speeds = vec![3.0, 5.0, 4.0, 6.0];
        let speeds = ArrayView1::from(&speeds);
        let result = validate_wind_speeds(&speeds);
        assert_eq!(
            result,
            Err(ValidationError::NotMonotonicallyIncreasing(
                "wind_speeds".to_string()
            ))
        );
    }

    #[test]
    fn test_validate_wind_speeds_not_evenly_spaced() {
        let speeds = vec![3.0, 4.5, 5.0, 6.0];
        let speeds = ArrayView1::from(&speeds);
        let result = validate_wind_speeds(&speeds);
        assert_eq!(
            result,
            Err(ValidationError::NotEvenlySpaced("wind_speeds".to_string()))
        );
    }

    #[test]
    fn test_validate_wind_directions_cyclic() {
        // 跨越360°的扇区
        let directions = vec![350.0, 0.0, 10.0];
        let directions = ArrayView1::from(&directions);
        let result = validate_wind_directions(&directions);
        assert!(result.is_ok());
        assert!((result.unwrap() - 10.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_validate_wind_directions_simple() {
        let directions = vec![0.0, 30.0, 60.0, 90.0];
        let directions = ArrayView1::from(&directions);
        let result = validate_wind_directions(&directions);
        assert!(result.is_ok());
        assert!((result.unwrap() - 30.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_validate_wind_directions_invalid_cyclic() {
        // 无效的循环模式
        let directions = vec![340.0, 350.0, 0.0, 20.0];
        let directions = ArrayView1::from(&directions);
        let result = validate_wind_directions(&directions);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_wind_arrays_mixed() {
        let directions = vec![0.0, 30.0, 60.0];
        let speeds = vec![4.0, 5.0, 6.0, 7.0];

        let directions = ArrayView1::from(&directions);
        let speeds = ArrayView1::from(&speeds);

        let result = validate_wind_arrays(&directions, &speeds);
        assert!(result.is_ok());

        let (wd_step, ws_step) = result.unwrap();
        assert!((wd_step - 30.0).abs() < f64::EPSILON);
        assert!((ws_step - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_single_element_arrays() {
        // 单个元素应该通过验证
        let directions = vec![180.0];
        let speeds = vec![8.0];
        let directions = ArrayView1::from(&directions);
        let speeds = ArrayView1::from(&speeds);
        let result = validate_wind_arrays(&directions, &speeds);
        assert!(result.is_ok());

        let (wd_step, ws_step) = result.unwrap();
        assert!((wd_step - 0.0).abs() < f64::EPSILON);
        assert!((ws_step - 0.0).abs() < f64::EPSILON);
    }
}
