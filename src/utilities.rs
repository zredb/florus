//! Utility functions for FLORIS-RS
//!
//! Corresponds to utilities.py in the Python implementation

use crate::types::{Float, Array1, Array2, Array4};
use anyhow::{Context, Result};
use serde_yaml::Value;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use ndarray::{Array, ArrayView1};

/// Small value for floating-point comparisons
const EPSILON: f64 = 1e-6;

/// Load YAML configuration from file
pub fn load_yaml<P: AsRef<Path>>(path: P) -> Result<Value> {
    let mut file = File::open(path.as_ref())
        .with_context(|| format!("Failed to open file: {:?}", path.as_ref()))?;
    
    let mut contents = String::new();
    file.read_to_string(&mut contents)
        .context("Failed to read file contents")?;
    
    let value: Value = serde_yaml::from_str(&contents)
        .context("Failed to parse YAML")?;
    
    Ok(value)
}

/// Get nested value from configuration dictionary
pub fn nested_get<'a>(dict: &'a Value, keys: &[&str]) -> Option<&'a Value> {
    let mut current = dict;
    
    for key in keys {
        current = current.get(key)?;
    }
    
    Some(current)
}

/// Set nested value in configuration dictionary
pub fn nested_set(dict: &mut Value, keys: &[&str], value: Value) -> Result<()> {
    if keys.is_empty() {
        anyhow::bail!("Keys cannot be empty");
    }
    
    let mut current = dict;
    
    // Navigate to the parent of the final key
    for key in &keys[..keys.len() - 1] {
        current = current
            .get_mut(key)
            .ok_or_else(|| anyhow::anyhow!("Key not found: {}", key))?;
    }
    
    // Set the value at the final key
    let final_key = keys[keys.len() - 1];
    if let Some(map) = current.as_mapping_mut() {
        map.insert(Value::String(final_key.to_string()), value);
        Ok(())
    } else {
        anyhow::bail!("Cannot set value: parent is not a mapping")
    }
}

/// Rotate coordinates relative to west (270 degrees)
/// 
/// # Arguments
/// * `wind_directions` - Wind directions in degrees (0 = North)
/// * `turbine_coordinates` - Turbine coordinates (n_turbines, 3)
/// 
/// # Returns
/// Tuple of (x, y, z, x_center, y_center) where coordinates are rotated
pub fn rotate_coordinates_rel_west(
    wind_directions: &Array1,
    turbine_coordinates: &Array2,
) -> Result<(Array2, Array2, Array2, Array2, Array2)> {
    let n_findex = wind_directions.len();
    let n_turbines = turbine_coordinates.shape()[0];
    
    // Calculate center of rotation (centroid of turbines)
    let x_center = turbine_coordinates.column(0).mean().unwrap();
    let y_center = turbine_coordinates.column(1).mean().unwrap();
    
    let mut x = Array::zeros((n_findex, n_turbines));
    let mut y = Array::zeros((n_findex, n_turbines));
    let mut z = Array::zeros((n_findex, n_turbines));
    let mut x_center_of_rotation = Array::zeros((n_findex, n_turbines));
    let mut y_center_of_rotation = Array::zeros((n_findex, n_turbines));
    
    for fi in 0..n_findex {
        let wd = wind_directions[fi];
        // Rotate coordinates relative to west (270 deg)
        let angle = (wd - 270.0).to_radians();
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        
        for ti in 0..n_turbines {
            let x_orig = turbine_coordinates[[ti, 0]] - x_center;
            let y_orig = turbine_coordinates[[ti, 1]] - y_center;
            
            x[[fi, ti]] = x_orig * cos_angle - y_orig * sin_angle + x_center;
            y[[fi, ti]] = x_orig * sin_angle + y_orig * cos_angle + y_center;
            z[[fi, ti]] = turbine_coordinates[[ti, 2]];
            
            x_center_of_rotation[[fi, ti]] = x_center;
            y_center_of_rotation[[fi, ti]] = y_center;
        }
    }
    
    Ok((x, y, z, x_center_of_rotation, y_center_of_rotation))
}

/// Reverse rotate coordinates from wind frame back to inertial frame
/// 
/// # Arguments
/// * `wind_directions` - Wind directions in degrees
/// * `grid_x`, `grid_y`, `grid_z` - Grid coordinates in wind frame
/// * `x_center_of_rotation`, `y_center_of_rotation` - Centers of rotation
pub fn reverse_rotate_coordinates_rel_west(
    wind_directions: &Array1,
    grid_x: &Array4,
    grid_y: &Array4,
    grid_z: &Array4,
    x_center_of_rotation: &Array2,
    y_center_of_rotation: &Array2,
) -> Result<(Array4, Array4, Array4)> {
    let shape = grid_x.shape();
    let n_findex = shape[0];
    
    let mut x_inertial = grid_x.clone();
    let mut y_inertial = grid_y.clone();
    let z_inertial = grid_z.clone();
    
    for fi in 0..n_findex {
        let wd = wind_directions[fi];
        // Reverse rotation
        let angle = -(wd - 270.0).to_radians();
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        
        for ti in 0..shape[1] {
            let x_center = x_center_of_rotation[[fi, ti]];
            let y_center = y_center_of_rotation[[fi, ti]];
            
            for i in 0..shape[2] {
                for j in 0..shape[3] {
                    let x_rel = grid_x[[fi, ti, i, j]] - x_center;
                    let y_rel = grid_y[[fi, ti, i, j]] - y_center;
                    
                    x_inertial[[fi, ti, i, j]] = x_rel * cos_angle - y_rel * sin_angle + x_center;
                    y_inertial[[fi, ti, i, j]] = x_rel * sin_angle + y_rel * cos_angle + y_center;
                }
            }
        }
    }
    
    Ok((x_inertial, y_inertial, z_inertial))
}

/// Rotate coordinates relative to west (270 degrees) - legacy version for Vec inputs
/// 
/// # Arguments
/// * `x` - X coordinates
/// * `y` - Y coordinates  
/// * `wind_direction` - Wind direction in degrees (0 = North)
pub fn reverse_rotate_coordinates_rel_west_vec(
    x: &[Float],
    y: &[Float],
    wind_direction: Float,
) -> Result<(Vec<Float>, Vec<Float>)> {
    if x.len() != y.len() {
        anyhow::bail!("x and y must have same length");
    }
    
    // Convert wind direction to radians and adjust for rotation relative to west
    let angle = (wind_direction - 270.0).to_radians();
    let cos_angle = angle.cos();
    let sin_angle = angle.sin();
    
    let mut x_rot = Vec::with_capacity(x.len());
    let mut y_rot = Vec::with_capacity(y.len());
    
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        x_rot.push(xi * cos_angle - yi * sin_angle);
        y_rot.push(xi * sin_angle + yi * cos_angle);
    }
    
    Ok((x_rot, y_rot))
}

/// Cosd: cosine of angle in degrees
#[inline]
pub fn cosd(degrees: Float) -> Float {
    degrees.to_radians().cos()
}

/// Sind: sine of angle in degrees
#[inline]
pub fn sind(degrees: Float) -> Float {
    degrees.to_radians().sin()
}

/// Tand: tangent of angle in degrees
#[inline]
pub fn tand(degrees: Float) -> Float {
    degrees.to_radians().tan()
}

/// Wrap angle to [0, 360) degrees
pub fn wrap_360(angle: Float) -> Float {
    let mut result = angle % 360.0;
    if result < 0.0 {
        result += 360.0;
    }
    result
}

/// Wrap angle to [-180, 180) degrees
pub fn wrap_180(angle: Float) -> Float {
    let mut result = (angle + 180.0) % 360.0;
    if result < 0.0 {
        result += 360.0;
    }
    result - 180.0
}


/// Identifies the step size in a series of wind directions.
/// Returns the step size if the wind directions are evenly spaced, otherwise returns an error.
/// Handles circular wind direction data (e.g., [330, 0, 30] wraps around 360)
///
/// # Arguments
///
/// * `wind_directions` - Array of wind directions (ndarray)
///
/// # Returns
///
/// The step size of the wind directions
///
/// # Errors
///
/// Returns a `String` error message if:
/// * Array contains less than 2 elements
/// * Wind directions are not evenly spaced (considering circular nature)
pub fn check_and_identify_step_size(wind_directions: ArrayView1<f64>) -> Result<f64, String> {
    if wind_directions.len() < 2 {
        return Err("Array must contain at least 2 elements".to_string());
    }

    let wind_dirs_slice = wind_directions.as_slice().unwrap();

    // Check for monotonicity first
    let negative_steps: Vec<usize> = wind_dirs_slice.windows(2)
        .enumerate()
        .filter(|(_, window)| window[1] - window[0] < -EPSILON)
        .map(|(i, _)| i)
        .collect();

    // If there are negative steps, this is circular data
    if !negative_steps.is_empty() {
        // For circular data, we need to handle the wrap-around
        // Compute steps treating wrap-around as positive
        let mut circular_steps: Vec<f64> = Vec::new();
        for i in 0..wind_dirs_slice.len() - 1 {
            let mut step = wind_dirs_slice[i + 1] - wind_dirs_slice[i];
            if step < -EPSILON {
                step += 360.0; // Handle wrap-around
            }
            circular_steps.push(step);
        }

        // For circular data, we only check that all internal steps are approximately equal
        // The wrap step (from last to first going backwards) doesn't need to match for circular data
        if circular_steps.is_empty() {
            return Err("Invalid circular wind directions".to_string());
        }

        let first_step = circular_steps[0];

        // Check if all internal steps are approximately equal
        if circular_steps.iter().all(|&step| (step - first_step).abs() < EPSILON) {
            return Ok(first_step);
        } else {
            return Err("wind_directions must be evenly spaced".to_string());
        }
    }

    // Non-circular case: all steps must be positive
    let steps: Vec<f64> = wind_dirs_slice.windows(2)
        .map(|window| window[1] - window[0])
        .collect();

    // Confirm that the steps are all positive
    if !steps.iter().all(|&step| step > 0.0) {
        return Err("wind_directions must be monotonically increasing".to_string());
    }

    // Check the step from the last to the first element
    let last_step = wind_directions[0] - wind_directions[wind_directions.len() - 1] + 360.0;

    // If len(wind_directions) == 2, then return whichever step is smaller
    if wind_directions.len() == 2 {
        return Ok(steps[0].min(last_step));
    }

    // If len(wind_directions) == 3 make some checks
    if wind_directions.len() == 3 {
        let first_step = steps[0];
        let second_step = steps[1];

        if (first_step - second_step).abs() < EPSILON {
            return Ok(first_step);
        } else if (first_step - last_step).abs() < EPSILON {
            return Ok(first_step);
        } else if (second_step - last_step).abs() < EPSILON {
            return Ok(second_step);
        } else {
            return Err("wind_directions must be evenly spaced".to_string());
        }
    }

    // For len > 3
    // Check if all steps are approximately equal
    let first_step = steps[0];
    if steps.iter().all(|&step| (step - first_step).abs() < EPSILON) {
        return Ok(first_step);
    }

    // Count the frequency of each step value
    let mut step_values: Vec<(f64, usize)> = Vec::new();
    for &step in &steps {
        let mut found = false;
        for (value, count) in &mut step_values {
            if (step - *value).abs() < EPSILON * 100.0 {
                *count += 1;
                found = true;
                break;
            }
        }
        if !found {
            step_values.push((step, 1));
        }
    }

    // Check for the case where there are more than two different step sizes
    if step_values.len() > 2 {
        return Err("wind_directions must be evenly spaced".to_string());
    }

    // Find the most common step size
    let most_common_step = step_values.iter()
        .max_by(|a, b| a.1.cmp(&b.1))
        .map(|(value, _)| *value)
        .unwrap();

    // In the case there are only two step sizes, ensure that one only happens once
    if step_values.len() == 2 {
        let min_count = step_values.iter().map(|(_, count)| *count).min().unwrap();
        if min_count > 1 {
            return Err("wind_directions must be evenly spaced".to_string());
        }
    }

    // If the last step equals the most common step, return the most common step
    if (last_step - most_common_step).abs() < EPSILON {
        return Ok(most_common_step);
    }

    Err("wind_directions must be evenly spaced".to_string())
}

/// Reorders the wind directions so that they are adjacent.
/// Returns the reordered wind directions if the wind directions are not adjacent,
/// otherwise returns the input wind directions.
/// Handles circular wind direction data.
///
/// # Arguments
///
/// * `wind_directions` - Array of wind directions (ndarray)
///
/// # Returns
///
/// A tuple containing:
/// * The reordered wind directions to be adjacent
/// * Sort indices to go from the original to the new array
pub fn make_wind_directions_adjacent(wind_directions: ArrayView1<f64>) -> (Array1, Vec<usize>) {
    // Check the step size of the wind directions
    let step_size = check_and_identify_step_size(wind_directions)
        .expect("wind_directions must have a valid step size");

    let wind_dirs_slice = wind_directions.as_slice().unwrap();
    let n = wind_directions.len();

    // Check if data has negative steps (non-monotonic, indicating circular data)
    let has_negative_steps = wind_dirs_slice.windows(2)
        .any(|window| window[1] - window[0] < -EPSILON);

    // For circular data detection, we need both:
    // 1. Negative steps (non-monotonic progression)
    // 2. Large wrap step (data wraps around 360)
    let wrap_step = (wind_dirs_slice[0] + 360.0) - wind_dirs_slice[n - 1];
    let is_circular = has_negative_steps && wrap_step > 180.0;

    if is_circular {
        // For circular data, subtract 360 from all elements to make them adjacent
        // Find the wrap point (where direction wraps from high to low)
        let wrap_point = wind_dirs_slice.windows(2)
            .position(|window| window[1] - window[0] < -EPSILON)
            .map(|i| i + 1)
            .unwrap_or(0);

        // Create reordered array: elements after wrap point minus 360, then elements up to wrap point
        let mut reordered = Array1::zeros(n);

        // Fill first part: elements after wrap point minus 360
        for (i, j) in (wrap_point..n).enumerate() {
            reordered[i] = wind_directions[j] - 360.0;
        }

        // Fill second part: elements up to wrap point
        for (i, j) in (0..wrap_point).enumerate() {
            reordered[n - wrap_point + i] = wind_directions[j];
        }

        // Sort indices: elements after wrap point come first, then elements up to wrap point
        let sort_indices: Vec<usize> = (wrap_point..n).chain(0..wrap_point).collect();

        (reordered, sort_indices)
    } else {
        // Non-circular case: original logic
        let steps: Vec<f64> = wind_dirs_slice.windows(2)
            .map(|window| window[1] - window[0])
            .collect();

        // There will be at most one step with a size larger than the step size
        // If there is one, find it
        let large_step_idx = steps.iter()
            .position(|&step| step > step_size + EPSILON);

        match large_step_idx {
            Some(idx) => {
                // Now change wind_directions such that for each direction after that index
                // subtract 360 and move that block to the front
                let mut reordered = Array1::zeros(n);

                // Fill the first part: wind_directions[idx+1..] - 360
                for (i, j) in ((idx + 1)..n).enumerate() {
                    reordered[i] = wind_directions[j] - 360.0;
                }

                // Fill the second part: wind_directions[..idx+1]
                for (i, j) in (0..=idx).enumerate() {
                    reordered[n - idx - 1 + i] = wind_directions[j];
                }

                // Return the wind directions and indices to go from the original to the new
                let sort_indices: Vec<usize> = (idx + 1..n).chain(0..=idx).collect();

                (reordered, sort_indices)
            }
            None => {
                // Return the original wind directions with sequential indices
                let sort_indices: Vec<usize> = (0..n).collect();

                (wind_directions.to_owned(), sort_indices)
            }
        }
    }
}

/// Calculates the deviation from West (270 degrees) for a single wind direction.
/// The result is always a positive value between 0 and 360.
///
/// # Arguments
///
/// * `wind_direction` - A single wind direction (can be any number, negative or positive)
///
/// # Returns
///
/// The delta between the given wind direction and 270, as a positive value between 0 and 360.
pub fn wind_delta_scalar(wind_direction: f64) -> f64 {
    let result = (wind_direction - 270.0) % 360.0;
    if result < 0.0 {
        result + 360.0
    } else {
        result
    }
}

/// Calculates the deviation from West (270 degrees) for an array of wind directions.
/// The result is always positive values between 0 and 360.
///
/// # Arguments
///
/// * `wind_directions` - Array of wind directions (ndarray)
///
/// # Returns
///
/// An array of deltas between each wind direction and 270, as positive values between 0 and 360.
pub fn wind_delta_array(wind_directions: ArrayView1<f64>) -> Array1 {
    wind_directions.mapv(|dir| {
        let result = (dir - 270.0) % 360.0;
        if result < 0.0 {
            result + 360.0
        } else {
            result
        }
    })
}


#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_cosd() {
        assert_relative_eq!(cosd(0.0), 1.0);
        assert_relative_eq!(cosd(90.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(cosd(180.0), -1.0, epsilon = 1e-10);
    }
    
    #[test]
    fn test_wrap_360() {
        assert_relative_eq!(wrap_360(0.0), 0.0);
        assert_relative_eq!(wrap_360(360.0), 0.0);
        assert_relative_eq!(wrap_360(-90.0), 270.0);
        assert_relative_eq!(wrap_360(450.0), 90.0);
    }
    
    #[test]
    fn test_wrap_180() {
        assert_relative_eq!(wrap_180(0.0), 0.0);
        assert_relative_eq!(wrap_180(180.0), -180.0);
        assert_relative_eq!(wrap_180(-270.0), 90.0);
    }

     #[test]
    fn test_valid_two_elements() {
        let directions = Array1::from_vec(vec![0.0, 90.0]);
        let result = check_and_identify_step_size(directions.view());
        assert_eq!(result.unwrap(), 90.0);
    }

    #[test]
    fn test_valid_three_elements() {
        let directions = Array1::from_vec(vec![0.0, 90.0, 180.0]);
        let result = check_and_identify_step_size(directions.view());
        assert_eq!(result.unwrap(), 90.0);
    }

    #[test]
    fn test_valid_three_elements_with_last_step() {
        let directions = Array1::from_vec(vec![0.0, 90.0, 270.0]);
        let result = check_and_identify_step_size(directions.view());
        assert_eq!(result.unwrap(), 90.0);
    }

    #[test]
    fn test_valid_many_elements() {
        let directions = Array1::from_vec(vec![0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]);
        let result = check_and_identify_step_size(directions.view());
        assert_eq!(result.unwrap(), 30.0);
    }

    #[test]
    fn test_error_less_than_two() {
        let directions = Array1::from_vec(vec![0.0]);
        let result = check_and_identify_step_size(directions.view());
        assert!(result.is_err());
    }

    #[test]
    fn test_error_not_monotonic() {
        let directions = Array1::from_vec(vec![0.0, 90.0, 45.0]);
        let result = check_and_identify_step_size(directions.view());
        assert!(result.is_err());
    }

    #[test]
    fn test_error_not_evenly_spaced() {
        let directions = Array1::from_vec(vec![0.0, 30.0, 90.0, 120.0]);
        let result = check_and_identify_step_size(directions.view());
        assert!(result.is_err());
    }

    #[test]
    fn test_make_adjacent_already_adjacent() {
        let directions = Array1::from_vec(vec![0.0, 30.0, 60.0, 90.0]);
        let (reordered, indices) = make_wind_directions_adjacent(directions.view());

        assert_eq!(reordered, Array1::from_vec(vec![0.0, 30.0, 60.0, 90.0]));
        assert_eq!(indices, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_make_adjacent_with_wrap() {
        // Wind directions that wrap around 360
        let directions = Array1::from_vec(vec![330.0, 0.0, 30.0, 60.0, 90.0]);
        let (reordered, indices) = make_wind_directions_adjacent(directions.view());

        // Expected: [0-360, 30-360, 60-360, 90-360, 330] = [-360, -330, -300, -270, 330]
        // But actually: subtract 360 from elements after the large step
        // Large step is between 330 and 0 (330), so we take [0, 30, 60, 90] - 360 = [-360, -330, -300, -270]
        // And prepend [330]
        assert_eq!(reordered, Array1::from_vec(vec![-360.0, -330.0, -300.0, -270.0, 330.0]));
        assert_eq!(indices, vec![1, 2, 3, 4, 0]);
    }

    #[test]
    fn test_make_adjacent_wrap_in_middle() {
        // Wind directions with wrap in the middle - using evenly spaced circular data
        // [300, 330, 0, 30] has step 30 going: 300→330=30, 330→0(360)=30, 0→30=30
        let directions = Array1::from_vec(vec![300.0, 330.0, 0.0, 30.0]);
        let (reordered, indices) = make_wind_directions_adjacent(directions.view());

        // Wrap point is between 330 and 0, so we take [0, 30] - 360 = [-360, -330]
        // And prepend [300, 330]
        assert_eq!(reordered, Array1::from_vec(vec![-360.0, -330.0, 300.0, 330.0]));
        assert_eq!(indices, vec![2, 3, 0, 1]);
    }
}
