//! Utility functions for FLORUS
//!
//! Corresponds to utilities.py in the Python implementation

use crate::types::{Array1, Array2, Array4, Float};
use anyhow::{Context, Result};
use ndarray::{Array, ArrayView1};
use serde_yaml::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// Small value for floating-point comparisons
const EPSILON: f64 = 1e-6;

// ============================================================================
// Debug/Print Utilities
// ============================================================================

/// Print shape of an array with optional label (for debugging)
pub fn pshape<D>(array: &ndarray::ArrayBase<ndarray::OwnedRepr<f64>, D>, label: &str)
where
    D: ndarray::Dimension,
{
    println!("{} {:?}", label, array.shape());
}

// ============================================================================
// Trigonometric Functions (Degrees)
// ============================================================================

/// Cosine of an angle with the angle given in degrees.
#[inline]
pub fn cosd(degrees: Float) -> Float {
    degrees.to_radians().cos()
}

/// Sine of an angle with the angle given in degrees.
#[inline]
pub fn sind(degrees: Float) -> Float {
    degrees.to_radians().sin()
}

/// Tangent of an angle with the angle given in degrees.
#[inline]
pub fn tand(degrees: Float) -> Float {
    degrees.to_radians().tan()
}

// ============================================================================
// Angle Wrapping Functions
// ============================================================================

/// Shift the given values to within the range (-180, 180].
pub fn wrap_180(angle: Float) -> Float {
    let mut result = (angle + 180.0) % 360.0;
    if result < 0.0 {
        result += 360.0;
    }
    result - 180.0
}

/// Shift the given values to within the range (0, 360].
pub fn wrap_360(angle: Float) -> Float {
    let mut result = angle % 360.0;
    if result < 0.0 {
        result += 360.0;
    }
    result
}

// ============================================================================
// Wind Direction Analysis Functions
// ============================================================================

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
/// * Wind directions are not monotonically increasing
/// * Wind directions are not evenly spaced (considering circular nature)
pub fn check_and_identify_step_size(wind_directions: ArrayView1<f64>) -> Result<f64, String> {
    if wind_directions.len() < 2 {
        return Err("Array must contain at least 2 elements".to_string());
    }

    let wind_dirs_slice = wind_directions.as_slice().unwrap();

    // First compute the steps between each wind direction
    let steps: Vec<f64> = wind_dirs_slice
        .windows(2)
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
    if steps
        .iter()
        .all(|&step| (step - first_step).abs() < EPSILON)
    {
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

    // In the case there are only two step sizes, ensure that one only happens once
    if step_values.len() == 2 {
        let min_count = step_values.iter().map(|(_, count)| *count).min().unwrap();
        if min_count > 1 {
            return Err("wind_directions must be evenly spaced".to_string());
        }
    }

    // Find the most common step size
    let most_common_step = step_values
        .iter()
        .max_by(|a, b| a.1.cmp(&b.1))
        .map(|(value, _)| *value)
        .unwrap();

    // If the last step equals the most common step, return the most common step
    if (last_step - most_common_step).abs() < EPSILON {
        return Ok(most_common_step);
    }

    Err("wind_directions must be evenly spaced".to_string())
}

/// Reorders the wind directions so that they are adjacent.
/// Returns the reordered wind directions if the wind directions are not adjacent,
/// otherwise returns the input wind directions.
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

    // Get a list of steps
    let steps: Vec<f64> = wind_dirs_slice
        .windows(2)
        .map(|window| window[1] - window[0])
        .collect();

    // There will be at most one step with a size larger than the step size
    // If there is one, find it
    let large_step_idx = steps.iter().position(|&step| step > step_size + EPSILON);

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

/// Calculates the deviation from West (270 degrees) for wind directions.
/// The result is always positive values between 0 and 360.
///
/// | Input | Result |
/// | ----- | ------ |
/// | 270.0 | 0.0    |
/// | 280.0 | 10.0   |
/// | 360.0 | 90.0   |
/// | 180.0 | 270.0  |
/// | -10.0 | 80.0   |
/// |-100.0 | 350.0  |
///
/// # Arguments
///
/// * `wind_directions` - Array or scalar of wind directions
///
/// # Returns
///
/// The delta between the given wind direction(s) and 270, as positive value(s) between 0 and 360.
pub fn wind_delta(wind_directions: ArrayView1<f64>) -> Array1 {
    wind_directions.mapv(|dir| {
        let result = (dir - 270.0) % 360.0;
        if result < 0.0 {
            result + 360.0
        } else {
            result
        }
    })
}

/// Scalar version of wind_delta for single values
pub fn wind_delta_scalar(wind_direction: f64) -> f64 {
    let result = (wind_direction - 270.0) % 360.0;
    if result < 0.0 {
        result + 360.0
    } else {
        result
    }
}

// ============================================================================
// Coordinate Rotation Functions
// ============================================================================

/// Rotate coordinates relative to west (270 degrees).
///
/// This function rotates the given coordinates so that they are aligned with West (270) rather
/// than North (0). The rotation happens about the centroid of the coordinates.
///
/// # Arguments
///
/// * `wind_directions` - Series of wind directions to base the rotation
/// * `coordinates` - Series of coordinates to rotate with shape (N coordinates, 3)
///
/// # Returns
///
/// Tuple of (x_rotated, y_rotated, z_rotated, x_center, y_center)
pub fn rotate_coordinates_rel_west(
    wind_directions: &Array1,
    coordinates: &Array2,
) -> Result<(Array2, Array2, Array2, Array2, Array2)> {
    rotate_coordinates_rel_west_with_centers(wind_directions, coordinates, None, None)
}

/// Rotate coordinates to the wind-aligned (west) frame, optionally about an
/// explicit center (Python `rotate_coordinates_rel_west`'s
/// `x_center_of_rotation` / `y_center_of_rotation` arguments, utilities.py).
///
/// The center matters when mixing grids in one rotated frame: a `PointsGrid`
/// built over sample points must rotate about the SAME center as the farm
/// solve (`core.solve_for_points` passes the turbine grid's centers,
/// core.py:246-247), otherwise relative turbine/point positions drift.
/// `None` falls back to the bounding-box center of the given coordinates.
pub fn rotate_coordinates_rel_west_with_centers(
    wind_directions: &Array1,
    coordinates: &Array2,
    x_center_opt: Option<Float>,
    y_center_opt: Option<Float>,
) -> Result<(Array2, Array2, Array2, Array2, Array2)> {
    let n_findex = wind_directions.len();
    let n_coords = coordinates.shape()[0];

    // Construct the arrays storing the coordinates
    let x_coordinates = coordinates.column(0);
    let y_coordinates = coordinates.column(1);
    let z_coordinates = if coordinates.shape()[1] > 2 {
        coordinates.column(2).to_owned()
    } else {
        Array1::zeros(n_coords)
    };

    // Find center of rotation - this is the center of box bounding all coordinates
    let x_min = x_coordinates.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = x_coordinates
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = y_coordinates.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = y_coordinates
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let x_center = x_center_opt.unwrap_or((x_min + x_max) / 2.0);
    let y_center = y_center_opt.unwrap_or((y_min + y_max) / 2.0);

    // Calculate the difference in given wind direction from 270 / West
    let wind_deviation_from_west = wind_delta(wind_directions.view());

    let mut x_rotated = Array2::zeros((n_findex, n_coords));
    let mut y_rotated = Array2::zeros((n_findex, n_coords));
    let mut z_rotated = Array2::zeros((n_findex, n_coords));
    let mut x_center_array = Array2::zeros((n_findex, n_coords));
    let mut y_center_array = Array2::zeros((n_findex, n_coords));

    for fi in 0..n_findex {
        let angle = wind_deviation_from_west[fi];
        let cos_angle = cosd(angle);
        let sin_angle = sind(angle);

        for ci in 0..n_coords {
            let x_offset = coordinates[[ci, 0]] - x_center;
            let y_offset = coordinates[[ci, 1]] - y_center;

            x_rotated[[fi, ci]] = x_offset * cos_angle - y_offset * sin_angle + x_center;
            y_rotated[[fi, ci]] = x_offset * sin_angle + y_offset * cos_angle + y_center;
            z_rotated[[fi, ci]] = z_coordinates[ci];

            x_center_array[[fi, ci]] = x_center;
            y_center_array[[fi, ci]] = y_center;
        }
    }

    Ok((
        x_rotated,
        y_rotated,
        z_rotated,
        x_center_array,
        y_center_array,
    ))
}

/// Reverse rotate coordinates from wind frame back to inertial frame.
///
/// This function reverses the rotation of the given grid so that the coordinates are aligned with
/// the original wind direction. The rotation happens about the centroid of the coordinates.
///
/// # Arguments
///
/// * `wind_directions` - Series of wind directions to base the rotation
/// * `grid_x`, `grid_y`, `grid_z` - Grid coordinates in wind frame (n_findex, ...)
/// * `x_center_of_rotation` - The x-coordinate for the rotation center
/// * `y_center_of_rotation` - The y-coordinate for the rotational center
///
/// # Returns
///
/// Tuple of (x_reversed, y_reversed, z_reversed)
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

    let mut x_reversed = grid_x.to_owned();
    let mut y_reversed = grid_y.to_owned();
    let z_reversed = grid_z.clone();

    for fi in 0..n_findex {
        // Calculate the difference in given wind direction from 270 / West
        // We are rotating in the other direction
        let angle = -wind_delta_scalar(wind_directions[fi]);
        let cos_angle = cosd(angle);
        let sin_angle = sind(angle);

        // Get center of rotation for this findex (use first turbine's center)
        let x_center = x_center_of_rotation[[fi, 0]];
        let y_center = y_center_of_rotation[[fi, 0]];

        for ti in 0..shape[1] {
            for i in 0..shape[2] {
                for j in 0..shape[3] {
                    let x_rot = grid_x[[fi, ti, i, j]];
                    let y_rot = grid_y[[fi, ti, i, j]];

                    let x_offset = x_rot - x_center;
                    let y_offset = y_rot - y_center;

                    x_reversed[[fi, ti, i, j]] =
                        x_offset * cos_angle - y_offset * sin_angle + x_center;
                    y_reversed[[fi, ti, i, j]] =
                        x_offset * sin_angle + y_offset * cos_angle + y_center;
                }
            }
        }
    }

    Ok((x_reversed, y_reversed, z_reversed))
}

// ============================================================================
// YAML Loading Functions
// ============================================================================

/// Load YAML configuration from file
pub fn load_yaml<P: AsRef<Path>>(path: P) -> Result<Value> {
    let mut file = File::open(path.as_ref())
        .with_context(|| format!("Failed to open file: {:?}", path.as_ref()))?;

    let mut contents = String::new();
    file.read_to_string(&mut contents)
        .context("Failed to read file contents")?;

    let value: Value = serde_yaml::from_str(&contents).context("Failed to parse YAML")?;

    Ok(value)
}

// ============================================================================
// Rounding Functions
// ============================================================================

/// Rounds a number (with a 0.5 buffer) up to the nearest integer divisible by 2 or 5.
///
/// # Arguments
///
/// * `x` - The number to be rounded
///
/// # Returns
///
/// The rounded number
pub fn round_nearest_2_or_5(x: f64) -> i64 {
    let base_2: i64 = 2;
    let base_5: i64 = 5;
    let rounded_2 = base_2 * (((x + 0.5) / (base_2 as f64)).ceil() as i64);
    let rounded_5 = base_5 * (((x + 0.5) / (base_5 as f64)).ceil() as i64);
    rounded_2.min(rounded_5)
}

/// Rounds a number (with a 0.5 buffer) up to the nearest integer divisible by base.
///
/// # Arguments
///
/// * `x` - The number to be rounded
/// * `base` - The base to round to (default: 5)
///
/// # Returns
///
/// The rounded number
pub fn round_nearest(x: f64, base: i64) -> i64 {
    base * ((x + 0.5) / base as f64).ceil() as i64
}

// ============================================================================
// Nested Dictionary Operations
// ============================================================================

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

/// Print a nested dictionary with indentation (for debugging)
pub fn print_nested_dict(dictionary: &Value, indent: usize) {
    if let Some(map) = dictionary.as_mapping() {
        for (key, value) in map {
            println!("{}{:?}", " ".repeat(indent), key);
            if value.is_mapping() {
                print_nested_dict(value, indent + 4);
            } else {
                println!("{}{:?}", " ".repeat(indent + 4), value);
            }
        }
    }
}

/// Check if all values in a dictionary are scalar (not arrays or nested dicts)
pub fn is_all_scalar_dict(dictionary: &Value) -> bool {
    if let Some(map) = dictionary.as_mapping() {
        for (_, value) in map {
            if value.is_sequence() || value.is_mapping() {
                return false;
            }
        }
    }
    true
}

// ============================================================================
// Tests
// ============================================================================

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
    fn test_sind() {
        assert_relative_eq!(sind(0.0), 0.0);
        assert_relative_eq!(sind(90.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(sind(180.0), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_tand() {
        assert_relative_eq!(tand(0.0), 0.0);
        assert_relative_eq!(tand(45.0), 1.0, epsilon = 1e-10);
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
    fn test_check_step_size_valid_two_elements() {
        let directions = Array1::from_vec(vec![0.0, 90.0]);
        let result = check_and_identify_step_size(directions.view());
        assert_eq!(result.unwrap(), 90.0);
    }

    #[test]
    fn test_check_step_size_valid_three_elements() {
        let directions = Array1::from_vec(vec![0.0, 90.0, 180.0]);
        let result = check_and_identify_step_size(directions.view());
        assert_eq!(result.unwrap(), 90.0);
    }

    #[test]
    fn test_check_step_size_valid_many_elements() {
        let directions = Array1::from_vec(vec![
            0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0,
        ]);
        let result = check_and_identify_step_size(directions.view());
        assert_eq!(result.unwrap(), 30.0);
    }

    #[test]
    fn test_check_step_size_error_less_than_two() {
        let directions = Array1::from_vec(vec![0.0]);
        let result = check_and_identify_step_size(directions.view());
        assert!(result.is_err());
    }

    #[test]
    fn test_check_step_size_error_not_monotonic() {
        let directions = Array1::from_vec(vec![0.0, 90.0, 45.0]);
        let result = check_and_identify_step_size(directions.view());
        assert!(result.is_err());
    }

    #[test]
    fn test_check_step_size_error_not_evenly_spaced() {
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
        // Wind directions with a large gap (but still monotonically increasing)
        // The gap between 90 and 330 is 240, which is larger than the step size of 30
        let directions = Array1::from_vec(vec![0.0, 30.0, 60.0, 90.0, 330.0]);
        let (reordered, indices) = make_wind_directions_adjacent(directions.view());

        // After reordering: [330-360, 0, 30, 60, 90] = [-30, 0, 30, 60, 90]
        assert_eq!(
            reordered,
            Array1::from_vec(vec![-30.0, 0.0, 30.0, 60.0, 90.0])
        );
        assert_eq!(indices, vec![4, 0, 1, 2, 3]);
    }

    #[test]
    fn test_wind_delta() {
        let directions = Array1::from_vec(vec![270.0, 280.0, 360.0, 180.0]);
        let deltas = wind_delta(directions.view());

        assert_relative_eq!(deltas[0], 0.0);
        assert_relative_eq!(deltas[1], 10.0);
        assert_relative_eq!(deltas[2], 90.0);
        assert_relative_eq!(deltas[3], 270.0);
    }

    #[test]
    fn test_round_nearest_2_or_5() {
        // Test with Python: min(2*ceil((x+0.5)/2), 5*ceil((x+0.5)/5))
        assert_eq!(round_nearest_2_or_5(3.0), 4); // min(2*2, 5*1) = min(4, 5) = 4
        assert_eq!(round_nearest_2_or_5(6.0), 8); // min(2*4, 5*2) = min(8, 10) = 8
        assert_eq!(round_nearest_2_or_5(7.0), 8); // min(2*4, 5*2) = min(8, 10) = 8
        assert_eq!(round_nearest_2_or_5(11.0), 12); // min(2*6, 5*3) = min(12, 15) = 12
                                                    // For 4.0: ceil((4.0+0.5)/2) = ceil(2.25) = 3, so 2*3 = 6
                                                    //          ceil((4.0+0.5)/5) = ceil(0.9) = 1, so 5*1 = 5
                                                    //          min(6, 5) = 5
        assert_eq!(round_nearest_2_or_5(4.0), 5);
        // For 5.0: ceil((5.0+0.5)/2) = ceil(2.75) = 3, so 2*3 = 6
        //          ceil((5.0+0.5)/5) = ceil(1.1) = 2, so 5*2 = 10
        //          min(6, 10) = 6
        assert_eq!(round_nearest_2_or_5(5.0), 6);
    }

    #[test]
    fn test_round_nearest() {
        assert_eq!(round_nearest(3.0, 5), 5);
        assert_eq!(round_nearest(7.0, 5), 10);
        assert_eq!(round_nearest(12.0, 5), 15);
    }

    #[test]
    fn test_round_nearest_2_or_5_debug() {
        let x = 4.0_f64;
        let base_2: i64 = 2;
        let base_5: i64 = 5;

        let val1 = (x + 0.5) / (base_2 as f64);
        let ceil1 = val1.ceil();
        let as_i64_1 = ceil1 as i64;
        let rounded_2 = base_2 * as_i64_1;

        let val2 = (x + 0.5) / (base_5 as f64);
        let ceil2 = val2.ceil();
        let as_i64_2 = ceil2 as i64;
        let rounded_5 = base_5 * as_i64_2;

        eprintln!("Manual calculation:");
        eprintln!(
            "val1={}, ceil1={}, as_i64_1={}, rounded_2={}",
            val1, ceil1, as_i64_1, rounded_2
        );
        eprintln!(
            "val2={}, ceil2={}, as_i64_2={}, rounded_5={}",
            val2, ceil2, as_i64_2, rounded_5
        );
        eprintln!("min={}", rounded_2.min(rounded_5));

        let result = round_nearest_2_or_5(4.0);
        eprintln!("\nFunction result: {}", result);

        assert_eq!(result, 5);
    }
}
