/// Utility functions for FLORIS-RS
/// 
/// Corresponds to utilities.py in the Python implementation

use crate::types::{Float, Array1, Array2, Array4};
use anyhow::{Context, Result};
use serde_yaml::Value;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use ndarray::Array;

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
}
