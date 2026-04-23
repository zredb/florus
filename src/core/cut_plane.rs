//! CutPlane module for FLORUS
//! 
//! This module provides the CutPlane structure and related functions
//! for extracting 2D slices from 3D flow fields.
//! Corresponds to `floris/cut_plane.py` in Python FLORIS.

use crate::types::{Array1, Array2, Float};
use ndarray::Array;

/// Represents a 2D slice through a 3D flow field
#[derive(Debug, Clone)]
pub struct CutPlane {
    /// DataFrame-like structure with columns: x1, x2, x3, u, v, w
    pub data: CutPlaneData,
    /// Normal vector direction ("x", "y", or "z")
    pub normal_vector: String,
    /// Resolution in (x1, x2) directions
    pub resolution: (usize, usize),
}

/// Data structure for cut plane points
#[derive(Debug, Clone)]
pub struct CutPlaneData {
    /// x1 coordinates
    pub x1: Array1,
    /// x2 coordinates
    pub x2: Array1,
    /// x3 coordinates (constant for a plane)
    pub x3: Array1,
    /// u velocity component
    pub u: Array1,
    /// v velocity component
    pub v: Array1,
    /// w velocity component
    pub w: Array1,
}

impl CutPlane {
    /// Create a new CutPlane from data arrays
    pub fn new(
        x1: Array1,
        x2: Array1,
        x3: Array1,
        u: Array1,
        v: Array1,
        w: Array1,
        normal_vector: &str,
        resolution: (usize, usize),
    ) -> Self {
        Self {
            data: CutPlaneData { x1, x2, x3, u, v, w },
            normal_vector: normal_vector.to_string(),
            resolution,
        }
    }

    /// Get unique x1 values
    pub fn unique_x1(&self) -> Vec<f64> {
        let mut values: Vec<f64> = self.data.x1.iter().cloned().collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        values.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
        values
    }

    /// Get unique x2 values
    pub fn unique_x2(&self) -> Vec<f64> {
        let mut values: Vec<f64> = self.data.x2.iter().cloned().collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        values.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
        values
    }

    /// Set origin of the cut plane
    pub fn set_origin(mut self, center_x1: f64, center_x2: f64) -> Self {
        self.data.x1 = self.data.x1 - center_x1;
        self.data.x2 = self.data.x2 - center_x2;
        self
    }

    /// Rescale axis
    pub fn rescale_axis(mut self, x1_factor: f64, x2_factor: f64) -> Self {
        self.data.x1 = self.data.x1 / x1_factor;
        self.data.x2 = self.data.x2 / x2_factor;
        self
    }
}

/// Extract a horizontal plane (z-normal) from flow field data at hub height
/// 
/// This is a simplified version that extracts data at a specific z-level.
/// Full implementation would support interpolation between grid points.
pub fn extract_horizontal_plane(
    u_field: &Array<f64, ndarray::Ix4>,
    v_field: &Array<f64, ndarray::Ix4>,
    w_field: &Array<f64, ndarray::Ix4>,
    x_coords: &Array<f64, ndarray::Ix2>,
    y_coords: &Array<f64, ndarray::Ix2>,
    z_coords: &Array<f64, ndarray::Ix2>,
    findex: usize,
    z_value: f64,
) -> Option<CutPlane> {
    // Find the closest z-index to the requested z_value
    let z_shape = z_coords.shape();
    if z_shape.len() < 2 {
        return None;
    }
    
    let n_y = z_shape[0];
    let n_z = z_shape[1];
    
    // For simplicity, use the middle turbine's hub height grid
    // In full implementation, would search across all turbines
    let mut closest_z_idx = 0;
    let mut min_diff = f64::INFINITY;
    
    for k in 0..n_z {
        let diff = (z_coords[[0, k]] - z_value).abs();
        if diff < min_diff {
            min_diff = diff;
            closest_z_idx = k;
        }
    }
    
    // Extract the slice at this z-level for all y points
    let mut x1_vals = Vec::new();
    let mut x2_vals = Vec::new();
    let mut x3_vals = Vec::new();
    let mut u_vals = Vec::new();
    let mut v_vals = Vec::new();
    let mut w_vals = Vec::new();
    
    // Average over turbines (for now, just use first turbine)
    let turbine_idx = 0;
    
    for j in 0..n_y {
        let x = x_coords[[j, closest_z_idx]];
        let y = y_coords[[j, closest_z_idx]];
        let z = z_coords[[j, closest_z_idx]];
        
        // Get velocity at this point (averaged over y-z plane for this turbine)
        // Note: This is simplified - proper implementation would handle the 4D array correctly
        let u = if u_field.shape().len() >= 4 && turbine_idx < u_field.shape()[1] {
            u_field[[findex, turbine_idx, j, closest_z_idx]]
        } else {
            0.0
        };
        
        let v = if v_field.shape().len() >= 4 && turbine_idx < v_field.shape()[1] {
            v_field[[findex, turbine_idx, j, closest_z_idx]]
        } else {
            0.0
        };
        
        let w = if w_field.shape().len() >= 4 && turbine_idx < w_field.shape()[1] {
            w_field[[findex, turbine_idx, j, closest_z_idx]]
        } else {
            0.0
        };
        
        x1_vals.push(x);
        x2_vals.push(y);
        x3_vals.push(z);
        u_vals.push(u);
        v_vals.push(v);
        w_vals.push(w);
    }
    
    let resolution = (n_y, 1); // Simplified - only one z-level
    
    Some(CutPlane::new(
        Array::from_vec(x1_vals),
        Array::from_vec(x2_vals),
        Array::from_vec(x3_vals),
        Array::from_vec(u_vals),
        Array::from_vec(v_vals),
        Array::from_vec(w_vals),
        "z",
        resolution,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr1;
    
    #[test]
    fn test_cutplane_creation() {
        let x1 = arr1(&[0.0, 1.0, 2.0]);
        let x2 = arr1(&[0.0, 0.0, 0.0]);
        let x3 = arr1(&[90.0, 90.0, 90.0]);
        let u = arr1(&[8.0, 7.5, 7.0]);
        let v = arr1(&[0.0, 0.1, 0.2]);
        let w = arr1(&[0.0, 0.0, 0.0]);
        
        let cp = CutPlane::new(x1.clone(), x2.clone(), x3.clone(), 
                               u.clone(), v.clone(), w.clone(), "z", (3, 1));
        
        assert_eq!(cp.normal_vector, "z");
        assert_eq!(cp.resolution, (3, 1));
        assert_eq!(cp.data.x1.len(), 3);
    }
    
    #[test]
    fn test_set_origin() {
        let x1 = arr1(&[10.0, 20.0, 30.0]);
        let x2 = arr1(&[5.0, 5.0, 5.0]);
        let x3 = arr1(&[90.0, 90.0, 90.0]);
        let u = arr1(&[8.0, 8.0, 8.0]);
        let v = arr1(&[0.0, 0.0, 0.0]);
        let w = arr1(&[0.0, 0.0, 0.0]);
        
        let cp = CutPlane::new(x1, x2, x3, u, v, w, "z", (3, 1));
        let cp_shifted = cp.set_origin(10.0, 5.0);
        
        assert!((cp_shifted.data.x1[0] - 0.0).abs() < 1e-10);
        assert!((cp_shifted.data.x2[0] - 0.0).abs() < 1e-10);
    }
}
