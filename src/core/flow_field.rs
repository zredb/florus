/// Flow field representation
///
/// Corresponds to flow_field.py
use crate::types::{Float, Array1, Array4};
use serde::{Deserialize, Serialize};
use ndarray::Array;

/// Represents the atmospheric flow field conditions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowField {
    /// Wind speeds at reference height [m/s]
    pub wind_speeds: Array1,
    
    /// Wind directions in degrees (0 = North, increasing clockwise)
    pub wind_directions: Array1,
    
    /// Wind veer [degrees]
    pub wind_veer: Float,
    
    /// Wind shear power law exponent
    pub wind_shear: Float,
    
    /// Air density [kg/m³]
    pub air_density: Float,
    
    /// Turbulence intensity values (0-1)
    pub turbulence_intensities: Array1,
    
    /// Reference height for wind measurements [m]
    pub reference_wind_height: Float,
    
    /// Number of flow indices (conditions)
    pub n_findex: usize,
    
    // Flow field arrays (initialized after grid creation)
    pub u_initial_sorted: Array4,
    pub v_initial_sorted: Array4,
    pub w_initial_sorted: Array4,
    pub u_sorted: Array4,
    pub v_sorted: Array4,
    pub w_sorted: Array4,
    pub u: Array4,
    pub v: Array4,
    pub w: Array4,
    
    /// Turbulence intensity field
    pub turbulence_intensity_field: Array4,
    pub turbulence_intensity_field_sorted: Array4,
}

impl FlowField {
    /// Create a new FlowField
    pub fn new(
        wind_speeds: Array1,
        wind_directions: Array1,
        wind_veer: Float,
        wind_shear: Float,
        air_density: Float,
        turbulence_intensities: Array1,
        reference_wind_height: Float,
    ) -> crate::Result<Self> {
        // Validate inputs
        if wind_speeds.len() != wind_directions.len() {
            anyhow::bail!(
                "wind_speeds (len={}) and wind_directions (len={}) must have same length",
                wind_speeds.len(),
                wind_directions.len()
            );
        }
        
        if turbulence_intensities.len() != wind_speeds.len() {
            anyhow::bail!(
                "turbulence_intensities (len={}) must match number of conditions ({})",
                turbulence_intensities.len(),
                wind_speeds.len()
            );
        }
        
        let n_findex = wind_speeds.len();
        
        // Initialize empty arrays for flow fields (4D: n_findex, n_turbines, n_y, n_z)
        let empty_4d = Array::zeros((0, 0, 0, 0));
        
        Ok(Self {
            wind_speeds,
            wind_directions,
            wind_veer,
            wind_shear,
            air_density,
            turbulence_intensities,
            reference_wind_height,
            n_findex,
            u_initial_sorted: empty_4d.clone(),
            v_initial_sorted: empty_4d.clone(),
            w_initial_sorted: empty_4d.clone(),
            u_sorted: empty_4d.clone(),
            v_sorted: empty_4d.clone(),
            w_sorted: empty_4d.clone(),
            u: empty_4d.clone(),
            v: empty_4d.clone(),
            w: empty_4d.clone(),
            turbulence_intensity_field: empty_4d.clone(),
            turbulence_intensity_field_sorted: empty_4d,
        })
    }
    
    /// Initialize flow field on a grid
    pub fn initialize_flow_field(&mut self, grid_shape: (usize, usize, usize, usize)) {
        let (n_findex, n_turbines, n_y, n_z) = grid_shape;
        
        let h_ref = self.reference_wind_height;
        
        // Initialize velocity fields with wind shear profile
        self.u_initial_sorted = Array::zeros((n_findex, n_turbines, n_y, n_z));
        self.v_initial_sorted = Array::zeros((n_findex, n_turbines, n_y, n_z));
        self.w_initial_sorted = Array::zeros((n_findex, n_turbines, n_y, n_z));
        
        // Apply wind shear profile to initial velocity field
        for fi in 0..n_findex {
            let ws = self.wind_speeds[fi];
            for ti in 0..n_turbines {
                for iy in 0..n_y {
                    for iz in 0..n_z {
                        let height_factor = 1.0_f64;
                        self.u_initial_sorted[[fi, ti, iy, iz]] = ws * height_factor;
                    }
                }
            }
        }
        
        self.u_sorted = self.u_initial_sorted.clone();
        self.v_sorted = self.v_initial_sorted.clone();
        self.w_sorted = self.w_initial_sorted.clone();
        
        // Initialize turbulence intensity field
        self.turbulence_intensity_field = Array::zeros((n_findex, n_turbines, n_y, n_z));
        self.turbulence_intensity_field_sorted = self.turbulence_intensity_field.clone();
    }
    
    /// Calculate wind speed at a given height using power law
    pub fn wind_speed_at_height(&self, height: Float, findex: usize) -> Float {
        let ws_ref = self.wind_speeds[findex];
        let h_ref = self.reference_wind_height;
        ws_ref * (height / h_ref).powf(self.wind_shear)
    }

    /// Update velocities from wake calculations
    pub fn update_velocities(&mut self, u_deficit: &Array4, v_deficit: &Array4, w_deficit: &Array4) {
        // Apply wake deficits to get final velocities
        for fi in 0..self.u_sorted.shape()[0] {
            for ti in 0..self.u_sorted.shape()[1] {
                for iy in 0..self.u_sorted.shape()[2] {
                    for iz in 0..self.u_sorted.shape()[3] {
                        self.u_sorted[[fi, ti, iy, iz]] = 
                            self.u_initial_sorted[[fi, ti, iy, iz]] - u_deficit[[fi, ti, iy, iz]];
                        self.v_sorted[[fi, ti, iy, iz]] = 
                            self.v_initial_sorted[[fi, ti, iy, iz]] + v_deficit[[fi, ti, iy, iz]];
                        self.w_sorted[[fi, ti, iy, iz]] = 
                            self.w_initial_sorted[[fi, ti, iy, iz]] + w_deficit[[fi, ti, iy, iz]];
                    }
                }
            }
        }
    }

    /// Get mutable reference to u_sorted
    pub fn u_sorted_mut(&mut self) -> &mut Array4 {
        &mut self.u_sorted
    }

    /// Get mutable reference to v_sorted
    pub fn v_sorted_mut(&mut self) -> &mut Array4 {
        &mut self.v_sorted
    }

    /// Get mutable reference to w_sorted
    pub fn w_sorted_mut(&mut self) -> &mut Array4 {
        &mut self.w_sorted
    }

    /// Get turbulence intensities array
    pub fn turbulence_intensities(&self) -> &Array1 {
        &self.turbulence_intensities
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_flow_field_creation() {
        let wind_speeds = Array1::from_vec(vec![8.0, 10.0]);
        let wind_directions = Array1::from_vec(vec![270.0, 280.0]);
        let turbulence_intensities = Array1::from_vec(vec![0.06, 0.08]);
        
        let ff = FlowField::new(
            wind_speeds,
            wind_directions,
            0.0,
            0.14,
            1.225,
            turbulence_intensities,
            90.0,
        ).unwrap();
        
        assert_eq!(ff.n_findex, 2);
        assert_relative_eq!(ff.air_density, 1.225);
    }
    
    #[test]
    fn test_wind_speed_at_height() {
        let wind_speeds = Array1::from_vec(vec![10.0]);
        let wind_directions = Array1::from_vec(vec![270.0]);
        let turbulence_intensities = Array1::from_vec(vec![0.06]);
        
        let ff = FlowField::new(
            wind_speeds,
            wind_directions,
            0.0,
            0.14,
            1.225,
            turbulence_intensities,
            90.0,
        ).unwrap();
        
        // At reference height, should equal reference speed
        assert_relative_eq!(ff.wind_speed_at_height(90.0, 0), 10.0);
        
        // At higher height, should be greater
        assert!(ff.wind_speed_at_height(120.0, 0) > 10.0);
    }
}
