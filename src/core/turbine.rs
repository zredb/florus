/// Turbine model and operations
///
/// Corresponds to core/turbine/ module in Python implementation
use crate::types::{Float, Array2, Array4};
use crate::core::rotor_velocity::{rotor_effective_velocity, AveragingMethod};
use ndarray::Array;

/// Wind turbine representation with power and thrust tables
#[derive(Debug, Clone)]
pub struct Turbine {
    /// Turbine type name
    pub turbine_type: String,
    
    /// Hub height [m]
    pub hub_height: Float,
    
    /// Rotor diameter [m]
    pub rotor_diameter: Float,
    
    /// Tip-speed ratio
    pub tsr: Float,
    
    /// Reference tilt angle [degrees]
    pub ref_tilt: Float,
    
    /// Rated power [W]
    pub rated_power: Float,
    
    /// Cut-in wind speed [m/s]
    pub cut_in_wind_speed: Float,
    
    /// Cut-out wind speed [m/s]
    pub cut_out_wind_speed: Float,
    
    /// Rated wind speed [m/s]
    pub rated_wind_speed: Float,
    
    /// Power curve - wind speeds [m/s]
    pub power_curve_wind_speeds: Vec<Float>,
    
    /// Power curve - power values [W]
    pub power_curve_powers: Vec<Float>,
    
    /// Thrust coefficient curve - wind speeds [m/s]
    pub thrust_coefficient_wind_speeds: Vec<Float>,
    
    /// Thrust coefficient curve - Ct values
    pub thrust_coefficient_values: Vec<Float>,
    
    /// Operation model type
    pub operation_model: String,
}

impl Turbine {
    /// Calculate power output for given velocities
    /// 
    /// # Arguments
    /// * `velocities` - Wind velocities at turbine rotor (findex, n_turbines, n_points, n_grid)
    /// * `air_density` - Air density [kg/m³]
    /// * `yaw_angles` - Yaw angles [degrees]
    /// * `tilt_angles` - Tilt angles [degrees]
    /// * `average_method` - Method for averaging velocities
    pub fn power(
        &self,
        velocities: &Array4,
        air_density: Float,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        average_method: AveragingMethod,
    ) -> crate::Result<Array2> {
        // Calculate rotor effective velocity
        let rotor_velocities = rotor_effective_velocity(
            velocities,
            average_method,
            None, // cubature_weights
            yaw_angles,
            tilt_angles,
            None, // ref_tilt
            None, // cosine_loss_exponent_yaw
            None, // cosine_loss_exponent_tilt
            None, // correct_cp_ct_for_tilt
        )?;
        
        let shape = rotor_velocities.shape();
        let mut power_output = Array::zeros((shape[0], shape[1]));
        
        // Calculate power for each turbine at each condition
        for fi in 0..shape[0] {
            for ti in 0..shape[1] {
                let v = rotor_velocities[[fi, ti, 0]];
                
                if v < self.cut_in_wind_speed || v > self.cut_out_wind_speed {
                    power_output[[fi, ti]] = 0.0;
                } else {
                    let cp = self.power_coefficient(v);
                    // P = 0.5 * ρ * A * v³ * Cp
                    let area = std::f64::consts::PI * (self.rotor_diameter / 2.0).powi(2);
                    power_output[[fi, ti]] = 0.5 * air_density * area * v.powi(3) * cp;
                }
            }
        }
        
        Ok(power_output)
    }
    
    /// Calculate thrust coefficient for given velocities
    pub fn thrust_coefficient(
        &self,
        velocities: &Array4,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        average_method: AveragingMethod,
    ) -> crate::Result<Array2> {
        let rotor_velocities = rotor_effective_velocity(
            velocities,
            average_method,
            None,
            yaw_angles,
            tilt_angles,
            None,
            None,
            None,
            None,
        )?;
        
        let shape = rotor_velocities.shape();
        let mut ct_output = Array::zeros((shape[0], shape[1]));
        
        for fi in 0..shape[0] {
            for ti in 0..shape[1] {
                let v = rotor_velocities[[fi, ti, 0]];
                ct_output[[fi, ti]] = self.ct_at_speed(v);
            }
        }
        
        Ok(ct_output)
    }
    
    /// Get power coefficient at given wind speed via interpolation
    pub fn power_coefficient(&self, wind_speed: Float) -> Float {
        if wind_speed < self.cut_in_wind_speed || wind_speed > self.cut_out_wind_speed {
            return 0.0;
        }
        
        let power = self.interpolate(
            &self.power_curve_wind_speeds,
            &self.power_curve_powers,
            wind_speed,
        );
        
        // Convert power to Cp: Cp = P / (0.5 * ρ * A * v³)
        // Using standard air density
        let area = std::f64::consts::PI * (self.rotor_diameter / 2.0).powi(2);
        let rho = 1.225;
        let max_power = 0.5 * rho * area * wind_speed.powi(3);
        
        if max_power > 0.0 {
            (power / max_power).min(0.59) // Betz limit
        } else {
            0.0
        }
    }
    
    /// Get thrust coefficient at given wind speed
    pub fn ct_at_speed(&self, wind_speed: Float) -> Float {
        if wind_speed < self.cut_in_wind_speed || wind_speed > self.cut_out_wind_speed {
            return 0.0;
        }
        
        self.interpolate(
            &self.thrust_coefficient_wind_speeds,
            &self.thrust_coefficient_values,
            wind_speed,
        )
    }
    
    /// Calculate axial induction factor from thrust coefficient
    pub fn axial_induction(&self, ct: Float) -> Float {
        // Using momentum theory relationship
        if ct < 0.96 {
            0.5 * (1.0 - (1.0 - ct).sqrt())
        } else {
            // High thrust region - empirical relationship
            0.143 + (0.0203 - 0.6427 * (0.889 - ct)).sqrt()
        }
    }
    
    /// Linear interpolation helper
    fn interpolate(&self, x_values: &[Float], y_values: &[Float], x: Float) -> Float {
        let n = x_values.len();
        
        // Handle edge cases
        if n == 0 {
            return 0.0;
        }
        if x <= x_values[0] {
            return y_values[0];
        }
        if x >= x_values[n - 1] {
            return y_values[n - 1];
        }
        
        // Find bracketing indices
        for i in 0..n - 1 {
            if x >= x_values[i] && x <= x_values[i + 1] {
                let x0 = x_values[i];
                let x1 = x_values[i + 1];
                let y0 = y_values[i];
                let y1 = y_values[i + 1];
                
                // Linear interpolation
                let t = (x - x0) / (x1 - x0);
                return y0 + t * (y1 - y0);
            }
        }
        
        // Should not reach here
        y_values[n - 1]
    }
}

/// Calculate power for array of turbines
pub fn power(
    velocities: &Array4,
    turbines: &[Turbine],
    air_density: Float,
    yaw_angles: Option<&Array2>,
    tilt_angles: Option<&Array2>,
    average_method: AveragingMethod,
) -> crate::Result<Array2> {
    let shape = velocities.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];
    
    let mut power_output = Array::zeros((n_findex, n_turbines));
    
    for ti in 0..n_turbines {
        if ti < turbines.len() {
            let turbine_velocities = velocities.slice(ndarray::s![.., ti..(ti+1), .., ..]).to_owned();
            let turbine_power = turbines[ti].power(
                &turbine_velocities,
                air_density,
                yaw_angles,
                tilt_angles,
                average_method,
            )?;
            
            for fi in 0..n_findex {
                power_output[[fi, ti]] = turbine_power[[fi, 0]];
            }
        }
    }
    
    Ok(power_output)
}

/// Calculate thrust coefficient for array of turbines
pub fn thrust_coefficient(
    velocities: &Array4,
    turbines: &[Turbine],
    yaw_angles: Option<&Array2>,
    tilt_angles: Option<&Array2>,
    average_method: AveragingMethod,
) -> crate::Result<Array2> {
    let shape = velocities.shape();
    let n_findex = shape[0];
    let n_turbines = shape[1];
    
    let mut ct_output = Array::zeros((n_findex, n_turbines));
    
    for ti in 0..n_turbines {
        if ti < turbines.len() {
            let turbine_velocities = velocities.slice(ndarray::s![.., ti..(ti+1), .., ..]).to_owned();
            let turbine_ct = turbines[ti].thrust_coefficient(
                &turbine_velocities,
                yaw_angles,
                tilt_angles,
                average_method,
            )?;
            
            for fi in 0..n_findex {
                ct_output[[fi, ti]] = turbine_ct[[fi, 0]];
            }
        }
    }
    
    Ok(ct_output)
}

/// Calculate axial induction from thrust coefficient
pub fn axial_induction(
    velocities: &Array4,
    turbines: &[Turbine],
    yaw_angles: Option<&Array2>,
    tilt_angles: Option<&Array2>,
    average_method: AveragingMethod,
) -> crate::Result<Array2> {
    let ct = thrust_coefficient(velocities, turbines, yaw_angles, tilt_angles, average_method)?;
    
    let mut ai = Array::zeros(ct.dim());
    
    for ((i, j), &ct_val) in ct.indexed_iter() {
        if j < turbines.len() {
            ai[[i, j]] = turbines[j].axial_induction(ct_val);
        }
    }
    
    Ok(ai)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    fn create_test_turbine() -> Turbine {
        Turbine {
            turbine_type: "test_turbine".to_string(),
            hub_height: 90.0,
            rotor_diameter: 126.0,
            tsr: 8.0,
            ref_tilt: 5.0,
            rated_power: 5_000_000.0,
            cut_in_wind_speed: 3.0,
            cut_out_wind_speed: 25.0,
            rated_wind_speed: 11.0,
            power_curve_wind_speeds: vec![0.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 25.0, 30.0],
            power_curve_powers: vec![
                0.0, 0.0, 100_000.0, 500_000.0, 1_500_000.0, 
                5_000_000.0, 5_000_000.0, 5_000_000.0, 5_000_000.0, 0.0
            ],
            thrust_coefficient_wind_speeds: vec![0.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 25.0, 30.0],
            thrust_coefficient_values: vec![0.0, 0.8, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.0],
            operation_model: "simple".to_string(),
        }
    }
    
    #[test]
    fn test_turbine_power_coefficient() {
        let turbine = create_test_turbine();
        
        // Below cut-in
        assert_eq!(turbine.power_coefficient(2.0), 0.0);
        
        // Above cut-out
        assert_eq!(turbine.power_coefficient(26.0), 0.0);
        
        // In operating range
        let cp = turbine.power_coefficient(8.0);
        assert!(cp > 0.0 && cp <= 0.59);
    }
    
    #[test]
    fn test_turbine_ct() {
        let turbine = create_test_turbine();
        
        assert_eq!(turbine.ct_at_speed(2.0), 0.0);
        assert_eq!(turbine.ct_at_speed(26.0), 0.0);
        
        let ct = turbine.ct_at_speed(8.0);
        assert!(ct > 0.0 && ct < 1.0);
    }
    
    #[test]
    fn test_axial_induction() {
        let turbine = create_test_turbine();
        
        // Low Ct regime
        let ai_low = turbine.axial_induction(0.5);
        assert!(ai_low > 0.0 && ai_low < 0.5);
        
        // High Ct regime
        let ai_high = turbine.axial_induction(0.99);
        assert!(ai_high > 0.0 && ai_high < 1.0);
    }
    
    #[test]
    fn test_interpolation() {
        let turbine = create_test_turbine();
        
        // Test interpolation
        let val = turbine.interpolate(
            &vec![0.0, 10.0, 20.0],
            &vec![0.0, 100.0, 200.0],
            15.0
        );
        assert_relative_eq!(val, 150.0);
        
        // Test edge cases
        let val_low = turbine.interpolate(
            &vec![10.0, 20.0],
            &vec![100.0, 200.0],
            5.0
        );
        assert_eq!(val_low, 100.0);
        
        let val_high = turbine.interpolate(
            &vec![10.0, 20.0],
            &vec![100.0, 200.0],
            25.0
        );
        assert_eq!(val_high, 200.0);
    }
}
