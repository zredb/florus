/// Turbine operation models
/// 
/// Different control modes for turbine operation

use crate::types::Float;

pub const POWER_SETPOINT_DEFAULT: Float = 1e12;
pub const POWER_SETPOINT_DISABLED: Float = 1e12;

/// Simple derating operation model
#[derive(Debug, Clone)]
pub struct SimpleDeratingModel {
    pub power_setpoint: Float,
}

impl SimpleDeratingModel {
    pub fn new(power_setpoint: Float) -> Self {
        Self { power_setpoint }
    }
    
    /// Calculate derated power
    pub fn calculate_power(&self, available_power: Float) -> Float {
        if self.power_setpoint >= POWER_SETPOINT_DISABLED {
            available_power
        } else {
            available_power.min(self.power_setpoint)
        }
    }
}

/// Cosine loss model for yaw misalignment
pub fn cosine_loss_model(yaw_angle: Float, p_p: Float) -> Float {
    let yaw_rad = yaw_angle.to_radians();
    yaw_rad.cos().powf(p_p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_simple_derating() {
        let model = SimpleDeratingModel::new(3e6);
        
        // Available power below setpoint
        assert_relative_eq!(model.calculate_power(2e6), 2e6);
        
        // Available power above setpoint
        assert_relative_eq!(model.calculate_power(4e6), 3e6);
        
        // Disabled setpoint
        let model_disabled = SimpleDeratingModel::new(POWER_SETPOINT_DISABLED);
        assert_relative_eq!(model_disabled.calculate_power(5e6), 5e6);
    }
    
    #[test]
    fn test_cosine_loss() {
        // No yaw misalignment
        assert_relative_eq!(cosine_loss_model(0.0, 1.88), 1.0);
        
        // 25 degree yaw
        let loss = cosine_loss_model(25.0, 1.88);
        assert!(loss < 1.0);
        assert!(loss > 0.5);
    }
}
