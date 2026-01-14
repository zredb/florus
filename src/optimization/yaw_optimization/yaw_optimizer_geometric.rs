use crate::types::Float;
use crate::floris_model::FlorisModel;
use crate::Array2;
use super::{YawOptimization, YawOptimizationResult, YawOptimizationConfig, TrapezoidBounds};

/// Geometric yaw optimizer
///
/// Computes yaw angles based on farm geometry and wind direction.
/// Based on trapezoid wake model and turbine positioning.
#[derive(Debug, Clone)]
pub struct YawOptimizationGeometric {
    config: YawOptimizationConfig,
}

impl YawOptimizationGeometric {
    pub fn new() -> Self {
        Self {
            config: YawOptimizationConfig::default(),
        }
    }
}

impl Default for YawOptimizationGeometric {
    fn default() -> Self {
        Self::new()
    }
}

impl YawOptimization for YawOptimizationGeometric {
    fn optimize(
        &mut self,
        fmodel: &mut FlorisModel,
        config: Option<YawOptimizationConfig>,
    ) -> anyhow::Result<YawOptimizationResult> {
        if let Some(cfg) = config {
            self.config = cfg;
        }

        if !fmodel.state.initialized {
            anyhow::bail!("FlorisModel must be run before optimization");
        }

        let n_findex = fmodel.flow_field.n_findex;
        let n_turbines = fmodel.farm.n_turbines();
        let rotor_diameter = fmodel.farm.rotor_diameters[0];

        let layout_x = fmodel.farm.layout_x.clone();
        let layout_y = fmodel.farm.layout_y.clone();
        let wind_directions = fmodel.flow_field.wind_directions.clone();

        let mut yaw_angles = Array2::zeros((n_findex, n_turbines));
        let bounds = TrapezoidBounds::new();

        for fi in 0..n_findex {
            let wind_direction = wind_directions[fi];
            let turbine_x: Vec<Float> = layout_x.iter().copied().collect();
            let turbine_y: Vec<Float> = layout_y.iter().copied().collect();

            let yaw_for_direction = super::geometric_yaw(
                &turbine_x,
                &turbine_y,
                wind_direction,
                rotor_diameter,
                Some(bounds.clone()),
            );

            for ti in 0..n_turbines {
                yaw_angles[[fi, ti]] = yaw_for_direction[ti].clamp(
                    self.config.minimum_yaw_angle,
                    self.config.maximum_yaw_angle,
                );
            }
        }

        let baseline_yaw_angles = self.config.yaw_angles_baseline.clone()
            .unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));
        fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
        fmodel.run()?;
        let baseline_power = fmodel.get_farm_power().sum();

        fmodel.set_yaw_angles(yaw_angles.clone())?;
        fmodel.run()?;
        let optimized_power = fmodel.get_farm_power().sum();

        let power_improvement = optimized_power - baseline_power;
        let improvement_percentage = if baseline_power > 0.0 {
            (power_improvement / baseline_power) * 100.0
        } else {
            0.0
        };

        let turbine_powers = fmodel.get_turbine_powers();

        Ok(YawOptimizationResult {
            yaw_angles,
            powers: turbine_powers,
            baseline_power,
            optimized_power,
            power_improvement,
            improvement_percentage,
        })
    }
}
