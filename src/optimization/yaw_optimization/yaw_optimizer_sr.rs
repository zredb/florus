use crate::types::Float;
use crate::floris_model::FlorisModel;
use crate::Array2;
use super::{YawOptimization, YawOptimizationResult, YawOptimizationConfig};

/// Configuration for Serial Refine optimizer
#[derive(Debug, Clone)]
pub struct SerialRefineConfig {
    pub ny_passes: usize,
    pub eval_grid_resolution: usize,
    pub use_memory: bool,
    pub step_size: Float,
}

impl Default for SerialRefineConfig {
    fn default() -> Self {
        Self {
            ny_passes: 3,
            eval_grid_resolution: 5,
            use_memory: true,
            step_size: 5.0,
        }
    }
}

/// Serial Refine yaw optimizer
///
/// Optimizes turbines one at a time from front to back with multiple passes.
#[derive(Debug, Clone)]
pub struct YawOptimizationSR {
    config: YawOptimizationConfig,
    sr_config: SerialRefineConfig,
    turbine_order: Vec<usize>,
}

impl YawOptimizationSR {
    pub fn new() -> Self {
        Self {
            config: YawOptimizationConfig::default(),
            sr_config: SerialRefineConfig::default(),
            turbine_order: Vec::new(),
        }
    }

    fn calculate_turbine_order(fmodel: &FlorisModel) -> Vec<usize> {
        let n_turbines = fmodel.farm().n_turbines();
        let layout_x = fmodel.farm().layout_x.clone();

        let mut indices: Vec<(Float, usize)> = (0..n_turbines)
            .map(|i| (layout_x[i], i))
            .collect();

        indices.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        indices.into_iter().map(|(_, i)| i).collect()
    }
}

impl Default for YawOptimizationSR {
    fn default() -> Self {
        Self::new()
    }
}

impl YawOptimization for YawOptimizationSR {
    fn optimize(
        &mut self,
        fmodel: &mut FlorisModel,
        config: Option<YawOptimizationConfig>,
    ) -> anyhow::Result<YawOptimizationResult> {
        if let Some(cfg) = config {
            self.config = cfg;
        }

        if !fmodel.state().initialized {
            anyhow::bail!("FlorisModel must be run before optimization");
        }

        let n_findex = fmodel.flow_field().n_findex;
        let n_turbines = fmodel.farm().n_turbines();

        self.turbine_order = Self::calculate_turbine_order(fmodel);

        let baseline_yaw_angles = self.config.yaw_angles_baseline.clone()
            .unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));

        fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
        fmodel.run()?;
        let _baseline_power = fmodel.get_farm_power().sum();

        let mut yaw_angles = baseline_yaw_angles.clone();

        for fi in 0..n_findex {
            for _ in 0..self.sr_config.ny_passes {
                for &ti in &self.turbine_order {
                    let yaw_range = self.config.maximum_yaw_angle - self.config.minimum_yaw_angle;
                    let step = yaw_range / self.sr_config.eval_grid_resolution as Float;

                    let mut best_power = fmodel.get_farm_power()[fi];
                    let mut best_yaw = yaw_angles[[fi, ti]];

                    for i in 0..=self.sr_config.eval_grid_resolution {
                        let test_yaw = self.config.minimum_yaw_angle + i as Float * step;
                        yaw_angles[[fi, ti]] = test_yaw;

                        fmodel.set_yaw_angles(yaw_angles.clone())?;
                        fmodel.run()?;

                        let power = fmodel.get_farm_power()[fi];
                        if power > best_power {
                            best_power = power;
                            best_yaw = test_yaw;
                        }
                    }

                    yaw_angles[[fi, ti]] = best_yaw;
                }
            }
        }

        fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
        fmodel.run()?;
        let baseline_power_verify = fmodel.get_farm_power().sum();

        fmodel.set_yaw_angles(yaw_angles.clone())?;
        fmodel.run()?;
        let optimized_power = fmodel.get_farm_power().sum();

        let power_improvement = optimized_power - baseline_power_verify;
        let improvement_percentage = if baseline_power_verify > 0.0 {
            (power_improvement / baseline_power_verify) * 100.0
        } else {
            0.0
        };

        let turbine_powers = fmodel.get_turbine_powers();

        Ok(YawOptimizationResult {
            yaw_angles,
            powers: turbine_powers,
            baseline_power: baseline_power_verify,
            optimized_power,
            power_improvement,
            improvement_percentage,
        })
    }
}
