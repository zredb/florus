use crate::types::Float;
use crate::floris_model::FlorisModel;
use crate::Array2;
use super::{YawOptimization, YawOptimizationResult, YawOptimizationConfig};

/// Configuration for SciPy optimizer
#[derive(Debug, Clone)]
pub struct SciPyOptimizationConfig {
    pub opt_method: String,
    pub maxiter: usize,
    pub ftol: Float,
    pub eps: Float,
}

impl Default for SciPyOptimizationConfig {
    fn default() -> Self {
        Self {
            opt_method: "SLSQP".to_string(),
            maxiter: 100,
            ftol: 1e-12,
            eps: 0.1,
        }
    }
}

/// SciPy-style yaw optimizer using coordinate descent
#[derive(Debug, Clone)]
pub struct YawOptimizationScipy {
    config: YawOptimizationConfig,
    scipy_config: SciPyOptimizationConfig,
}

impl YawOptimizationScipy {
    pub fn new() -> Self {
        Self {
            config: YawOptimizationConfig::default(),
            scipy_config: SciPyOptimizationConfig::default(),
        }
    }
}

impl Default for YawOptimizationScipy {
    fn default() -> Self {
        Self::new()
    }
}

impl YawOptimization for YawOptimizationScipy {
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

        let baseline_yaw_angles = self.config.yaw_angles_baseline.clone()
            .unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));

        fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
        fmodel.run()?;
        let baseline_power = fmodel.get_farm_power().sum();

        let mut yaw_angles = baseline_yaw_angles.clone();
        let step_size = self.scipy_config.eps;
        let tolerance = self.scipy_config.ftol;
        let max_iter = self.scipy_config.maxiter;

        let turbine_weights = self.config.turbine_weights.clone()
            .unwrap_or_else(|| vec![1.0; n_turbines]);

        let downstream_indices = if self.config.exclude_downstream_turbines {
            super::derive_downstream_turbines(fmodel, 0.3, false)?
        } else {
            vec![]
        };

        let active_turbines: Vec<usize> = (0..n_turbines)
            .filter(|ti| !downstream_indices.contains(ti))
            .collect();

        for _ in 0..max_iter {
            let mut max_change: Float = 0.0;

            for &ti in &active_turbines {
                for fi in 0..n_findex {
                    let current_yaw = yaw_angles[[fi, ti]];

                    let yaw_plus = (current_yaw + step_size)
                        .clamp(self.config.minimum_yaw_angle, self.config.maximum_yaw_angle);
                    let yaw_minus = (current_yaw - step_size)
                        .clamp(self.config.minimum_yaw_angle, self.config.maximum_yaw_angle);

                    yaw_angles[[fi, ti]] = yaw_plus;
                    fmodel.set_yaw_angles(yaw_angles.clone())?;
                    fmodel.run()?;
                    let power_plus = fmodel.get_farm_power()[fi];

                    yaw_angles[[fi, ti]] = yaw_minus;
                    fmodel.set_yaw_angles(yaw_angles.clone())?;
                    fmodel.run()?;
                    let power_minus = fmodel.get_farm_power()[fi];

                    let gradient = (power_plus - power_minus) / (2.0 * step_size);
                    let weight = turbine_weights[ti];

                    let new_yaw = (current_yaw - gradient * weight * step_size)
                        .clamp(self.config.minimum_yaw_angle, self.config.maximum_yaw_angle);

                    let change = (new_yaw - current_yaw).abs();
                    max_change = max_change.max(change);

                    yaw_angles[[fi, ti]] = new_yaw;
                }
            }

            if max_change < tolerance {
                break;
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
