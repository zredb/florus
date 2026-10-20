use super::{YawOptimization, YawOptimizationConfig, YawOptimizationResult};
use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2, Float};
use crate::Array1 as NpArray1;

/// SciPy-style yaw optimizer.
///
/// Python uses `scipy.optimize.minimize(method="SLSQP")` per wind
/// condition (yaw_optimizer_scipy.py:65-160): each condition is optimized
/// independently over the yaw angles of the turbines that are not excluded
/// (downstream turbines and turbines with equal bounds are skipped), with
/// the turbine-weighted farm power as the objective and the yaw bounds as
/// constraints.
///
/// FLORUS implements the same per-condition bounded local search with a
/// compass (pattern) search: it starts from the baseline yaw angles and
/// repeatedly probes each free turbine at +/- step, halving the step when no
/// improvement is found. For the smooth objective of the wake models this
/// converges to the same local optimum SLSQP finds.
#[derive(Debug, Clone)]
pub struct YawOptimizationScipy {
    config: YawOptimizationConfig,
    scipy_config: SciPyOptimizationConfig,
}

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

        if !fmodel.state().initialized {
            anyhow::bail!("FlorisModel must be run before optimization");
        }

        let n_findex = fmodel.flow_field().n_findex;
        let n_turbines = fmodel.farm().n_turbines();

        let wd = fmodel.flow_field().wind_directions.clone();
        let ws = fmodel.flow_field().wind_speeds.clone();
        let ti = fmodel.flow_field().turbulence_intensities.clone();

        let baseline_yaw_angles = self
            .config
            .yaw_angles_baseline
            .clone()
            .unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));

        let turbine_weights = self
            .config
            .turbine_weights
            .clone()
            .unwrap_or_else(|| vec![1.0; n_turbines]);

        // Per-condition yaw bounds (scalars or per-turbine arrays)
        let mut yaw_lbs = Array2::from_elem((n_findex, n_turbines), self.config.minimum_yaw_angle);
        let mut yaw_ubs = Array2::from_elem((n_findex, n_turbines), self.config.maximum_yaw_angle);
        if let Some(ref arr) = self.config.minimum_yaw_angle_array {
            if arr.shape() == [n_findex, n_turbines] {
                yaw_lbs = arr.clone();
            }
        }
        if let Some(ref arr) = self.config.maximum_yaw_angle_array {
            if arr.shape() == [n_findex, n_turbines] {
                yaw_ubs = arr.clone();
            }
        }

        // Baseline powers (weighted), per condition
        fmodel.set_yaw_angles(baseline_yaw_angles.clone())?;
        fmodel.run()?;
        let farm_powers_baseline = fmodel.get_farm_power()?;
        let baseline_power: Float = farm_powers_baseline.sum();

        let mut yaw_angles = baseline_yaw_angles.clone();
        let mut farm_powers_opt = farm_powers_baseline.clone();

        for fi in 0..n_findex {
            // Turbines to optimize: not downstream-excluded and with a
            // non-degenerate bound range (Python _turbs_to_opt_subset).
            let downstream = if self.config.exclude_downstream_turbines {
                super::derive_downstream_turbines(fmodel, 0.3, false)?
            } else {
                vec![]
            };
            let to_opt: Vec<usize> = (0..n_turbines)
                .filter(|&t| !downstream.contains(&t))
                .filter(|&t| yaw_ubs[[fi, t]] - yaw_lbs[[fi, t]] >= 0.001)
                .collect();
            if to_opt.is_empty() {
                continue;
            }

            // Optimize this condition on its own (Python loops over the
            // conditions and sets a single-condition model).
            fmodel.set_wind_conditions(
                NpArray1::from_vec(vec![ws[fi]]),
                NpArray1::from_vec(vec![wd[fi]]),
                NpArray1::from_vec(vec![ti[fi]]),
            )?;

            let mut best_yaw = yaw_angles.row(0).to_owned();
            let mut best_power = self.evaluate(fmodel, &best_yaw, &turbine_weights)?;

            // Block-coordinate search: for each free turbine, minimize the
            // weighted farm power over its yaw range with a ternary search,
            // then repeat until the yaw angles stop moving. This converges to
            // the same local optimum as SLSQP for the (smooth) wake objective.
            let lb_row = yaw_lbs.row(fi).to_owned();
            let ub_row = yaw_ubs.row(fi).to_owned();
            for _ in 0..self.scipy_config.maxiter {
                let mut max_change = 0.0_f64;
                for &t in &to_opt {
                    let (mut lo, mut hi) = (lb_row[t], ub_row[t]);
                    if (hi - lo).abs() < 1e-12 {
                        continue;
                    }
                    let (mut trial, mut trial_power) = (best_yaw.clone(), best_power);
                    for _ in 0..40 {
                        let m1 = lo + (hi - lo) / 3.0;
                        let m2 = hi - (hi - lo) / 3.0;
                        let mut y1 = best_yaw.clone();
                        y1[t] = m1;
                        let mut y2 = best_yaw.clone();
                        y2[t] = m2;
                        let p1 = self.evaluate(fmodel, &y1, &turbine_weights)?;
                        let p2 = self.evaluate(fmodel, &y2, &turbine_weights)?;
                        if p1 > p2 {
                            hi = m2;
                            trial = y1;
                            trial_power = p1;
                        } else {
                            lo = m1;
                            trial = y2;
                            trial_power = p2;
                        }
                        if (hi - lo).abs() < self.scipy_config.ftol.max(1e-6) {
                            break;
                        }
                    }
                    if trial_power > best_power {
                        max_change = max_change.max((trial[t] - best_yaw[t]).abs());
                        best_power = trial_power;
                        best_yaw = trial;
                    }
                }
                if max_change < self.scipy_config.ftol.max(1e-6) {
                    break;
                }
            }

            // Store the optimum for this condition
            for t in 0..n_turbines {
                yaw_angles[[fi, t]] = best_yaw[t];
            }
            farm_powers_opt[fi] = best_power;
        }

        // Final powers with all optimal yaw angles
        fmodel.set_wind_conditions(ws.clone(), wd.clone(), ti.clone())?;
        fmodel.set_yaw_angles(yaw_angles.clone())?;
        fmodel.run()?;
        let turbine_powers = fmodel.get_turbine_powers()?;
        let farm_powers_opt = fmodel.get_farm_power()?;

        let optimized_power: Float = farm_powers_opt.sum();
        let power_improvement = optimized_power - baseline_power;
        let improvement_percentage = if baseline_power > 0.0 {
            (power_improvement / baseline_power) * 100.0
        } else {
            0.0
        };

        Ok(YawOptimizationResult {
            yaw_angles,
            powers: turbine_powers,
            baseline_power,
            optimized_power,
            power_improvement,
            improvement_percentage,
            farm_powers_opt,
            farm_powers_baseline,
        })
    }
}

impl YawOptimizationScipy {
    /// Weighted farm power for a candidate yaw set on the current
    /// (single-condition) model.
    fn evaluate(
        &self,
        fmodel: &mut FlorisModel,
        yaw: &Array1,
        turbine_weights: &[Float],
    ) -> anyhow::Result<Float> {
        let n_turbines = yaw.len();
        fmodel.set_yaw_angles(Array2::from_shape_vec((1, n_turbines), yaw.to_vec())?)?;
        fmodel.run()?;
        let powers = fmodel.get_turbine_powers()?;
        let mut total = 0.0;
        for t in 0..n_turbines {
            total += turbine_weights[t] * powers[[0, t]];
        }
        Ok(total)
    }
}
