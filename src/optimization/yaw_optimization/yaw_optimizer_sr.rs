//! Serial-Refine (SR) yaw optimizer — port of Python
//! `floris.optimization.yaw_optimization.yaw_optimizer_sr.YawOptimizationSR`.
//!
//! Algorithm (Python yaw_optimizer_sr.py):
//! 1. Baseline weighted farm power per condition at the baseline yaw angles
//!    (`calc_baseline_power` is always true for SR).
//! 2. For each pass (`Ny_passes`, default `[5, 4]`), for each turbine from
//!    front to back (per wind direction, by rotated x):
//!    - build an evaluation grid of `Ny` yaw values for that turbine
//!      (pass 0: `linspace(lb, ub, Ny)`; later passes:
//!      `linspace(lb, ub, Ny + 1)` minus the middle point, which was
//!      already evaluated in the previous pass);
//!    - evaluate the weighted farm power for every grid row and condition in
//!      one batched FLORIS run (`_calc_powers_with_memory`, without the
//!      memory shortcut — same results);
//!    - take the per-condition argmax (NaN-ignoring) and update the optimum
//!      where it improves;
//!    - shrink that turbine's per-condition bounds to
//!      `[opt - dx/2, opt + dx/2]` (clipped to the global yaw limits), where
//!      `dx` is the grid spacing.
//! 3. Final per-turbine powers at the optimal yaw angles.

use crate::floris_model::FlorisModel;
use crate::types::{Array1, Array2, Float};

use super::{YawOptimization, YawOptimizationConfig, YawOptimizationResult};

/// Serial-Refine configuration (Python `Ny_passes=[5, 4]`).
///
/// The first entry may be odd or even; further entries must be even so the
/// refined grid never re-evaluates the previous pass's points.
#[derive(Debug, Clone)]
pub struct SerialRefineConfig {
    pub ny_passes: Vec<usize>,
}

impl Default for SerialRefineConfig {
    fn default() -> Self {
        Self {
            ny_passes: vec![5, 4],
        }
    }
}

impl SerialRefineConfig {
    /// Python yaw_optimizer_sr.py:83-93 validation.
    pub fn validate(&self) -> anyhow::Result<()> {
        for (nii, &ny) in self.ny_passes.iter().enumerate() {
            if ny < 2 {
                anyhow::bail!("Each entry in Ny_passes must have a value of at least 2.");
            }
            if nii > 0 && (ny + 1) % 2 == 0 {
                anyhow::bail!(
                    "The second and further entries of Ny_passes must be even numbers. \
                     This is to ensure the same yaw angles are not evaluated twice \
                     between passes."
                );
            }
        }
        Ok(())
    }
}

/// Serial-Refine yaw optimizer
///
/// Optimizes turbines one at a time from front to back with multiple
/// refinement passes, exactly mirroring Python `YawOptimizationSR`.
#[derive(Debug, Clone)]
pub struct YawOptimizationSR {
    pub config: YawOptimizationConfig,
    pub sr_config: SerialRefineConfig,
    pub print_progress: bool,
}

impl YawOptimizationSR {
    pub fn new() -> Self {
        Self {
            config: YawOptimizationConfig::default(),
            sr_config: SerialRefineConfig::default(),
            print_progress: true,
        }
    }

    /// Python `_get_turbine_orders` (yaw_optimizer_sr.py:101-112): for each
    /// wind direction, turbines sorted by x in the wind-aligned frame
    /// (`np.argsort` of the rotated layout x; numpy's small-array sort is
    /// stable for ties, matching Rust's stable sort).
    fn turbine_orders(fmodel: &FlorisModel, wind_directions: &Array1) -> Vec<Vec<usize>> {
        let n_turbines = fmodel.farm().n_turbines();
        let layout_x = fmodel.farm().layout_x.clone();
        let layout_y = fmodel.farm().layout_y.clone();

        wind_directions
            .iter()
            .map(|&wd| {
                let angle = (wd - 270.0) * std::f64::consts::PI / 180.0;
                let (c, s) = (angle.cos(), angle.sin());
                let mut order: Vec<(usize, Float)> = (0..n_turbines)
                    .map(|i| (i, c * layout_x[i] - s * layout_y[i]))
                    .collect();
                // Stable sort: ties keep index order (numpy argsort on small
                // arrays uses a stable insertion sort).
                order.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
                order.into_iter().map(|(i, _)| i).collect()
            })
            .collect()
    }

    /// Python `_calculate_farm_power` (yaw_optimization_base.py:306-375):
    /// one batched FLORIS run for the given conditions and yaw angles, then
    /// the turbine-weighted farm power per condition.
    ///
    /// `set_wind_conditions` resets the control arrays (yaw, setpoints, AWC)
    /// to their reference defaults; the optimizer must preserve the
    /// power setpoints that were set by `set_operation(disable_turbines=...)`
    /// across evaluations (Python deep-copies the model, keeping them).
    fn calculate_farm_powers(
        fmodel: &mut FlorisModel,
        wd: &Array1,
        ws: &Array1,
        ti: &Array1,
        yaw_angles: &Array2,
        turbine_weights: &[Float],
    ) -> anyhow::Result<Array1> {
        let n = wd.len();

        // Preserve power setpoints across the condition reset. The optimizer
        // stacks Ny copies of the conditions per evaluation, so the saved
        // (n_findex, n_turbines) table must be tiled to the batch size
        // (Python's deep copy + set() handles this implicitly).
        let saved_setpoints = fmodel.core().farm.power_setpoints.clone();
        let n_saved = saved_setpoints.shape()[0];
        let n_turbines = saved_setpoints.shape()[1];

        fmodel.set_wind_conditions(ws.clone(), wd.clone(), ti.clone())?;

        // Tile the saved setpoints if the batch has more conditions
        let n_batch = ws.len();
        let setpoints_to_restore = if n_batch == n_saved {
            saved_setpoints
        } else {
            // Tile row i -> rows i, i+n_saved, i+2*n_saved, ...
            // This matches Python's np.tile(setpoints, (Ny, 1))
            Array2::from_shape_fn((n_batch, n_turbines), |(r, c)| {
                saved_setpoints[[r % n_saved, c]]
            })
        };
        fmodel
            .core_mut()
            .farm
            .set_power_setpoints(setpoints_to_restore);
        fmodel.set_yaw_angles(yaw_angles.clone())?;
        fmodel.run()?;
        let powers = fmodel.get_turbine_powers()?;

        let mut farm_power = Array1::zeros(n);
        for fi in 0..n {
            let mut sum = 0.0;
            for (ti_idx, &w) in turbine_weights.iter().enumerate() {
                sum += w * powers[[fi, ti_idx]];
            }
            farm_power[fi] = sum;
        }
        Ok(farm_power)
    }

    fn linspace(a: Float, b: Float, n: usize) -> Vec<Float> {
        if n == 1 {
            return vec![a];
        }
        let step = (b - a) / (n as Float - 1.0);
        (0..n).map(|i| a + step * i as Float).collect()
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
        self.sr_config.validate()?;

        if !fmodel.state().initialized {
            anyhow::bail!("FlorisModel must be run before optimization");
        }

        let n_findex = fmodel.flow_field().n_findex;
        let n_turbines = fmodel.farm().n_turbines();

        // Ambient conditions of the (subset) model
        let wd = fmodel.flow_field().wind_directions.clone();
        let ws = fmodel.flow_field().wind_speeds.clone();
        let ti = fmodel.flow_field().turbulence_intensities.clone();
        let turbine_weights = self
            .config
            .turbine_weights
            .clone()
            .unwrap_or_else(|| vec![1.0; n_turbines]);

        // Baseline yaw angles (Python defaults to 0.0)
        let yaw_angles_baseline = self
            .config
            .yaw_angles_baseline
            .clone()
            .unwrap_or_else(|| Array2::zeros((n_findex, n_turbines)));

        // Front-to-back turbine order per wind direction
        let orders = Self::turbine_orders(fmodel, &wd);

        // Baseline weighted farm powers (Python _calculate_baseline_farm_power)
        let farm_powers_baseline = Self::calculate_farm_powers(
            fmodel,
            &wd,
            &ws,
            &ti,
            &yaw_angles_baseline,
            &turbine_weights,
        )?;

        // Optimal solutions start at the baseline
        let mut yaw_opt = yaw_angles_baseline.clone();
        let mut farm_power_opt = farm_powers_baseline.clone();

        // Per-condition yaw bounds (Python __init__:180-181). The config may
        // provide per-(findex, turbine) arrays (Python's `minimum_yaw_angle`
        // accepts scalars or arrays), which override the scalar bounds.
        // `yaw_limit_min`/`yaw_limit_max` keep the ORIGINAL limits (used for
        // clipping, Python yaw_optimizer_sr.py:208-221 and 311-320), while
        // `yaw_lbs`/`yaw_ubs` hold the shrinking search bounds.
        let mut yaw_limit_min =
            Array2::from_elem((n_findex, n_turbines), self.config.minimum_yaw_angle);
        let mut yaw_limit_max =
            Array2::from_elem((n_findex, n_turbines), self.config.maximum_yaw_angle);
        if let Some(ref arr) = self.config.minimum_yaw_angle_array {
            if arr.shape() == [n_findex, n_turbines] {
                for fi in 0..n_findex {
                    for ti in 0..n_turbines {
                        yaw_limit_min[[fi, ti]] = arr[[fi, ti]];
                    }
                }
            }
        }
        if let Some(ref arr) = self.config.maximum_yaw_angle_array {
            if arr.shape() == [n_findex, n_turbines] {
                for fi in 0..n_findex {
                    for ti in 0..n_turbines {
                        yaw_limit_max[[fi, ti]] = arr[[fi, ti]];
                    }
                }
            }
        }
        let mut yaw_lbs = yaw_limit_min.clone();
        let mut yaw_ubs = yaw_limit_max.clone();

        let total_steps = self.sr_config.ny_passes.len() * n_turbines;
        let mut ii = 0usize;
        // The depth index selects a different turbine per findex (front-to-back
        // order depends on the wind direction), so it cannot be replaced by
        // iterating the turbine ids directly.
        let depths: Vec<usize> = (0..n_turbines).collect();

        for (pass_depth, &ny) in self.sr_config.ny_passes.iter().enumerate() {
            for &turbine_depth in &depths {
                let p = 100.0 * ii as Float / total_steps as Float;
                ii += 1;
                if self.print_progress {
                    println!(
                        "[Serial Refine] Processing pass={pass_depth}, turbine_depth={turbine_depth} ({p:.1}%)"
                    );
                }

                // Python _generate_evaluation_grid: (Ny, n_findex, n_turbines),
                // starting from Ny copies of the current optimum, with the
                // manipulated turbine's column replaced per condition.
                let mut grid = vec![yaw_opt.clone(); ny];
                for fi in 0..n_findex {
                    let turbid = orders[fi][turbine_depth];

                    // Saturate to the allowable yaw limits (Python
                    // yaw_optimizer_sr.py:208-221): clip the current search
                    // bounds to the ORIGINAL per-turbine limit interval.
                    let lb = yaw_lbs[[fi, turbid]]
                        .clamp(yaw_limit_min[[fi, turbid]], yaw_limit_max[[fi, turbid]]);
                    let ub = yaw_ubs[[fi, turbid]]
                        .clamp(yaw_limit_min[[fi, turbid]], yaw_limit_max[[fi, turbid]]);

                    let values: Vec<Float> = if pass_depth == 0 {
                        // Evaluate all possible coordinates
                        Self::linspace(lb, ub, ny)
                    } else {
                        // Remove the middle point: evaluated in the previous pass
                        let c = ny / 2;
                        Self::linspace(lb, ub, ny + 1)
                            .iter()
                            .enumerate()
                            .filter(|(k, _)| *k != c)
                            .map(|(_, &v)| v)
                            .collect()
                    };

                    for (k, &v) in values.iter().enumerate() {
                        grid[k][[fi, turbid]] = v;
                    }
                }

                // Evaluate the grid in one batched run: conditions tiled Ny
                // times (Python _calc_powers_with_memory), yaw stacked row-major.
                let stacked_wd = Array1::from_iter((0..ny).flat_map(|_| wd.iter().copied()));
                let stacked_ws = Array1::from_iter((0..ny).flat_map(|_| ws.iter().copied()));
                let stacked_ti = Array1::from_iter((0..ny).flat_map(|_| ti.iter().copied()));
                let stacked_yaw = Array2::from_shape_fn((ny * n_findex, n_turbines), |(row, t)| {
                    grid[row / n_findex][[row % n_findex, t]]
                });

                let farm_powers = Self::calculate_farm_powers(
                    fmodel,
                    &stacked_wd,
                    &stacked_ws,
                    &stacked_ti,
                    &stacked_yaw,
                    &turbine_weights,
                )?;

                // Per-condition argmax (np.nanargmax) and monotone update
                for fi in 0..n_findex {
                    let mut best_k = 0usize;
                    let mut best_p = Float::NEG_INFINITY;
                    for k in 0..ny {
                        let v = farm_powers[k * n_findex + fi];
                        if v.is_finite() && v > best_p {
                            best_p = v;
                            best_k = k;
                        }
                    }

                    if best_p > farm_power_opt[fi] {
                        let turbid = orders[fi][turbine_depth];
                        let new_yaw = grid[best_k][[fi, turbid]];
                        let dx = grid[1][[fi, turbid]] - grid[0][[fi, turbid]];

                        farm_power_opt[fi] = best_p;
                        yaw_opt[[fi, turbid]] = new_yaw;

                        // Shrink the bounds around the new optimum, clipped
                        // to the original limits (Python
                        // yaw_optimizer_sr.py:311-320).
                        yaw_lbs[[fi, turbid]] = (new_yaw - 0.5 * dx)
                            .clamp(yaw_limit_min[[fi, turbid]], yaw_limit_max[[fi, turbid]]);
                        yaw_ubs[[fi, turbid]] = (new_yaw + 0.5 * dx)
                            .clamp(yaw_limit_min[[fi, turbid]], yaw_limit_max[[fi, turbid]]);
                    }
                }
            }
        }

        // Final per-turbine powers at the optimal yaw angles
        let powers = {
            fmodel.set_wind_conditions(ws.clone(), wd.clone(), ti.clone())?;
            fmodel.set_yaw_angles(yaw_opt.clone())?;
            fmodel.run()?;
            fmodel.get_turbine_powers()?
        };

        let baseline_power: Float = farm_powers_baseline.sum();
        let optimized_power: Float = farm_power_opt.sum();
        let power_improvement = optimized_power - baseline_power;
        let improvement_percentage = if baseline_power > 0.0 {
            power_improvement / baseline_power * 100.0
        } else {
            0.0
        };

        Ok(YawOptimizationResult {
            yaw_angles: yaw_opt,
            powers,
            baseline_power,
            optimized_power,
            power_improvement,
            improvement_percentage,
            farm_powers_opt: farm_power_opt,
            farm_powers_baseline,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sr_config_validation() {
        assert!(SerialRefineConfig::default().validate().is_ok());
        assert!(SerialRefineConfig {
            ny_passes: vec![5, 4, 6]
        }
        .validate()
        .is_ok());
        assert!(SerialRefineConfig {
            ny_passes: vec![5, 5]
        }
        .validate()
        .is_err());
        assert!(SerialRefineConfig {
            ny_passes: vec![1, 4]
        }
        .validate()
        .is_err());
    }

    #[test]
    fn test_linspace() {
        let v = YawOptimizationSR::linspace(0.0, 25.0, 5);
        assert_eq!(v.len(), 5);
        assert!((v[0] - 0.0).abs() < 1e-12);
        assert!((v[4] - 25.0).abs() < 1e-12);
        assert!((v[2] - 12.5).abs() < 1e-12);
    }
}
