use crate::types::Float;
use crate::floris_model::FlorisModel;
use crate::Array2;
use crate::Array4;
use crate::optimization::YawAngleBounds;

/// Calculate cosine loss factor for yaw misalignment
pub fn yaw_cosine_loss(yaw_angle: Float, exponent: Float) -> Float {
    let yaw_rad = yaw_angle.to_radians();
    let cos_yaw = yaw_rad.cos();
    (cos_yaw.powf(exponent)).max(0.0)
}

/// Calculate derivative of cosine loss with respect to yaw angle
pub fn yaw_cosine_loss_derivative(yaw_angle: Float, exponent: Float) -> Float {
    let yaw_rad = yaw_angle.to_radians();
    let cos_yaw = yaw_rad.cos();
    let sin_yaw = yaw_rad.sin();
    
    if cos_yaw.abs() < 1e-10 {
        return 0.0;
    }
    
    -exponent * cos_yaw.powf(exponent - 1.0) * sin_yaw
}

/// Estimate wake deflection from yaw angle
pub fn estimate_wake_deflection_angle(
    yaw_angle: Float,
    thrust_coefficient: Float,
    rotor_diameter: Float,
    downstream_distance: Float,
    kd: Float,
    ad: Float,
) -> Float {
    if thrust_coefficient <= 0.0 || yaw_angle == 0.0 {
        return 0.0;
    }

    let yaw_rad = yaw_angle.to_radians();
    let axial_induction = if thrust_coefficient < 0.96 {
        0.5 * (1.0 - (1.0 - thrust_coefficient).sqrt())
    } else {
        0.143 + (0.0203 - 0.6427 * (0.889 - thrust_coefficient)).sqrt().max(0.0)
    };

    let c = ad / kd;
    let exp_term = (-kd * downstream_distance / rotor_diameter).exp();
    c * axial_induction * (1.0 - exp_term) * yaw_rad * downstream_distance
}

/// Check if a turbine is in the wake of another turbine
pub fn is_turbine_in_wake(
    upstream_idx: usize,
    downstream_idx: usize,
    x_sorted: &Array4,
    y_sorted: &Array4,
    _rotor_diameter: Float,
    wake_slope: Float,
) -> bool {
    let dx = x_sorted[[0, 0, 0, downstream_idx]] - x_sorted[[0, 0, 0, upstream_idx]];
    let dy = y_sorted[[0, 0, 0, downstream_idx]] - y_sorted[[0, 0, 0, upstream_idx]];
    let distance = (dx * dx + dy * dy).sqrt();

    if dx <= 0.0 {
        return false;
    }

    let wake_width = wake_slope * distance;
    dy.abs() < wake_width / 2.0
}

/// Derive downstream turbine indices for exclusion
pub fn derive_downstream_turbines(
    fmodel: &FlorisModel,
    wake_slope: Float,
    _sort_turbines: bool,
) -> anyhow::Result<Vec<usize>> {
    if !fmodel.state.initialized {
        anyhow::bail!("FlorisModel must be run before deriving downstream turbines");
    }

    let grid = fmodel.grid.as_ref()
        .ok_or_else(|| anyhow::anyhow!("Grid not initialized"))?;

    let n_turbines = fmodel.farm.n_turbines();
    let x_sorted = grid.x_sorted();
    let y_sorted = grid.y_sorted();

    let mut downstream_indices = Vec::new();

    for ui in 0..n_turbines {
        for ti in 0..n_turbines {
            if ui != ti {
                if is_turbine_in_wake(ui, ti, x_sorted, y_sorted, fmodel.farm.rotor_diameters[0], wake_slope) {
                    if !downstream_indices.contains(&ti) {
                        downstream_indices.push(ti);
                    }
                }
            }
        }
    }

    Ok(downstream_indices)
}

/// Calculate turbine weights for optimization
pub fn calculate_turbine_weights(n_turbines: usize, weight: Option<Float>) -> Vec<Float> {
    match weight {
        Some(w) => vec![w; n_turbines],
        None => vec![1.0; n_turbines],
    }
}

/// Optimize yaw for a single findex using gradient-based search
pub fn optimize_yaw_single_findex(
    fmodel: &mut FlorisModel,
    yaw_angles: &mut Array2,
    findex: usize,
    turbine_weights: &[Float],
    yaw_bounds: (Float, Float),
) -> anyhow::Result<()> {
    let n_turbines = fmodel.farm.n_turbines();
    let (min_yaw, max_yaw) = yaw_bounds;

    for ti in 0..n_turbines {
        let current_yaw = yaw_angles[[findex, ti]];
        let weight = turbine_weights[ti];

        let delta_yaw = 1.0;
        let yaw_plus = (current_yaw + delta_yaw).clamp(min_yaw, max_yaw);
        let yaw_minus = (current_yaw - delta_yaw).clamp(min_yaw, max_yaw);

        yaw_angles[[findex, ti]] = yaw_plus;
        fmodel.set_yaw_angles(yaw_angles.clone())?;
        fmodel.run()?;
        let power_plus = fmodel.get_farm_power()[findex];

        yaw_angles[[findex, ti]] = yaw_minus;
        fmodel.set_yaw_angles(yaw_angles.clone())?;
        fmodel.run()?;
        let power_minus = fmodel.get_farm_power()[findex];

        let gradient = (power_plus - power_minus) / (2.0 * delta_yaw);
        let learning_rate = 2.0 * weight;
        let step = gradient * learning_rate;

        let new_yaw = (current_yaw - step).clamp(min_yaw, max_yaw);
        yaw_angles[[findex, ti]] = new_yaw;
    }

    Ok(())
}

/// Golden section search for yaw angle optimization
pub fn golden_section_search_yaw<F>(
    f: F,
    mut a: Float,
    mut b: Float,
    tol: Float,
    max_iter: usize,
) -> (Float, Float)
where
    F: Fn(Float) -> Float,
{
    let golden_ratio = 1.618033988749895;
    let inv_golden_ratio = 1.0 / golden_ratio;

    let mut c = b - inv_golden_ratio * (b - a);
    let mut d = a + inv_golden_ratio * (b - a);

    let mut fc = f(c);
    let mut fd = f(d);

    for _ in 0..max_iter {
        if (b - a) < tol {
            break;
        }

        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - inv_golden_ratio * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + inv_golden_ratio * (b - a);
            fd = f(d);
        }
    }

    let midpoint = (a + b) / 2.0;
    let max_value = f(midpoint);

    (midpoint, max_value)
}

/// Coordinate descent yaw optimization
pub fn coordinate_descent_yaw<F>(
    yaw_angles: &mut Array2,
    get_power_fn: F,
    bounds: &YawAngleBounds,
    max_iter: usize,
    tolerance: Float,
) -> Float
where
    F: Fn(&Array2) -> Float,
{
    let mut prev_power = get_power_fn(yaw_angles);

    for _ in 0..max_iter {
        let n_turbines = yaw_angles.shape()[1];
        let mut improved = false;

        for ti in 0..n_turbines {
            let current_yaw = yaw_angles[[0, ti]];
            let perturbations = [-15.0, -10.0, -5.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0, 10.0, 15.0];
            let mut best_yaw = current_yaw;
            let mut best_power = prev_power;

            for &delta in &perturbations {
                let new_yaw = (current_yaw + delta).clamp(bounds.min_yaw, bounds.max_yaw);
                yaw_angles[[0, ti]] = new_yaw;
                let power = get_power_fn(yaw_angles);

                if power > best_power {
                    best_power = power;
                    best_yaw = new_yaw;
                    improved = true;
                }
            }

            yaw_angles[[0, ti]] = best_yaw;
        }

        let new_power = get_power_fn(yaw_angles);

        if (new_power - prev_power).abs() < tolerance * prev_power.abs() {
            break;
        }

        if !improved {
            break;
        }

        prev_power = new_power;
    }

    prev_power
}

/// Calculate derivative of power with respect to yaw angle
pub fn yaw_angle_derivative(power_plus: Float, power_minus: Float, dx: Float) -> Float {
    if dx == 0.0 {
        return 0.0;
    }
    (power_plus - power_minus) / (2.0 * dx)
}
