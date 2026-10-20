//! Gridded layout optimizer — port of Python
//! `floris.optimization.layout_optimization.layout_optimization_gridded.LayoutOptimizationGridded`.
//!
//! Places as many turbines as possible inside a boundary by sweeping a
//! (square or hexagonal) grid over the boundary with rotations and
//! translations, and keeping the layout that fits the most turbines.

use crate::optimization::layout_optimization::layout_optimization_base::{
    is_point_in_boundary, Boundary, LayoutOptimizationResult,
};
use crate::types::{Array1, Float};

/// Gridded layout optimizer.
#[derive(Debug, Clone)]
pub struct LayoutOptimizationGridded {
    /// Boundary the turbines must lie inside.
    pub boundaries: Boundary,
    /// Minimum spacing between turbines (meters).
    pub min_dist: Float,
    /// Grid rotation step (degrees).
    pub rotation_step: Float,
    /// Grid rotation range (degrees).
    pub rotation_range: (Float, Float),
    /// Grid translation step (meters).
    pub translation_step: Float,
    /// Grid translation range (meters).
    pub translation_range: (Float, Float),
    /// Use hexagonal instead of square packing.
    pub hexagonal_packing: bool,
    /// Rotor diameter of the first turbine (meters).
    pub rotor_diameter: Float,
    /// Number of turbines in the best layout.
    pub n_turbines_max: usize,
    /// Optimal layout x coordinates.
    pub x_opt: Array1,
    /// Optimal layout y coordinates.
    pub y_opt: Array1,
}

impl LayoutOptimizationGridded {
    /// Create a gridded layout optimizer.
    ///
    /// * `min_dist_d` — minimum spacing in rotor diameters (negative means 5D).
    /// * `translation_step_d` — translation step in rotor diameters
    ///   (negative means 1D).
    pub fn new(
        fmodel: &crate::floris_model::FlorisModel,
        boundaries: Boundary,
        min_dist_d: Option<Float>,
        translation_step_d: Option<Float>,
    ) -> crate::Result<Self> {
        // Rotor diameter of the first turbine (Python
        // `fmodel.core.farm.rotor_diameters.flat[0]`). The farm's runtime
        // array can be uninitialized before the first run, so take the value
        // from the turbine definition.
        let rotor_diameter = fmodel
            .farm()
            .turbines
            .first()
            .map(|t| t.turbine_type.rotor_diameter)
            .unwrap_or(125.88);

        let min_dist = match min_dist_d {
            Some(d) if d >= 0.0 => d * rotor_diameter,
            _ => 5.0 * rotor_diameter,
        };
        let translation_step = match translation_step_d {
            Some(d) if d >= 0.0 => d * rotor_diameter,
            _ => 1.0 * rotor_diameter,
        };

        Ok(Self {
            boundaries,
            min_dist,
            rotation_step: 5.0,
            rotation_range: (0.0, 360.0),
            translation_step,
            translation_range: (0.0, min_dist),
            hexagonal_packing: false,
            rotor_diameter,
            n_turbines_max: 0,
            x_opt: Array1::zeros(0),
            y_opt: Array1::zeros(0),
        })
    }

    /// Build the (trimmed) grid of candidate positions.
    ///
    /// Python `LayoutOptimizationGridded.__init__` (grid.py:104-151).
    fn build_grid(&self) -> (Vec<Float>, Vec<Float>) {
        let (xmin, xmax, ymin, ymax) = self.boundaries.bounds();
        let d = 1.1 * ((xmax - xmin).powi(2) + (ymax - ymin).powi(2)).sqrt();
        let min_dist = self.min_dist;

        // np.arange(0, d + min_dist, min_dist)
        let n_grid = ((d + min_dist) / min_dist).floor() as usize + 1;
        let grid_1d: Vec<Float> = (0..n_grid).map(|i| i as Float * min_dist).collect();

        let mut x_locs: Vec<Float> = Vec::new();
        let mut y_locs: Vec<Float> = Vec::new();
        if self.hexagonal_packing {
            // tile(grid_1D, (len,1)) with the odd rows offset by 0.5*min_dist
            for (i, &g) in grid_1d.iter().enumerate() {
                let offset = if i % 2 == 1 { 0.5 * min_dist } else { 0.0 };
                for &g2 in grid_1d.iter() {
                    x_locs.push(g2 + offset);
                    y_locs.push(3.0_f64.sqrt() / 2.0 * g);
                }
            }
        } else {
            // np.meshgrid(grid_1D, grid_1D) flattened
            for &gy in grid_1d.iter() {
                for &gx in grid_1d.iter() {
                    x_locs.push(gx);
                    y_locs.push(gy);
                }
            }
        }

        // Center the grid on the boundary
        let mean_x = x_locs.iter().sum::<Float>() / x_locs.len() as Float;
        let mean_y = y_locs.iter().sum::<Float>() / y_locs.len() as Float;
        let cx = 0.5 * (xmax + xmin);
        let cy = 0.5 * (ymax + ymin);
        for v in x_locs.iter_mut() {
            *v = *v - mean_x + cx;
        }
        for v in y_locs.iter_mut() {
            *v = *v - mean_y + cy;
        }

        // Trim to a circle to avoid wasted computation
        // (radius = (grid.max() - grid.min() + min_dist) / 2)
        let gmax = grid_1d.last().copied().unwrap_or(0.0);
        let gmin = grid_1d.first().copied().unwrap_or(0.0);
        let radius = (gmax - gmin + min_dist) / 2.0;
        let mut x_max = f64::NEG_INFINITY;
        let mut x_min = f64::INFINITY;
        for &v in x_locs.iter() {
            if v > x_max {
                x_max = v;
            }
            if v < x_min {
                x_min = v;
            }
        }
        let mut y_max = f64::NEG_INFINITY;
        let mut y_min = f64::INFINITY;
        for &v in y_locs.iter() {
            if v > y_max {
                y_max = v;
            }
            if v < y_min {
                y_min = v;
            }
        }
        let center_x = 0.5 * (x_max + x_min);
        let center_y = 0.5 * (y_max + y_min);

        let mut kept_x = Vec::new();
        let mut kept_y = Vec::new();
        for (x, y) in x_locs.iter().zip(y_locs.iter()) {
            let dx = x - center_x;
            let dy = y - center_y;
            if (dx * dx + dy * dy).sqrt() <= radius {
                kept_x.push(*x);
                kept_y.push(*y);
            }
        }
        (kept_x, kept_y)
    }

    /// Sweep rotations and translations, keeping the layout that fits the
    /// most turbines (Python `optimize`, grid.py:153-200).
    pub fn optimize(&mut self) -> LayoutOptimizationResult {
        let (x_grid, y_grid) = self.build_grid();

        // Rotation range: square grid has 90 degree symmetry, hexagonal 60
        let symmetry = if self.hexagonal_packing { 60.0 } else { 90.0 };
        let rot_end = self.rotation_range.1.min(self.rotation_range.0 + symmetry);
        let n_rots =
            (((rot_end - self.rotation_range.0) / self.rotation_step).ceil() as usize).max(1);
        let rotations: Vec<Float> = (0..n_rots)
            .map(|i| self.rotation_range.0 + i as Float * self.rotation_step)
            .collect();

        let n_trans = (((self.translation_range.1 - self.translation_range.0)
            / self.translation_step)
            .ceil() as usize)
            .max(1);
        let mut translations: Vec<Float> = (0..n_trans)
            .map(|i| self.translation_range.0 + i as Float * self.translation_step)
            .collect();
        // np.concatenate([-np.flip(translations), translations])
        let flipped: Vec<Float> = translations.iter().rev().map(|v| -v).collect();
        translations.splice(0..0, flipped);

        // Candidate layouts. Python builds them with
        //   translations_all = repeat(translation_matrices, n_rots, axis=0)
        //   rotations_all   = tile(rotation_matrices, (n_trans**2, 1, 1))
        // so the enumeration order is: translation pair (ty slowest, tx
        // fastest), then the rotations innermost. The first maximizing
        // candidate wins (np.argmax), so the order matters.
        let mut best_count = 0usize;
        let mut best_x: Vec<Float> = Vec::new();
        let mut best_y: Vec<Float> = Vec::new();

        for &ty in &translations {
            for &tx in &translations {
                for &rot_deg in &rotations {
                    let rad = rot_deg.to_radians();
                    let (c, s) = (rad.cos(), rad.sin());
                    let mut candidate: Vec<(Float, Float)> = Vec::with_capacity(x_grid.len());
                    for (x, y) in x_grid.iter().zip(y_grid.iter()) {
                        candidate.push((c * x - s * y + tx, s * x + c * y + ty));
                    }
                    let count = candidate
                        .iter()
                        .filter(|(x, y)| is_point_in_boundary(*x, *y, &self.boundaries))
                        .count();
                    if count > best_count {
                        best_count = count;
                        best_x = candidate.iter().map(|(x, _)| *x).collect();
                        best_y = candidate.iter().map(|(_, y)| *y).collect();
                    }
                }
            }
        }

        // Keep only the in-bounds turbines
        let mut x_opt = Vec::new();
        let mut y_opt = Vec::new();
        for (x, y) in best_x.iter().zip(best_y.iter()) {
            if is_point_in_boundary(*x, *y, &self.boundaries) {
                x_opt.push(*x);
                y_opt.push(*y);
            }
        }

        self.n_turbines_max = best_count;
        self.x_opt = Array1::from_vec(x_opt);
        self.y_opt = Array1::from_vec(y_opt);

        LayoutOptimizationResult {
            x: self.x_opt.clone(),
            y: self.y_opt.clone(),
            value: best_count as Float,
            iterations: 0,
            improvement_pct: 0.0,
        }
    }
}
