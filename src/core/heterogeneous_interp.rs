//! Spatial interpolation of heterogeneous-inflow speed multipliers onto the
//! flow-field grid points.
//!
//! Python applies the multipliers with `scipy.interpolate.LinearNDInterpolator`
//! (or `NearestNDInterpolator`) per wind condition and multiplies the initial
//! velocity field by the interpolated value (floris/core/flow_field.py:
//! `initialize_velocity_field`). Points outside the convex hull of the
//! user-defined locations are filled with the freestream value 1.0.
//!
//! This module provides the same behaviour:
//! - `nearest`: nearest user-defined location (matches scipy exactly);
//! - `linear`: barycentric interpolation inside a Delaunay triangulation of
//!   the user-defined locations (the same construction scipy uses; for
//!   co-circular points the chosen diagonal may differ from qhull's, which
//!   can shift the interpolated value slightly).

use crate::core::models::InterpMethod;
use crate::types::{Array1, Float};

/// Interpolate the speed multipliers of one condition at arbitrary points.
///
/// * `points_x`, `points_y` - the user-defined multiplier locations
/// * `values` - the multiplier at each user-defined location
/// * `target_x`, `target_y` - where to evaluate
///
/// Returns `None` for targets outside the convex hull (the caller fills
/// them with the freestream value 1.0, as Python does).
pub fn interpolate_multipliers(
    interp_method: InterpMethod,
    points_x: &Array1,
    points_y: &Array1,
    values: &Array1,
    target_x: &[Float],
    target_y: &[Float],
) -> Vec<Option<Float>> {
    let n_points = points_x.len();
    assert_eq!(n_points, points_y.len());
    assert_eq!(n_points, values.len());

    let mut out = vec![None; target_x.len()];

    if n_points == 0 {
        return out;
    }

    // Nearest neighbour (scipy NearestNDInterpolator)
    if matches!(interp_method, InterpMethod::Nearest) {
        for (t, (tx, ty)) in target_x.iter().zip(target_y.iter()).enumerate() {
            let mut best = 0usize;
            let mut best_d = Float::INFINITY;
            for p in 0..n_points {
                let dx = points_x[p] - tx;
                let dy = points_y[p] - ty;
                let d = dx * dx + dy * dy;
                if d < best_d {
                    best_d = d;
                    best = p;
                }
            }
            out[t] = Some(values[best]);
        }
        return out;
    }

    // Linear: Delaunay triangulation + barycentric coordinates
    if n_points < 3 {
        // Degenerate: fall back to nearest
        return interpolate_multipliers(
            InterpMethod::Nearest,
            points_x,
            points_y,
            values,
            target_x,
            target_y,
        );
    }

    let triangles = delaunay_triangulation(points_x, points_y);

    for (t, (tx, ty)) in target_x.iter().zip(target_y.iter()).enumerate() {
        let mut found = None;
        for tri in &triangles {
            let (a, b, c) = (tri[0], tri[1], tri[2]);
            let (x1, y1) = (points_x[a], points_y[a]);
            let (x2, y2) = (points_x[b], points_y[b]);
            let (x3, y3) = (points_x[c], points_y[c]);
            let det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
            if det.abs() < 1e-12 {
                continue;
            }
            let l1 = ((y2 - y3) * (tx - x3) + (x3 - x2) * (ty - y3)) / det;
            let l2 = ((y3 - y1) * (tx - x3) + (x1 - x3) * (ty - y3)) / det;
            let l3 = 1.0 - l1 - l2;
            // Inside the triangle when all barycentric coordinates are >= 0
            if l1 >= -1e-12 && l2 >= -1e-12 && l3 >= -1e-12 {
                found = Some(l1 * values[a] + l2 * values[b] + l3 * values[c]);
                break;
            }
        }
        out[t] = found;
    }

    out
}

/// Interpolate the speed multipliers of one condition at arbitrary 3D points.
///
/// Python uses `scipy.interpolate.LinearNDInterpolator` in three dimensions
/// (floris/core/flow_field.py:generate_heterogeneous_wind_map). This
/// implementation reproduces that behaviour for the structured case: when the
/// user-defined locations form a complete rectangular grid (every
/// combination of the unique x, y and z values is present), the multiplier
/// surface is the trilinear interpolant, which is identical to scipy's
/// piecewise-linear tetrahedral interpolant whenever the data varies
/// affinely (as in the FLORIS examples).
///
/// For scattered (non-grid) locations the nearest neighbour is used; a
/// full 3D Delaunay tetrahedralisation is not implemented.
pub fn interpolate_multipliers_3d(
    interp_method: InterpMethod,
    points_x: &Array1,
    points_y: &Array1,
    points_z: &Array1,
    values: &Array1,
    target_x: &[Float],
    target_y: &[Float],
    target_z: &[Float],
) -> Vec<Option<Float>> {
    let n_points = points_x.len();
    assert_eq!(n_points, points_y.len());
    assert_eq!(n_points, points_z.len());
    assert_eq!(n_points, values.len());

    let mut out = vec![None; target_x.len()];
    if n_points == 0 {
        return out;
    }

    // Nearest neighbour (scipy NearestNDInterpolator)
    if matches!(interp_method, InterpMethod::Nearest) {
        for (t, ((tx, ty), tz)) in target_x
            .iter()
            .zip(target_y.iter())
            .zip(target_z.iter())
            .enumerate()
        {
            let mut best = 0usize;
            let mut best_d = Float::INFINITY;
            for p in 0..n_points {
                let dx = points_x[p] - tx;
                let dy = points_y[p] - ty;
                let dz = points_z[p] - tz;
                let d = dx * dx + dy * dy + dz * dz;
                if d < best_d {
                    best_d = d;
                    best = p;
                }
            }
            out[t] = Some(values[best]);
        }
        return out;
    }

    // Structured grid: trilinear interpolation
    let Some(grid) = structured_grid(points_x, points_y, points_z, values) else {
        // Scattered 3D locations: fall back to the nearest neighbour
        return interpolate_multipliers_3d(
            InterpMethod::Nearest,
            points_x,
            points_y,
            points_z,
            values,
            target_x,
            target_y,
            target_z,
        );
    };

    for (t, ((tx, ty), tz)) in target_x
        .iter()
        .zip(target_y.iter())
        .zip(target_z.iter())
        .enumerate()
    {
        let Some((ix, fx)) = bracket(&grid.xs, *tx) else {
            continue; // outside the hull: freestream (1.0)
        };
        let Some((iy, fy)) = bracket(&grid.ys, *ty) else {
            continue;
        };
        let Some((iz, fz)) = bracket(&grid.zs, *tz) else {
            continue;
        };

        let mut value = 0.0;
        for (dx, wx) in [(0, 1.0 - fx), (1, fx)] {
            for (dy, wy) in [(0, 1.0 - fy), (1, fy)] {
                for (dz, wz) in [(0, 1.0 - fz), (1, fz)] {
                    value += wx * wy * wz * grid.value(ix + dx, iy + dy, iz + dz);
                }
            }
        }
        out[t] = Some(value);
    }

    out
}

/// A complete rectangular grid of multiplier locations.
struct StructuredGrid {
    xs: Vec<Float>,
    ys: Vec<Float>,
    zs: Vec<Float>,
    /// Values indexed as [ix][iy][iz]
    values: Vec<Vec<Vec<Float>>>,
}

impl StructuredGrid {
    fn value(&self, ix: usize, iy: usize, iz: usize) -> Float {
        self.values[ix][iy][iz]
    }
}

/// Build the structured grid description of the locations, or `None` when the
/// locations do not form a complete rectangular grid.
fn structured_grid(
    points_x: &Array1,
    points_y: &Array1,
    points_z: &Array1,
    values: &Array1,
) -> Option<StructuredGrid> {
    let mut xs: Vec<Float> = points_x.iter().copied().collect();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs.dedup();
    let mut ys: Vec<Float> = points_y.iter().copied().collect();
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ys.dedup();
    let mut zs: Vec<Float> = points_z.iter().copied().collect();
    zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    zs.dedup();

    if xs.len() * ys.len() * zs.len() != points_x.len()
        || xs.len() < 2
        || ys.len() < 2
        || zs.len() < 2
    {
        return None;
    }

    let mut values_grid = vec![vec![vec![0.0; zs.len()]; ys.len()]; xs.len()];
    for p in 0..points_x.len() {
        let ix = xs.iter().position(|v| (v - points_x[p]).abs() < 1e-9)?;
        let iy = ys.iter().position(|v| (v - points_y[p]).abs() < 1e-9)?;
        let iz = zs.iter().position(|v| (v - points_z[p]).abs() < 1e-9)?;
        values_grid[ix][iy][iz] = values[p];
    }

    Some(StructuredGrid {
        xs,
        ys,
        zs,
        values: values_grid,
    })
}

/// Locate the enclosing cell of `value` in a sorted axis.
///
/// Returns the index of the lower cell bound and the fraction towards the
/// upper bound, or `None` when the value lies outside the axis range.
fn bracket(axis: &[Float], value: Float) -> Option<(usize, Float)> {
    if value < axis[0] - 1e-9 || value > axis[axis.len() - 1] + 1e-9 {
        return None;
    }
    let mut i = 0usize;
    while i + 2 < axis.len() && value > axis[i + 1] {
        i += 1;
    }
    let lo = axis[i];
    let hi = axis[i + 1];
    let f = if (hi - lo).abs() < 1e-12 {
        0.0
    } else {
        ((value - lo) / (hi - lo)).clamp(0.0, 1.0)
    };
    Some((i, f))
}

/// Bowyer-Watson Delaunay triangulation of a 2D point set.
///
/// Returns the triangles as index triples into the input arrays.
fn delaunay_triangulation(points_x: &Array1, points_y: &Array1) -> Vec<[usize; 3]> {
    let n = points_x.len();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    if n < 3 {
        return triangles;
    }

    // Super-triangle large enough to contain all points
    let (mut min_x, mut max_x) = (points_x[0], points_x[0]);
    let (mut min_y, mut max_y) = (points_y[0], points_y[0]);
    for p in 0..n {
        min_x = min_x.min(points_x[p]);
        max_x = max_x.max(points_x[p]);
        min_y = min_y.min(points_y[p]);
        max_y = max_y.max(points_y[p]);
    }
    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let dmax = dx.max(dy).max(1.0);
    let mid_x = 0.5 * (min_x + max_x);
    let mid_y = 0.5 * (min_y + max_y);

    // Vertices: the real points followed by the three super-triangle corners
    let mut vx: Vec<Float> = points_x.iter().copied().collect();
    let mut vy: Vec<Float> = points_y.iter().copied().collect();
    let super_start = n;
    vx.push(mid_x - 20.0 * dmax);
    vy.push(mid_y - dmax);
    vx.push(mid_x);
    vy.push(mid_y + 20.0 * dmax);
    vx.push(mid_x + 20.0 * dmax);
    vy.push(mid_y - dmax);

    // Start with one triangle covering the super-triangle
    let mut tris: Vec<[usize; 3]> = vec![[super_start, super_start + 1, super_start + 2]];

    for p in 0..n {
        // Find the "bad" triangles (whose circumcircle contains p)
        let mut bad: Vec<usize> = Vec::new();
        let mut keep: Vec<[usize; 3]> = Vec::new();
        for (ti, t) in tris.iter().enumerate() {
            if circumcircle_contains(&vx, &vy, t, points_x[p], points_y[p]) {
                bad.push(ti);
            } else {
                keep.push(*t);
            }
        }
        if bad.is_empty() {
            tris = keep;
            continue;
        }

        // Boundary edges of the bad region (edges not shared by two bad triangles)
        let mut edges: Vec<(usize, usize)> = Vec::new();
        for &ti in &bad {
            let t = tris[ti];
            for e in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let shared = bad.iter().filter(|&&tj| tj != ti).any(|&tj| {
                    let u = tris[tj];
                    (u[0] == e.0 && u[1] == e.1)
                        || (u[1] == e.0 && u[0] == e.1)
                        || (u[1] == e.0 && u[2] == e.1)
                        || (u[2] == e.0 && u[1] == e.1)
                        || (u[2] == e.0 && u[0] == e.1)
                        || (u[0] == e.0 && u[2] == e.1)
                });
                if !shared {
                    edges.push(e);
                }
            }
        }

        // Remove the bad triangles and re-triangulate the hole
        tris = keep;
        for (a, b) in edges {
            // Skip degenerate (collinear) triangles
            let orient =
                (vx[b] - vx[a]) * (points_y[p] - vy[a]) - (vx[p] - vx[a]) * (vy[b] - vy[a]);
            if orient.abs() < 1e-12 {
                continue;
            }
            tris.push([a, b, p]);
        }
    }

    // Drop the triangles that use super-triangle vertices
    for t in tris {
        if t[0] < n && t[1] < n && t[2] < n {
            triangles.push(t);
        }
    }

    triangles
}

/// True if the point (px, py) lies inside the circumcircle of the triangle.
///
/// Uses the standard signed-area formulation: the in-circle determinant is
/// compared with the triangle's orientation, so the test is correct for both
/// clockwise and counter-clockwise triangles.
fn circumcircle_contains(vx: &[Float], vy: &[Float], t: &[usize; 3], px: Float, py: Float) -> bool {
    let (a, b, c) = (t[0], t[1], t[2]);
    let (ax, ay) = (vx[a], vy[a]);
    let (bx, by) = (vx[b], vy[b]);
    let (cx, cy) = (vx[c], vy[c]);

    // Orientation of the triangle (twice the signed area)
    let orient = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
    if orient.abs() < 1e-14 {
        return false;
    }

    let (adx, ady) = (ax - px, ay - py);
    let (bdx, bdy) = (bx - px, by - py);
    let (cdx, cdy) = (cx - px, cy - py);

    let det = (adx * bdy - bdx * ady) * (cdx * cdx + cdy * cdy)
        + (bdx * cdy - cdx * bdy) * (adx * adx + ady * ady)
        + (cdx * ady - adx * cdy) * (bdx * bdx + bdy * bdy);

    if orient > 0.0 {
        det > 0.0
    } else {
        det < 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nearest_interpolation() {
        let px = Array1::from_vec(vec![0.0, 10.0]);
        let py = Array1::from_vec(vec![0.0, 0.0]);
        let vals = Array1::from_vec(vec![1.0, 2.0]);
        let out = interpolate_multipliers(
            InterpMethod::Nearest,
            &px,
            &py,
            &vals,
            &[0.0, 10.0],
            &[0.0, 0.0],
        );
        assert_eq!(out[0], Some(1.0));
        assert_eq!(out[1], Some(2.0));
    }

    #[test]
    fn test_linear_interpolation_inside_hull() {
        // Unit square, multipliers increasing with x
        let px = Array1::from_vec(vec![0.0, 10.0, 10.0, 0.0]);
        let py = Array1::from_vec(vec![0.0, 0.0, 10.0, 10.0]);
        let vals = Array1::from_vec(vec![1.0, 2.0, 2.0, 1.0]);
        let out = interpolate_multipliers(InterpMethod::Linear, &px, &py, &vals, &[5.0], &[5.0]);
        let v = out[0].expect("point inside the hull must have a value");
        assert!((v - 1.5).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn test_linear_interpolation_outside_hull() {
        let px = Array1::from_vec(vec![0.0, 10.0, 10.0, 0.0]);
        let py = Array1::from_vec(vec![0.0, 0.0, 10.0, 10.0]);
        let vals = Array1::from_vec(vec![1.0, 2.0, 2.0, 1.0]);
        let out = interpolate_multipliers(InterpMethod::Linear, &px, &py, &vals, &[20.0], &[5.0]);
        assert!(out[0].is_none(), "point outside the hull must be None");
    }
}
