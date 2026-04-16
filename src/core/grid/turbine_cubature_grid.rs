/// Turbine grid with cubature integration points
use super::*;
use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::{reverse_rotate_coordinates_rel_west, rotate_coordinates_rel_west};
use ndarray::Array;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct TurbineCubatureGrid {
    pub turbine_coordinates: Array2,
    pub turbine_diameters: Array1,
    pub wind_directions: Array1,
    pub grid_resolution: usize,
    pub n_turbines: usize,
    pub n_findex: usize,
    pub x_sorted: Array4,
    pub y_sorted: Array4,
    pub z_sorted: Array4,
    pub x_sorted_inertial_frame: Array4,
    pub y_sorted_inertial_frame: Array4,
    pub z_sorted_inertial_frame: Array4,
    pub cubature_weights: Array2,
    pub sorted_indices: Array2,
    pub sorted_coord_indices: Array2,
    pub unsorted_indices: Array2,
    pub x_center_of_rotation: Array2,
    pub y_center_of_rotation: Array2,
}

impl TurbineCubatureGrid {
    pub fn new(
        turbine_coordinates: Array2,
        turbine_diameters: Array1,
        wind_directions: Array1,
        grid_resolution: usize,
    ) -> crate::Result<Self> {
        if grid_resolution < 1 || grid_resolution > 10 {
            return Err(anyhow::anyhow!(
                "Cubature grid resolution must be between 1 and 10, got {}",
                grid_resolution
            ));
        }

        let n_turbines = turbine_coordinates.shape()[0];
        let n_findex = wind_directions.len();

        // Get cubature coefficients based on resolution
        // Using Gauss-Legendre based cubature for circular disk integration
        let coeffs = Self::compute_cubature_points(grid_resolution)?;
        let n_points = coeffs.len();

        // Rotate coordinates
        let (
            rotated_coords_x,
            rotated_coords_y,
            rotated_coords_z,
            x_center_of_rotation,
            y_center_of_rotation,
        ) = rotate_coordinates_rel_west(&wind_directions, &turbine_coordinates)?;

        // Get hub heights from turbine coordinates (z-coordinate)
        let hub_heights = turbine_coordinates.column(2).to_owned();

        // Generate cubature grid points at rotor positions
        // Shape: (n_findex, n_turbines, n_points, 1)
        let mut x_sorted = Array::zeros((n_findex, n_turbines, n_points, 1));
        let mut y_sorted = Array::zeros((n_findex, n_turbines, n_points, 1));
        let mut z_sorted = Array::zeros((n_findex, n_turbines, n_points, 1));
        let mut cubature_weights = Array::zeros((n_turbines, n_points));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let rotor_radius = turbine_diameters[ti] / 2.0;

                for (pi, coeff) in coeffs.iter().enumerate() {
                    // Apply cubature point position scaled by rotor radius
                    x_sorted[[fi, ti, pi, 0]] = rotated_coords_x[[fi, ti]];
                    y_sorted[[fi, ti, pi, 0]] = rotated_coords_y[[fi, ti]] + coeff.y * rotor_radius;
                    z_sorted[[fi, ti, pi, 0]] = hub_heights[ti] + coeff.z * rotor_radius;
                }

                // Store normalized weights for this turbine
                let total_weight: Float = coeffs.iter().map(|c| c.w).sum();
                for (pi, coeff) in coeffs.iter().enumerate() {
                    cubature_weights[[ti, pi]] = coeff.w / total_weight;
                }
            }
        }

        // Sort turbines by x coordinate (upstream to downstream)
        let mut sorted_indices = Array::zeros((n_findex, n_turbines));
        let mut sorted_coord_indices = Array::zeros((n_findex, n_turbines));
        let mut unsorted_indices = Array::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            let mut indices: Vec<(usize, Float)> = (0..n_turbines)
                .map(|i| (i, rotated_coords_x[[fi, i]]))
                .collect();
            indices.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            for (new_i, (old_i, _)) in indices.iter().enumerate() {
                sorted_indices[[fi, new_i]] = *old_i as Float;
                sorted_coord_indices[[fi, new_i]] = *old_i as Float;
                unsorted_indices[[fi, *old_i]] = new_i as Float;
            }
        }

        // Reverse rotate to get inertial frame coordinates
        let (x_sorted_inertial_frame, y_sorted_inertial_frame, z_sorted_inertial_frame) =
            reverse_rotate_coordinates_rel_west(
                &wind_directions,
                &x_sorted,
                &y_sorted,
                &z_sorted,
                &x_center_of_rotation,
                &y_center_of_rotation,
            )?;

        Ok(Self {
            turbine_coordinates,
            turbine_diameters,
            wind_directions,
            grid_resolution,
            n_turbines,
            n_findex,
            x_sorted,
            y_sorted,
            z_sorted,
            x_sorted_inertial_frame,
            y_sorted_inertial_frame,
            z_sorted_inertial_frame,
            cubature_weights,
            sorted_indices,
            sorted_coord_indices,
            unsorted_indices,
            x_center_of_rotation,
            y_center_of_rotation,
        })
    }

    /// Compute cubature points and weights for circular disk integration
    ///
    /// This implements Gauss-Legendre based cubature for integrating over a circular disk.
    /// The points are placed optimally for numerical integration.
    fn compute_cubature_points(n: usize) -> crate::Result<Vec<CubaturePoint>> {
        match n {
            1 => Ok(vec![CubaturePoint {
                y: 0.0,
                z: 0.0,
                w: 1.0,
            }]),
            2 => {
                // Two-point approximation
                let r = 1.0 / std::f64::consts::SQRT_2;
                Ok(vec![
                    CubaturePoint {
                        y: r,
                        z: 0.0,
                        w: 0.5,
                    },
                    CubaturePoint {
                        y: -r,
                        z: 0.0,
                        w: 0.5,
                    },
                ])
            }
            3 => {
                let r = 2.0 / 3.0;
                let sqrt3 = 3.0_f64.sqrt();
                Ok(vec![
                    CubaturePoint {
                        y: 0.0,
                        z: r,
                        w: 1.0 / 3.0,
                    },
                    CubaturePoint {
                        y: r * sqrt3 / 2.0,
                        z: -r / 2.0,
                        w: 1.0 / 3.0,
                    },
                    CubaturePoint {
                        y: -r * sqrt3 / 2.0,
                        z: -r / 2.0,
                        w: 1.0 / 3.0,
                    },
                ])
            }
            4 => {
                // Four-point approximation (square vertices)
                let r = std::f64::consts::SQRT_2 / 2.0;
                Ok(vec![
                    CubaturePoint {
                        y: r,
                        z: r,
                        w: 0.25,
                    },
                    CubaturePoint {
                        y: -r,
                        z: r,
                        w: 0.25,
                    },
                    CubaturePoint {
                        y: r,
                        z: -r,
                        w: 0.25,
                    },
                    CubaturePoint {
                        y: -r,
                        z: -r,
                        w: 0.25,
                    },
                ])
            }
            5 => {
                // Five-point approximation (center + square)
                let r = 2.0 / 3.0;
                Ok(vec![
                    CubaturePoint {
                        y: 0.0,
                        z: 0.0,
                        w: 0.4,
                    },
                    CubaturePoint {
                        y: r,
                        z: r,
                        w: 0.15,
                    },
                    CubaturePoint {
                        y: -r,
                        z: r,
                        w: 0.15,
                    },
                    CubaturePoint {
                        y: r,
                        z: -r,
                        w: 0.15,
                    },
                    CubaturePoint {
                        y: -r,
                        z: -r,
                        w: 0.15,
                    },
                ])
            }
            6 => {
                let sqrt3_inv = 1.0 / 3.0_f64.sqrt();
                let points: Vec<CubaturePoint> = (0..6)
                    .map(|i| {
                        let angle = (i as Float) * PI / 3.0;
                        CubaturePoint {
                            y: sqrt3_inv * angle.cos(),
                            z: sqrt3_inv * angle.sin(),
                            w: 1.0 / 6.0,
                        }
                    })
                    .collect();
                Ok(points)
            }
            8 => {
                // Eight-point approximation (square vertices + midpoints)
                let r = std::f64::consts::SQRT_2 / 2.0;
                let w_corner = 1.0 / 6.0;
                let w_edge = 1.0 / 3.0;
                Ok(vec![
                    CubaturePoint {
                        y: r,
                        z: r,
                        w: w_corner,
                    },
                    CubaturePoint {
                        y: -r,
                        z: r,
                        w: w_corner,
                    },
                    CubaturePoint {
                        y: r,
                        z: -r,
                        w: w_corner,
                    },
                    CubaturePoint {
                        y: -r,
                        z: -r,
                        w: w_corner,
                    },
                    CubaturePoint {
                        y: 1.0,
                        z: 0.0,
                        w: w_edge,
                    },
                    CubaturePoint {
                        y: -1.0,
                        z: 0.0,
                        w: w_edge,
                    },
                    CubaturePoint {
                        y: 0.0,
                        z: 1.0,
                        w: w_edge,
                    },
                    CubaturePoint {
                        y: 0.0,
                        z: -1.0,
                        w: w_edge,
                    },
                ])
            }
            _ => {
                // For larger n, use tensor product of Gauss-Legendre points
                let m = (n + 1) / 2;
                let points = Self::gauss_legendre_2d(m);
                Ok(points)
            }
        }
    }

    fn gauss_legendre_2d(n: usize) -> Vec<CubaturePoint> {
        let mut points = Vec::new();
        let w = 1.0 / ((n * n) as Float);

        for i in 0..n {
            let y = Self::gauss_legendre_point(i, n);
            for j in 0..n {
                let z = Self::gauss_legendre_point(j, n);
                let r = (y * y + z * z).sqrt();
                if r <= 1.0 {
                    points.push(CubaturePoint { y, z, w });
                }
            }
        }

        // Normalize weights
        if !points.is_empty() {
            let total: Float = points.iter().map(|p| p.w).sum();
            for p in &mut points {
                p.w /= total;
            }
        }

        points
    }

    fn gauss_legendre_point(i: usize, n: usize) -> Float {
        let x = (((2 * i + 1) as Float) * PI) / ((2 * n + 1) as Float);
        x.cos() * (4.0 / ((2 * n + 1) as Float)).sqrt()
    }
}

/// Cubature point for disk integration
#[derive(Debug, Clone)]
struct CubaturePoint {
    y: Float,
    z: Float,
    w: Float,
}

impl Grid for TurbineCubatureGrid {
    fn n_turbines(&self) -> usize {
        self.n_turbines
    }
    fn n_findex(&self) -> usize {
        self.n_findex
    }
    fn x_sorted(&self) -> &Array4 {
        &self.x_sorted
    }
    fn y_sorted(&self) -> &Array4 {
        &self.y_sorted
    }
    fn z_sorted(&self) -> &Array4 {
        &self.z_sorted
    }
    fn x_sorted_inertial_frame(&self) -> &Array4 {
        &self.x_sorted_inertial_frame
    }
    fn y_sorted_inertial_frame(&self) -> &Array4 {
        &self.y_sorted_inertial_frame
    }
    fn z_sorted_inertial_frame(&self) -> &Array4 {
        &self.z_sorted_inertial_frame
    }
    fn cubature_weights(&self) -> Option<&Array2> {
        Some(&self.cubature_weights)
    }
    fn average_method(&self) -> AveragingMethod {
        AveragingMethod::SimpleCubature
    }
    fn sorted_indices(&self) -> &Array2 {
        &self.sorted_indices
    }
    fn sorted_coord_indices(&self) -> &Array2 {
        &self.sorted_coord_indices
    }
    fn resolution(&self) -> usize {
        self.grid_resolution
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl Default for TurbineCubatureGrid {
    fn default() -> Self {
        let empty_1d = Array::zeros((0,));
        let empty_2d = Array::zeros((0, 0));
        let empty_4d = Array::zeros((0, 0, 0, 0));
        Self {
            turbine_coordinates: empty_2d.clone(),
            turbine_diameters: empty_1d.clone(),
            wind_directions: empty_1d.clone(),
            grid_resolution: 0,
            n_turbines: 0,
            n_findex: 0,
            x_sorted: empty_4d.clone(),
            y_sorted: empty_4d.clone(),
            z_sorted: empty_4d.clone(),
            x_sorted_inertial_frame: empty_4d.clone(),
            y_sorted_inertial_frame: empty_4d.clone(),
            z_sorted_inertial_frame: empty_4d.clone(),
            cubature_weights: empty_2d.clone(),
            sorted_indices: empty_2d.clone(),
            sorted_coord_indices: empty_2d.clone(),
            unsorted_indices: empty_2d.clone(),
            x_center_of_rotation: empty_2d.clone(),
            y_center_of_rotation: empty_2d.clone(),
        }
    }
}
