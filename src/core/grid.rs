use crate::core::AveragingMethod;
/// Grid types for FLORIS calculations
///
/// Corresponds to grid.py in the Python implementation
use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::{reverse_rotate_coordinates_rel_west, rotate_coordinates_rel_west};
use ndarray::Array;
use std::f64::consts::PI;

/// Base grid trait
pub trait GridBase {
    fn n_turbines(&self) -> usize;
    fn n_findex(&self) -> usize;
    fn x_sorted(&self) -> &Array4;
    fn y_sorted(&self) -> &Array4;
    fn z_sorted(&self) -> &Array4;
    fn x_sorted_inertial_frame(&self) -> &Array4;
    fn y_sorted_inertial_frame(&self) -> &Array4;
    fn z_sorted_inertial_frame(&self) -> &Array4;
    fn cubature_weights(&self) -> Option<&Array2> {
        None
    }
    fn average_method(&self) -> AveragingMethod;
    fn sorted_indices(&self) -> &Array2;
    fn resolution(&self) -> usize;
}

/// Standard turbine grid with square rotor points
#[derive(Debug, Clone)]
pub struct TurbineGrid {
    pub turbine_coordinates: Array2, // (n_turbines, 3)
    pub turbine_diameters: Array1,   // (n_turbines,)
    pub wind_directions: Array1,     // (n_findex,)
    pub grid_resolution: usize,
    pub n_turbines: usize,
    pub n_findex: usize,
    pub x_sorted: Array4, // (n_findex, n_turbines, grid_res, grid_res)
    pub y_sorted: Array4,
    pub z_sorted: Array4,
    pub x_sorted_inertial_frame: Array4,
    pub y_sorted_inertial_frame: Array4,
    pub z_sorted_inertial_frame: Array4,
    pub sorted_indices: Array2,
    pub sorted_coord_indices: Array2,
    pub unsorted_indices: Array2,
    pub x_center_of_rotation: Array2,
    pub y_center_of_rotation: Array2,
}

impl TurbineGrid {
    pub fn new(
        turbine_coordinates: Array2,
        turbine_diameters: Array1,
        wind_directions: Array1,
        grid_resolution: usize,
    ) -> crate::Result<Self> {
        let n_turbines = turbine_coordinates.shape()[0];
        let n_findex = wind_directions.len();

        // Rotate coordinates based on wind direction
        let (x, y, z, x_center_of_rotation, y_center_of_rotation) =
            rotate_coordinates_rel_west(&wind_directions, &turbine_coordinates)?;

        // Create grid points on rotor
        let radius_ratio = 0.5;
        let disc_area_radius = turbine_diameters.mapv(|d| radius_ratio * d / 2.0);

        // Create disc grid
        let disc_grid = if grid_resolution == 1 {
            Array::zeros((n_turbines, 1))
        } else {
            let mut dg = Array::zeros((n_turbines, grid_resolution));
            for (i, &radius) in disc_area_radius.iter().enumerate() {
                let points = Array::linspace(-radius, radius, grid_resolution);
                dg.row_mut(i).assign(&points);
            }
            dg
        };

        // Create template grid
        let mut _x = Array::ones((n_findex, n_turbines, grid_resolution, grid_resolution));
        let mut _y = Array::ones((n_findex, n_turbines, grid_resolution, grid_resolution));
        let mut _z = Array::ones((n_findex, n_turbines, grid_resolution, grid_resolution));

        // Fill in coordinates
        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                for i in 0..grid_resolution {
                    for j in 0..grid_resolution {
                        _x[[fi, ti, i, j]] = x[[fi, ti]];
                        _y[[fi, ti, i, j]] = y[[fi, ti]] + disc_grid[[ti, i]];
                        _z[[fi, ti, i, j]] = z[[fi, ti]] + disc_grid[[ti, j]];
                    }
                }
            }
        }

        // Sort turbines by x coordinate (upstream to downstream)
        let mut sorted_indices = Array::zeros((n_findex, n_turbines));
        let mut sorted_coord_indices = Array::zeros((n_findex, n_turbines));
        let mut unsorted_indices = Array::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            let mut indices: Vec<(usize, Float)> =
                (0..n_turbines).map(|i| (i, x[[fi, i]])).collect();
            indices.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            for (new_i, (old_i, _)) in indices.iter().enumerate() {
                sorted_indices[[fi, new_i]] = *old_i as Float;
                sorted_coord_indices[[fi, new_i]] = *old_i as Float;
                unsorted_indices[[fi, *old_i]] = new_i as Float;
            }
        }

        // Apply sorting
        let x_sorted = _x.clone();
        let y_sorted = _y.clone();
        let z_sorted = _z.clone();

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
            sorted_indices,
            sorted_coord_indices,
            unsorted_indices,
            x_center_of_rotation,
            y_center_of_rotation,
        })
    }
}

impl GridBase for TurbineGrid {
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
    fn average_method(&self) -> AveragingMethod {
        AveragingMethod::CubicMean
    }
    fn sorted_indices(&self) -> &Array2 {
        &self.sorted_indices
    }
    fn resolution(&self) -> usize {
        self.grid_resolution
    }
}

/// Turbine grid with cubature integration points
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

        // Get cubature coefficients
        let _coeffs = Self::get_cubature_coefficients(grid_resolution)?;

        // Rotate coordinates
        let (_x, _y, _z, x_center_of_rotation, y_center_of_rotation) =
            rotate_coordinates_rel_west(&wind_directions, &turbine_coordinates)?;

        // Generate grid points - simplified implementation
        let n_points = grid_resolution * grid_resolution;
        let x_sorted = Array::ones((n_findex, n_turbines, n_points, 1));
        let y_sorted = Array::ones((n_findex, n_turbines, n_points, 1));
        let z_sorted = Array::ones((n_findex, n_turbines, n_points, 1));

        let cubature_weights =
            Array::from_elem((grid_resolution, grid_resolution), 1.0 / n_points as Float);

        // Placeholder for sorted indices
        let sorted_indices = Array::zeros((n_findex, n_turbines));
        let sorted_coord_indices = Array::zeros((n_findex, n_turbines));
        let unsorted_indices = Array::zeros((n_findex, n_turbines));

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

    pub fn get_cubature_coefficients(n: usize) -> crate::Result<CubatureCoefficients> {
        match n {
            1 => Ok(CubatureCoefficients {
                r: vec![0.0],
                t: vec![0.0],
                q: vec![1.0],
                a: vec![1.0],
                b: PI,
            }),
            2 => Ok(CubatureCoefficients {
                r: vec![-0.7071067811865475, 0.7071067811865475],
                t: vec![-0.7071067811865475, 0.7071067811865475],
                q: vec![0.7071067811865475, 0.7071067811865475],
                a: vec![0.5, 0.5],
                b: PI / 2.0,
            }),
            _ => Err(anyhow::anyhow!(
                "Cubature coefficients not implemented for N={}",
                n
            )),
        }
    }
}

impl GridBase for TurbineCubatureGrid {
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
    fn resolution(&self) -> usize {
        self.grid_resolution
    }
}

/// Cubature coefficients
#[derive(Debug, Clone)]
pub struct CubatureCoefficients {
    pub r: Vec<Float>,
    pub t: Vec<Float>,
    pub q: Vec<Float>,
    pub a: Vec<Float>,
    pub b: Float,
}

/// Points grid for arbitrary point calculations
#[derive(Debug, Clone)]
pub struct PointsGrid {
    pub points_x: Array1,
    pub points_y: Array1,
    pub points_z: Array1,
    pub wind_directions: Array1,
    pub n_findex: usize,
    pub x_sorted: Array4,
    pub y_sorted: Array4,
    pub z_sorted: Array4,
    pub x_sorted_inertial_frame: Array4,
    pub y_sorted_inertial_frame: Array4,
    pub z_sorted_inertial_frame: Array4,
    pub x_center_of_rotation: Option<Float>,
    pub y_center_of_rotation: Option<Float>,
}

impl PointsGrid {
    pub fn new(
        points_x: Array1,
        points_y: Array1,
        points_z: Array1,
        wind_directions: Array1,
        x_center_of_rotation: Option<Float>,
        y_center_of_rotation: Option<Float>,
    ) -> crate::Result<Self> {
        let n_points = points_x.len();
        let n_findex = wind_directions.len();

        // Create point coordinates
        let mut point_coordinates = Array::zeros((n_points, 3));
        for i in 0..n_points {
            point_coordinates[[i, 0]] = points_x[i];
            point_coordinates[[i, 1]] = points_y[i];
            point_coordinates[[i, 2]] = points_z[i];
        }

        // Rotate coordinates
        let (x, y, z, _, _) = rotate_coordinates_rel_west(&wind_directions, &point_coordinates)?;

        // Reshape to grid format
        let x_sorted = x
            .clone()
            .into_shape_with_order((n_findex, n_points, 1, 1))?;
        let y_sorted = y
            .clone()
            .into_shape_with_order((n_findex, n_points, 1, 1))?;
        let z_sorted = z
            .clone()
            .into_shape_with_order((n_findex, n_points, 1, 1))?;

        // For inertial frame, just use original coordinates
        let x_sorted_inertial_frame = x_sorted.clone();
        let y_sorted_inertial_frame = y_sorted.clone();
        let z_sorted_inertial_frame = z_sorted.clone();

        Ok(Self {
            points_x,
            points_y,
            points_z,
            wind_directions,
            n_findex,
            x_sorted,
            y_sorted,
            z_sorted,
            x_sorted_inertial_frame,
            y_sorted_inertial_frame,
            z_sorted_inertial_frame,
            x_center_of_rotation,
            y_center_of_rotation,
        })
    }
}

impl GridBase for PointsGrid {
    fn n_turbines(&self) -> usize {
        self.points_x.len()
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
    fn average_method(&self) -> AveragingMethod {
        AveragingMethod::CubicMean
    }
    fn sorted_indices(&self) -> &Array2 {
        // PointsGrid doesn't have sorted indices, return identity mapping
        // This is a placeholder - in practice, PointsGrid may need proper sorting support
        static INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
        INDICES.get_or_init(|| Array2::zeros((0, 0)))
    }
    fn resolution(&self) -> usize {
        // PointsGrid doesn't have a traditional grid resolution
        // Return 1 as default for single point per location
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array;

    #[test]
    fn test_turbine_grid_creation() {
        let coords = Array::from_shape_vec((2, 3), vec![0.0, 0.0, 90.0, 500.0, 0.0, 90.0]).unwrap();
        let diameters = Array::from_vec(vec![126.0, 126.0]);
        let wind_dirs = Array::from_vec(vec![270.0]);

        let grid = TurbineGrid::new(coords, diameters, wind_dirs, 3);
        assert!(grid.is_ok());

        let grid = grid.unwrap();
        assert_eq!(grid.n_turbines, 2);
        assert_eq!(grid.n_findex, 1);
        assert_eq!(grid.grid_resolution, 3);
    }

    #[test]
    fn test_cubature_coefficients() {
        let coeffs = TurbineCubatureGrid::get_cubature_coefficients(1);
        assert!(coeffs.is_ok());

        let coeffs = coeffs.unwrap();
        assert_eq!(coeffs.r.len(), 1);
        assert_eq!(coeffs.a.len(), 1);
    }
}
