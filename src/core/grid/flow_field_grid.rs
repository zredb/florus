/// Flow field grid for full domain visualization
use super::*;
use crate::utilities::{reverse_rotate_coordinates_rel_west, rotate_coordinates_rel_west};
use ndarray::Array;
use std::any::Any;

#[derive(Debug, Clone)]
pub struct FlowFieldGrid {
    pub turbine_coordinates: crate::Array2,
    pub turbine_diameters: crate::Array1,
    pub wind_directions: crate::Array1,
    pub grid_resolution: [usize; 3], // [nx, ny, nz]

    pub n_turbines: usize,
    pub n_findex: usize,
    pub x_sorted: crate::Array4,
    pub y_sorted: crate::Array4,
    pub z_sorted: crate::Array4,
    pub x_sorted_inertial_frame: crate::Array4,
    pub y_sorted_inertial_frame: crate::Array4,
    pub z_sorted_inertial_frame: crate::Array4,
    pub x_center_of_rotation: crate::Array2,
    pub y_center_of_rotation: crate::Array2,
}

impl FlowFieldGrid {
    pub fn new(
        turbine_coordinates: crate::Array2,
        turbine_diameters: crate::Array1,
        wind_directions: crate::Array1,
        grid_resolution: usize,
    ) -> crate::Result<Self> {
        // Convert single resolution to [usize; 3] tuple
        let grid_resolution = [grid_resolution, grid_resolution, grid_resolution];
        Self::new_with_resolution(turbine_coordinates, turbine_diameters, wind_directions, grid_resolution)
    }

    pub fn new_with_resolution(
        turbine_coordinates: crate::Array2,
        turbine_diameters: crate::Array1,
        wind_directions: crate::Array1,
        grid_resolution: [usize; 3],
    ) -> crate::Result<Self> {
        let n_turbines = turbine_coordinates.shape()[0];
        let n_findex = wind_directions.len();

        // Rotate coordinates based on wind direction
        let (x, y, z, x_center_of_rotation, y_center_of_rotation) =
            rotate_coordinates_rel_west(&wind_directions, &turbine_coordinates)?;

        // Calculate domain bounds
        let eps = 0.01;
        let max_diameter = turbine_diameters.iter().cloned().fold(Float::NEG_INFINITY, Float::max);

        // Use first findex for bounds calculation (assuming all turbines visible in all conditions)
        let xmin = x.column(0).iter().cloned().fold(Float::INFINITY, Float::min) - 2.0 * max_diameter;
        let xmax = x.column(0).iter().cloned().fold(Float::NEG_INFINITY, Float::max) + 10.0 * max_diameter;
        let ymin = y.column(0).iter().cloned().fold(Float::INFINITY, Float::min) - 2.0 * max_diameter;
        let ymax = y.column(0).iter().cloned().fold(Float::NEG_INFINITY, Float::max) + 2.0 * max_diameter;
        let zmin = eps;
        let zmax = 6.0 * z.column(0).iter().cloned().fold(Float::NEG_INFINITY, Float::max);

        let nx = grid_resolution[0];
        let ny = grid_resolution[1];
        let nz = grid_resolution[2];

        // Create meshgrid using linspace
        let x_vals = Array::linspace(xmin, xmax, nx);
        let y_vals = Array::linspace(ymin, ymax, ny);
        let z_vals = Array::linspace(zmin, zmax, nz);

        // Create 3D meshgrid with indexing="ij"
        // Shape: (nx, ny, nz)
        let mut x_points = Array::zeros((nx, ny, nz));
        let mut y_points = Array::zeros((nx, ny, nz));
        let mut z_points = Array::zeros((nx, ny, nz));

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    x_points[[i, j, k]] = x_vals[i];
                    y_points[[i, j, k]] = y_vals[j];
                    z_points[[i, j, k]] = z_vals[k];
                }
            }
        }

        // Reshape to add batch dimensions: (1, 1, nx, ny, nz) -> but we use Array4
        // Python uses [None, None, :, :, :] which adds two dimensions at the front
        // We'll store as (1, nx*ny*nz, 1, 1) or keep as is and reshape when needed
        // For compatibility with Grid trait, we reshape to (n_findex, n_points, 1, 1)
        // But flow field grid doesn't have turbine dimension, so we use (1, total_points, 1, 1)
        
        let total_points = nx * ny * nz;
        
        // Reshape to (1, total_points, 1, 1) for consistency with Grid trait
        // This represents: (dummy_findex, point_index, dummy_dim1, dummy_dim2)
        let x_sorted = x_points.clone().into_shape_with_order((1, total_points, 1, 1))?;
        let y_sorted = y_points.clone().into_shape_with_order((1, total_points, 1, 1))?;
        let z_sorted = z_points.clone().into_shape_with_order((1, total_points, 1, 1))?;

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
            x_center_of_rotation,
            y_center_of_rotation,
        })
    }
}

impl Grid for FlowFieldGrid {
    fn n_turbines(&self) -> usize {
        self.n_turbines
    }
    fn n_findex(&self) -> usize {
        self.n_findex
    }
    fn x_sorted(&self) -> &crate::Array4 {
        &self.x_sorted
    }
    fn y_sorted(&self) -> &crate::Array4 {
        &self.y_sorted
    }
    fn z_sorted(&self) -> &crate::Array4 {
        &self.z_sorted
    }
    fn x_sorted_inertial_frame(&self) -> &crate::Array4 {
        &self.x_sorted_inertial_frame
    }
    fn y_sorted_inertial_frame(&self) -> &crate::Array4 {
        &self.y_sorted_inertial_frame
    }
    fn z_sorted_inertial_frame(&self) -> &crate::Array4 {
        &self.z_sorted_inertial_frame
    }
    fn average_method(&self) -> AveragingMethod {
        AveragingMethod::CubicMean
    }
    fn sorted_indices(&self) -> &crate::Array2 {
        // FlowFieldGrid doesn't use turbine sorting
        static INDICES: std::sync::OnceLock<crate::Array2> = std::sync::OnceLock::new();
        INDICES.get_or_init(|| crate::Array2::zeros((0, 0)))
    }
    fn sorted_coord_indices(&self) -> &crate::Array2 {
        static COORD_INDICES: std::sync::OnceLock<crate::Array2> = std::sync::OnceLock::new();
        COORD_INDICES.get_or_init(|| crate::Array2::zeros((0, 0)))
    }
    fn resolution(&self) -> usize {
        // Return the product of all resolutions or just the first one
        self.grid_resolution[0]
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
}
