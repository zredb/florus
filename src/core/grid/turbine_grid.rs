// Standard turbine grid with square rotor points
use super::*;
use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::{reverse_rotate_coordinates_rel_west, rotate_coordinates_rel_west};
use ndarray::Array;

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
    pub average_method: AveragingMethod,
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

        let mut x_sorted = Array::zeros((n_findex, n_turbines, grid_resolution, grid_resolution));
        let mut y_sorted = Array::zeros((n_findex, n_turbines, grid_resolution, grid_resolution));
        let mut z_sorted = Array::zeros((n_findex, n_turbines, grid_resolution, grid_resolution));

        for fi in 0..n_findex {
            for ti in 0..n_turbines {
                let original_idx = sorted_indices[[fi, ti]] as usize;
                for i in 0..grid_resolution {
                    for j in 0..grid_resolution {
                        x_sorted[[fi, ti, i, j]] = _x[[fi, original_idx, i, j]];
                        y_sorted[[fi, ti, i, j]] = _y[[fi, original_idx, i, j]];
                        z_sorted[[fi, ti, i, j]] = _z[[fi, original_idx, i, j]];
                    }
                }
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
            sorted_indices,
            sorted_coord_indices,
            unsorted_indices,
            x_center_of_rotation,
            y_center_of_rotation,
            average_method: AveragingMethod::CubicMean,
        })
    }
}

impl Grid for TurbineGrid {
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
    fn sorted_coord_indices(&self) -> &Array2 {
        &self.sorted_coord_indices
    }
    fn unsorted_indices(&self) -> Option<&Array2> {
        Some(&self.unsorted_indices)
    }
    fn resolution(&self) -> usize {
        self.grid_resolution
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
}
