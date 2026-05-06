/// Flow field planar grid for 2D slice visualization
use super::*;
use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::{reverse_rotate_coordinates_rel_west, rotate_coordinates_rel_west};
use ndarray::Array;
use std::any::Any;

#[derive(Debug, Clone)]
pub struct FlowFieldPlanarGrid {
    pub turbine_coordinates: Array2,
    pub turbine_diameters: Array1,
    pub turbine_hub_heights: Array1,  // Added to store hub heights
    pub wind_directions: Array1,
    pub grid_resolution: [usize; 2], // [nx1, nx2]
    pub normal_vector: String, // "x", "y", or "z"
    pub planar_coordinate: Float,
    pub x1_bounds: Option<(Float, Float)>,
    pub x2_bounds: Option<(Float, Float)>,

    pub n_turbines: usize,
    pub n_findex: usize,
    pub x_sorted: Array4,
    pub y_sorted: Array4,
    pub z_sorted: Array4,
    pub x_sorted_inertial_frame: Array4,
    pub y_sorted_inertial_frame: Array4,
    pub z_sorted_inertial_frame: Array4,
    pub x_center_of_rotation: Array2,
    pub y_center_of_rotation: Array2,
    pub sorted_indices: Array2,
    pub unsorted_indices: Array2,
}

impl FlowFieldPlanarGrid {
    pub fn new(
        turbine_coordinates: Array2,
        turbine_diameters: Array1,
        turbine_hub_heights: Array1,
        wind_directions: Array1,
        grid_resolution: [usize; 2],
        normal_vector: String,
        planar_coordinate: Float,
        x1_bounds: Option<(Float, Float)>,
        x2_bounds: Option<(Float, Float)>,
    ) -> crate::Result<Self> {
        let n_turbines = turbine_coordinates.shape()[0];
        let n_findex = wind_directions.len();

        // Rotate coordinates based on wind direction
        let (x, y, z, x_center_of_rotation, y_center_of_rotation) =
            rotate_coordinates_rel_west(&wind_directions, &turbine_coordinates)?;

        let max_diameter = turbine_diameters.iter().cloned().fold(Float::NEG_INFINITY, Float::max);

        // Determine bounds based on normal vector
        let (x1_min, x1_max) = if let Some(bounds) = x1_bounds {
            bounds
        } else {
            match normal_vector.as_str() {
                "z" => {
                    // Horizontal plane: x1 is x-direction
                    let min_x = x.iter().cloned().fold(Float::INFINITY, Float::min);
                    let max_x = x.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
                    (min_x - 2.0 * max_diameter, max_x + 10.0 * max_diameter)
                }
                "x" => {
                    // Cross plane: x1 is y-direction
                    let min_y = y.iter().cloned().fold(Float::INFINITY, Float::min);
                    let max_y = y.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
                    (min_y - 2.0 * max_diameter, max_y + 2.0 * max_diameter)
                }
                "y" => {
                    // Y-plane: x1 is x-direction
                    let min_x = x.iter().cloned().fold(Float::INFINITY, Float::min);
                    let max_x = x.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
                    (min_x - 2.0 * max_diameter, max_x + 10.0 * max_diameter)
                }
                _ => {
                    return Err(anyhow::anyhow!(
                        "Invalid normal_vector: {}. Must be 'x', 'y', or 'z'",
                        normal_vector
                    ));
                }
            }
        };

        let (x2_min, x2_max) = if let Some(bounds) = x2_bounds {
            bounds
        } else {
            match normal_vector.as_str() {
                "z" => {
                    // Horizontal plane: x2 is y-direction
                    let min_y = y.iter().cloned().fold(Float::INFINITY, Float::min);
                    let max_y = y.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
                    (min_y - 2.0 * max_diameter, max_y + 2.0 * max_diameter)
                }
                "x" => {
                    // Cross plane: x2 is z-direction
                    // Use turbine_hub_heights to determine z range
                    let min_height = turbine_hub_heights.iter().cloned().fold(Float::INFINITY, Float::min);
                    let max_height = turbine_hub_heights.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
                    (min_height - 2.0 * max_diameter, max_height + 2.0 * max_diameter)
                }
                "y" => {
                    // Y-plane: x2 is z-direction
                    // Use turbine_hub_heights to determine z range
                    let min_height = turbine_hub_heights.iter().cloned().fold(Float::INFINITY, Float::min);
                    let max_height = turbine_hub_heights.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
                    (min_height - 2.0 * max_diameter, max_height + 2.0 * max_diameter)
                }
                _ => {
                    return Err(anyhow::anyhow!(
                        "Invalid normal_vector: {}. Must be 'x', 'y', or 'z'",
                        normal_vector
                    ));
                }
            }
        };

        let nx1 = grid_resolution[0];
        let nx2 = grid_resolution[1];

        // Create meshgrid based on normal vector
        let (mut x_points, mut y_points, mut z_points) = match normal_vector.as_str() {
            "z" => {
                // Horizontal plane (x-y plane at fixed z)
                // Python creates 3 z-planes: [planar_coord - 10, planar_coord, planar_coord + 10]
                let nz = 3;
                let x1_vals = Array::linspace(x1_min, x1_max, nx1);
                let x2_vals = Array::linspace(x2_min, x2_max, nx2);
                let z_vals = vec![
                    planar_coordinate - 10.0,
                    planar_coordinate,
                    planar_coordinate + 10.0,
                ];

                let mut x_grid = Array::zeros((nx1, nx2, nz));
                let mut y_grid = Array::zeros((nx1, nx2, nz));
                let mut z_grid = Array::zeros((nx1, nx2, nz));

                for i in 0..nx1 {
                    for j in 0..nx2 {
                        for k in 0..nz {
                            x_grid[[i, j, k]] = x1_vals[i];
                            y_grid[[i, j, k]] = x2_vals[j];
                            z_grid[[i, j, k]] = z_vals[k];
                        }
                    }
                }

                (x_grid, y_grid, z_grid)
            }
            "x" => {
                // Cross plane (y-z plane at fixed x)
                let x1_vals = Array::linspace(x1_min, x1_max, nx1); // y values
                let x2_vals = Array::linspace(x2_min, x2_max, nx2); // z values

                let mut x_grid = Array::zeros((1, nx1, nx2));
                let mut y_grid = Array::zeros((1, nx1, nx2));
                let mut z_grid = Array::zeros((1, nx1, nx2));

                for j in 0..nx1 {
                    for k in 0..nx2 {
                        x_grid[[0, j, k]] = planar_coordinate;
                        y_grid[[0, j, k]] = x1_vals[j];
                        z_grid[[0, j, k]] = x2_vals[k];
                    }
                }

                (x_grid, y_grid, z_grid)
            }
            "y" => {
                // Y-plane (x-z plane at fixed y)
                let x1_vals = Array::linspace(x1_min, x1_max, nx1); // x values
                let x2_vals = Array::linspace(x2_min, x2_max, nx2); // z values

                let mut x_grid = Array::zeros((nx1, 1, nx2));
                let mut y_grid = Array::zeros((nx1, 1, nx2));
                let mut z_grid = Array::zeros((nx1, 1, nx2));

                for i in 0..nx1 {
                    for k in 0..nx2 {
                        x_grid[[i, 0, k]] = x1_vals[i];
                        y_grid[[i, 0, k]] = planar_coordinate;
                        z_grid[[i, 0, k]] = x2_vals[k];
                    }
                }

                (x_grid, y_grid, z_grid)
            }
            _ => {
                return Err(anyhow::anyhow!(
                    "Invalid normal_vector: {}. Must be 'x', 'y', or 'z'",
                    normal_vector
                ));
            }
        };

        // Reshape to add batch dimension for findex
        // Python uses [None, :, :, :] which adds one dimension at the front
        // Shape becomes (n_findex, nx1, nx2, nz_or_1)
        
        // For simplicity and consistency with Grid trait, reshape to (n_findex, total_points, 1, 1)
        let total_points = x_points.len();
        
        // Broadcast to n_findex
        let x_sorted = if n_findex == 1 {
            x_points.into_shape_with_order((1, total_points, 1, 1))?
        } else {
            // Repeat for each findex
            let mut expanded = Array::zeros((n_findex, total_points, 1, 1));
            for fi in 0..n_findex {
                let flat = x_points.clone().into_shape_with_order((total_points,))?;
                for p in 0..total_points {
                    expanded[[fi, p, 0, 0]] = flat[p];
                }
            }
            expanded
        };

        let y_sorted = if n_findex == 1 {
            y_points.into_shape_with_order((1, total_points, 1, 1))?
        } else {
            let mut expanded = Array::zeros((n_findex, total_points, 1, 1));
            for fi in 0..n_findex {
                let flat = y_points.clone().into_shape_with_order((total_points,))?;
                for p in 0..total_points {
                    expanded[[fi, p, 0, 0]] = flat[p];
                }
            }
            expanded
        };

        let z_sorted = if n_findex == 1 {
            z_points.into_shape_with_order((1, total_points, 1, 1))?
        } else {
            let mut expanded = Array::zeros((n_findex, total_points, 1, 1));
            for fi in 0..n_findex {
                let flat = z_points.clone().into_shape_with_order((total_points,))?;
                for p in 0..total_points {
                    expanded[[fi, p, 0, 0]] = flat[p];
                }
            }
            expanded
        };

        // Sort turbines by x coordinate (for consistency with other grids)
        let mut sorted_indices = Array::zeros((n_findex, n_turbines));
        let mut unsorted_indices = Array::zeros((n_findex, n_turbines));

        for fi in 0..n_findex {
            let mut indices: Vec<(usize, Float)> = (0..n_turbines)
                .map(|i| (i, x[[fi, i]]))
                .collect();
            indices.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            for (new_i, (old_i, _)) in indices.iter().enumerate() {
                sorted_indices[[fi, new_i]] = *old_i as Float;
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
            turbine_hub_heights,
            wind_directions,
            grid_resolution,
            normal_vector,
            planar_coordinate,
            x1_bounds: Some((x1_min, x1_max)),
            x2_bounds: Some((x2_min, x2_max)),
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
            sorted_indices,
            unsorted_indices,
        })
    }
}

impl Grid for FlowFieldPlanarGrid {
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
        // Return identity mapping for planar grid
        static COORD_INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
        COORD_INDICES.get_or_init(|| Array2::zeros((0, 0)))
    }
    fn resolution(&self) -> usize {
        self.grid_resolution[0]
    }
    fn hub_heights(&self) -> Array1 {
        // Return the stored turbine hub heights
        self.turbine_hub_heights.clone()
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
}
