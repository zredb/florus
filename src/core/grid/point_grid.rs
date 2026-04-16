/// Points grid for arbitrary point calculations
use super::*;
use crate::types::{Array1, Array2, Array4, Float};
use crate::utilities::rotate_coordinates_rel_west;
use ndarray::Array;

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

impl Grid for PointsGrid {
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
        // PointsGrid is used for scattered point evaluations
        static INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
        INDICES.get_or_init(|| Array2::zeros((0, 0)))
    }
    fn sorted_coord_indices(&self) -> &Array2 {
        // PointsGrid doesn't use turbine sorting, return a static empty array
        static COORD_INDICES: std::sync::OnceLock<Array2> = std::sync::OnceLock::new();
        COORD_INDICES.get_or_init(|| Array2::zeros((0, 0)))
    }
    fn resolution(&self) -> usize {
        1
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
}
