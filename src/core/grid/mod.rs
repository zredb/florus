pub mod flow_field_grid;
pub mod flow_field_planar_grid;
pub mod point_grid;
pub mod turbine_cubature_grid;
pub mod turbine_grid;

use crate::core::AveragingMethod;
use crate::types::{Array1, Array2, Array4, Float};
use std::any::Any;

pub use flow_field_grid::FlowFieldGrid;
pub use flow_field_planar_grid::FlowFieldPlanarGrid;
pub use point_grid::PointsGrid;
pub use turbine_cubature_grid::TurbineCubatureGrid;
pub use turbine_grid::TurbineGrid;

pub trait Grid {
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
    fn sorted_coord_indices(&self) -> &Array2;
    fn unsorted_indices(&self) -> Option<&Array2> {
        None
    }
    fn unsorted_coord_indices(&self) -> Option<&Array2> {
        None
    }
    fn resolution(&self) -> usize;
    fn as_any(&self) -> &dyn Any;
    fn is_turbine_grid(&self) -> bool {
        self.as_any().is::<TurbineGrid>()
    }
    fn grid_shape(&self) -> (usize, usize, usize, usize) {
        let x = self.x_sorted();
        (x.shape()[0], x.shape()[1], x.shape()[2], x.shape()[3])
    }
    fn hub_heights(&self) -> Array1 {
        let n_turbines = self.n_turbines();
        let z = self.z_sorted();
        let mut heights = Array1::zeros(n_turbines);
        for ti in 0..n_turbines {
            heights[ti] = z[[0, ti, 0, 0]];
        }
        heights
    }
}
