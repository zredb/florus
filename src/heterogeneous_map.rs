use crate::core::models::InterpMethod;
use crate::types::{Array1, Array2, Float};
use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HeterogeneousMap {
    pub x: Array1,
    pub y: Array1,
    pub speed_multipliers: Array2, // [n_directions, n_points] or [n_speeds, n_points] or [1, n_points]
    pub z: Option<Array1>,
    pub wind_directions: Option<Array1>,
    pub wind_speeds: Option<Array1>,
    pub interp_method: InterpMethod,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MultidimConditions {
    pub tp: Array1,
    pub hs: Option<Array1>,
}

#[derive(Clone, Debug)]
pub struct HeterogeneousInflowConfig {
    pub x: Array1,
    pub y: Array1,
    pub z: Option<Array1>,
    pub wind_speeds: Option<Array1>,
    pub wind_directions: Option<Array1>,
    pub speed_multipliers: Array2, // shape: (n_conditions, n_points)
}

impl HeterogeneousMap {
    pub fn new(
        x: Array1,
        y: Array1,
        speed_multipliers: Array2,
        z: Option<Array1>,
        wind_directions: Option<Array1>,
        wind_speeds: Option<Array1>,
        interp_method: InterpMethod,
    ) -> Result<Self> {
        let n_points = x.len();
        let n_speeds = speed_multipliers.shape()[0];

        if x.len() != n_points || y.len() != n_points {
            bail!("Length of x and y must match the number of points in speed_multipliers");
        }
        if let Some(z_vals) = &z {
            if z_vals.len() != n_points {
                bail!("Length of z must match the number of points in speed_multipliers");
            }
        }
        if let Some(directions) = &wind_directions {
            if directions.len() != n_speeds {
                bail!(
                    "Length of wind_directions must match the first dimension of speed_multipliers"
                );
            }
        }
        if let Some(speeds) = &wind_speeds {
            if speeds.len() != n_speeds {
                bail!("Length of wind_speeds must match the first dimension of speed_multipliers");
            }
        }
        if wind_directions.is_some() && wind_speeds.is_some() {
            if speed_multipliers.shape()[0] != 1 {
                bail!("If both wind_directions and wind_speeds are specified, speed_multipliers should have length 1 in 0th dimension");
            }
        }

        Ok(Self {
            x,
            y,
            speed_multipliers,
            z,
            wind_directions,
            wind_speeds,
            interp_method,
        })
    }

    /// Get the heterogeneous inflow configuration for given wind directions and wind speeds.
    ///
    /// Returns a HeterogeneousInflowConfig containing x, y, and speed_multipliers
    /// for the given wind conditions.
    pub fn get_heterogeneous_inflow_config(
        &self,
        wind_directions: Array1,
        wind_speeds: Array1,
    ) -> Result<HeterogeneousInflowConfig> {
        if wind_directions.len() != wind_speeds.len() {
            bail!("wind_directions and wind_speeds must be the same length");
        }

        let n_conditions = wind_directions.len();
        let n_points = self.x.len();

        // Initialize output speed multipliers [n_conditions, n_points]
        let mut speed_multipliers_by_findex = Array2::zeros((n_conditions, n_points));

        // Select for wind direction first
        if let Some(ref config_directions) = self.wind_directions {
            // Calculate angle differences between requested and configured directions
            let angle_diffs: Array2 =
                Array2::from_shape_fn((n_conditions, config_directions.len()), |(i, j)| {
                    let diff = (wind_directions[i] - config_directions[j]).abs();
                    diff.min(360.0 - diff)
                });

            // If wind_speeds is None, return value by wind direction only
            if self.wind_speeds.is_none() {
                for i in 0..n_conditions {
                    // Find the index of minimum angle difference
                    let min_angle_idx = angle_diffs
                        .row(i)
                        .iter()
                        .enumerate()
                        .min_by(|(_, &a), (_, &b)| a.partial_cmp(&b).unwrap())
                        .map(|(idx, _)| idx)
                        .unwrap();

                    // Check if angle difference is within tolerance
                    if angle_diffs[[i, min_angle_idx]] > 1e-5 {
                        bail!("Provided wind_directions do not match those in heterogeneous map");
                    }

                    // Copy the speed multipliers for the closest direction
                    for j in 0..n_points {
                        speed_multipliers_by_findex[[i, j]] =
                            self.speed_multipliers[[min_angle_idx, j]];
                    }
                }
            } else {
                // Both wind_directions and wind_speeds are defined
                let config_speeds = self.wind_speeds.as_ref().unwrap();

                for i in 0..n_conditions {
                    // Find all indices in angle_diffs[i] that have minimum value
                    let row = angle_diffs.row(i);
                    let min_angle: Float =
                        row.iter().copied().fold(Float::INFINITY, |a, b| a.min(b));
                    let closest_wd_indices: Vec<usize> = row
                        .iter()
                        .enumerate()
                        .filter(|(_, &val)| (val - min_angle).abs() < 1e-9)
                        .map(|(idx, _)| idx)
                        .collect();

                    // Calculate speed differences for the closest wind direction indices
                    let speed_diffs: Vec<Float> = closest_wd_indices
                        .iter()
                        .map(|&wd_idx| {
                            let config_speed = config_speeds[wd_idx];
                            (wind_speeds[i] - config_speed).abs()
                        })
                        .collect();

                    // Find the index with minimum speed difference
                    let min_speed_idx = speed_diffs
                        .iter()
                        .enumerate()
                        .min_by(|(_, &a), (_, &b)| a.partial_cmp(&b).unwrap())
                        .map(|(idx, _)| idx)
                        .unwrap();

                    let closest_config_idx = closest_wd_indices[min_speed_idx];

                    // Check combined difference is within tolerance
                    let combined_diff = (angle_diffs[[i, closest_config_idx]].powi(2)
                        + speed_diffs[min_speed_idx].powi(2))
                    .sqrt();
                    if combined_diff > 1e-5 {
                        bail!("Provided wind_directions and wind_speeds do not match those in heterogeneous map");
                    }

                    // Copy the speed multipliers
                    for j in 0..n_points {
                        speed_multipliers_by_findex[[i, j]] =
                            self.speed_multipliers[[closest_config_idx, j]];
                    }
                }
            }
        } else if let Some(ref config_speeds) = self.wind_speeds {
            // Wind speeds are defined without wind direction
            let speed_diffs: Array2 =
                Array2::from_shape_fn((n_conditions, config_speeds.len()), |(i, j)| {
                    (wind_speeds[i] - config_speeds[j]).abs()
                });

            for i in 0..n_conditions {
                let min_speed_idx = speed_diffs
                    .row(i)
                    .iter()
                    .enumerate()
                    .min_by(|(_, &a), (_, &b)| a.partial_cmp(&b).unwrap())
                    .map(|(idx, _)| idx)
                    .unwrap();

                for j in 0..n_points {
                    speed_multipliers_by_findex[[i, j]] =
                        self.speed_multipliers[[min_speed_idx, j]];
                }
            }
        } else {
            // Both wind_directions and wind_speeds are None
            // Repeat the single row until length of wind_directions
            if self.speed_multipliers.shape()[0] != 1 {
                bail!("If both wind_directions and wind_speeds are None, speed_multipliers should have length 1 in 0th dimension");
            }

            for i in 0..n_conditions {
                for j in 0..n_points {
                    speed_multipliers_by_findex[[i, j]] = self.speed_multipliers[[0, j]];
                }
            }
        }

        Ok(HeterogeneousInflowConfig {
            x: self.x.clone(),
            y: self.y.clone(),
            z: self.z.clone(),
            wind_speeds: Some(wind_speeds),
            wind_directions: Some(wind_directions),
            speed_multipliers: speed_multipliers_by_findex,
        })
    }

    /// Return a HeterogeneousMap with only x and y coordinates and a constant z value.
    /// This selects from x, y and speed_multipliers where z is nearest to the given value.
    pub fn get_heterogeneous_map_2d(&self, z: Float) -> Result<HeterogeneousMap> {
        let z_values = self
            .z
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("No z values defined in HeterogeneousMap"))?;

        // Find the value in z that is closest to the given z value
        let closest_z_index = z_values
            .iter()
            .enumerate()
            .min_by(|(_, &a), (_, &b)| {
                (a - z)
                    .abs()
                    .partial_cmp(&((b - z).abs()))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(idx, _)| idx)
            .unwrap();

        // Get the value at that index
        let closest_z_value = z_values[closest_z_index];

        // Get all indices where z equals the closest value (within tolerance)
        let closest_z_indices: Vec<usize> = z_values
            .iter()
            .enumerate()
            .filter(|(_, &z_val)| (z_val - closest_z_value).abs() < 1e-9)
            .map(|(idx, _)| idx)
            .collect();

        // Get versions of x, y and speed_multipliers that include only the closest z values
        let x: Array1 = Array1::from_iter(closest_z_indices.iter().map(|&idx| self.x[idx]));
        let y: Array1 = Array1::from_iter(closest_z_indices.iter().map(|&idx| self.y[idx]));
        let speed_multipliers: Array2 = Array2::from_shape_fn(
            (self.speed_multipliers.shape()[0], closest_z_indices.len()),
            |(i, j)| self.speed_multipliers[[i, closest_z_indices[j]]],
        );

        Ok(HeterogeneousMap {
            x,
            y,
            speed_multipliers,
            z: None,
            wind_directions: self.wind_directions.clone(),
            wind_speeds: self.wind_speeds.clone(),
            interp_method: self.interp_method,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heterogeneous_map_new() {
        let x = Array1::from_vec(vec![0.0, 100.0, 200.0]);
        let y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let speed_multipliers = Array2::from_shape_vec((1, 3), vec![1.0, 1.1, 1.2]).unwrap();

        let het_map = HeterogeneousMap::new(
            x,
            y,
            speed_multipliers,
            None,
            None,
            None,
            InterpMethod::Linear,
        )
        .unwrap();

        assert_eq!(het_map.x.len(), 3);
        assert_eq!(het_map.speed_multipliers.shape(), &[1, 3]);
    }

    #[test]
    fn test_get_heterogeneous_inflow_config_no_direction_speed() {
        let x = Array1::from_vec(vec![0.0, 100.0, 200.0]);
        let y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let speed_multipliers = Array2::from_shape_vec((1, 3), vec![1.0, 1.1, 1.2]).unwrap();

        let het_map = HeterogeneousMap::new(
            x,
            y,
            speed_multipliers,
            None,
            None,
            None,
            InterpMethod::Linear,
        )
        .unwrap();

        let wind_directions = Array1::from_vec(vec![270.0, 280.0]);
        let wind_speeds = Array1::from_vec(vec![8.0, 10.0]);

        let config = het_map
            .get_heterogeneous_inflow_config(wind_directions, wind_speeds)
            .unwrap();

        assert_eq!(config.speed_multipliers.shape(), &[2, 3]);
        // Should have repeated the single row
        assert!((config.speed_multipliers[[0, 0]] - 1.0).abs() < 1e-10);
        assert!((config.speed_multipliers[[1, 1]] - 1.1).abs() < 1e-10);
    }

    #[test]
    fn test_get_heterogeneous_inflow_config_with_wind_directions() {
        let x = Array1::from_vec(vec![0.0, 100.0, 200.0]);
        let y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let speed_multipliers = Array2::from_shape_vec(
            (3, 3),
            vec![
                1.0, 1.1, 1.2, // for 270 deg
                1.05, 1.15, 1.25, // for 280 deg
                1.1, 1.2, 1.3, // for 290 deg
            ],
        )
        .unwrap();
        let wind_directions = Some(Array1::from_vec(vec![270.0, 280.0, 290.0]));

        let het_map = HeterogeneousMap::new(
            x,
            y,
            speed_multipliers,
            None,
            wind_directions,
            None,
            InterpMethod::Linear,
        )
        .unwrap();

        let wind_directions = Array1::from_vec(vec![270.0, 280.0]);
        let wind_speeds = Array1::from_vec(vec![8.0, 10.0]);

        let config = het_map
            .get_heterogeneous_inflow_config(wind_directions, wind_speeds)
            .unwrap();

        assert_eq!(config.speed_multipliers.shape(), &[2, 3]);
        // First row should match first row of speed_multipliers
        assert!((config.speed_multipliers[[0, 0]] - 1.0).abs() < 1e-10);
        // Second row should match second row of speed_multipliers
        assert!((config.speed_multipliers[[1, 0]] - 1.05).abs() < 1e-10);
    }

    #[test]
    fn test_get_heterogeneous_map_2d() {
        let x = Array1::from_vec(vec![0.0, 100.0, 0.0, 100.0]);
        let y = Array1::from_vec(vec![0.0, 0.0, 100.0, 100.0]);
        let z = Some(Array1::from_vec(vec![90.0, 90.0, 90.0, 90.0]));
        let speed_multipliers =
            Array2::from_shape_vec((2, 4), vec![1.0, 1.1, 1.2, 1.3, 1.05, 1.15, 1.25, 1.35])
                .unwrap();
        let wind_directions = Some(Array1::from_vec(vec![270.0, 280.0]));

        let het_map = HeterogeneousMap::new(
            x.clone(),
            y.clone(),
            speed_multipliers,
            z,
            wind_directions,
            None,
            InterpMethod::Linear,
        )
        .unwrap();

        let het_map_2d = het_map.get_heterogeneous_map_2d(90.0).unwrap();

        assert!(het_map_2d.z.is_none());
        assert_eq!(het_map_2d.x.len(), 4);
        assert_eq!(het_map_2d.y.len(), 4);
        assert_eq!(het_map_2d.speed_multipliers.shape(), &[2, 4]);
    }

    #[test]
    fn test_get_heterogeneous_map_2d_no_z() {
        let x = Array1::from_vec(vec![0.0, 100.0]);
        let y = Array1::from_vec(vec![0.0, 0.0]);
        let speed_multipliers = Array2::from_shape_vec((1, 2), vec![1.0, 1.1]).unwrap();

        let het_map = HeterogeneousMap::new(
            x,
            y,
            speed_multipliers,
            None,
            None,
            None,
            InterpMethod::Linear,
        )
        .unwrap();

        let result = het_map.get_heterogeneous_map_2d(90.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_heterogeneous_inflow_config_length_mismatch() {
        let x = Array1::from_vec(vec![0.0, 100.0]);
        let y = Array1::from_vec(vec![0.0, 0.0]);
        let speed_multipliers = Array2::from_shape_vec((1, 2), vec![1.0, 1.1]).unwrap();

        let het_map = HeterogeneousMap::new(
            x,
            y,
            speed_multipliers,
            None,
            None,
            None,
            InterpMethod::Linear,
        )
        .unwrap();

        let wind_directions = Array1::from_vec(vec![270.0, 280.0]);
        let wind_speeds = Array1::from_vec(vec![8.0]); // Different length

        let result = het_map.get_heterogeneous_inflow_config(wind_directions, wind_speeds);
        assert!(result.is_err());
    }
}
