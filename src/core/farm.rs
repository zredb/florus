use crate::core::turbines::turbine_type::TurbineType;
use crate::core::turbines::TurbineLibrary;
use crate::core::{AveragingMethod, State, Turbine, TurbineGrid};
use crate::types::{Array1, Array2, Float, NumericDict};
use crate::utilities::load_yaml;
use crate::Array4;
use ndarray::{Array, Array2 as NdArray2, Array3 as NdArray3};
use serde_yaml::Value;
use std::collections::HashMap;
use std::fmt;
use std::path::Path;

const POWER_SETPOINT_DEFAULT: Float = 5000.0; // 假设默认功率为5MW

#[derive(Clone)]
pub struct Farm {
    pub layout_x: Array1,
    pub layout_y: Array1,
    pub turbine_types: Vec<TurbineType>,

    pub yaw_angles: Array2,
    pub yaw_angles_sorted: Array2,
    pub tilt_angles: Array2,
    pub tilt_angles_sorted: Array2,
    pub power_setpoints: Array2,
    pub power_setpoints_sorted: Array2,
    pub awc_modes: NdArray2<String>,
    pub awc_modes_sorted: NdArray2<String>,
    pub awc_amplitudes: Array2,
    pub awc_amplitudes_sorted: Array2,
    pub awc_frequencies: Array2,
    pub awc_frequencies_sorted: Array2,
    pub hub_heights: Array1,
    pub hub_heights_sorted: Array2,
    pub turbine_map: Vec<Turbine>,
    pub turbine_type_map: NdArray2<String>,
    pub turbine_type_map_sorted: NdArray2<String>,
   

    pub rotor_diameters: Array1,
    pub rotor_diameters_sorted: Array2,
    pub tsrs: Array1,
    pub tsrs_sorted: Array2,
    pub ref_tilts: Array1,
    pub ref_tilts_sorted: Array2,
    pub correct_cp_ct_for_tilt: Vec<bool>,
    pub correct_cp_ct_for_tilt_sorted: Array2,
    pub state: State,
}

impl fmt::Debug for Farm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Farm")
            .field("layout_x", &self.layout_x)
            .field("layout_y", &self.layout_y)
            .field("turbine_type", &self.turbine_types)
            .field("yaw_angles", &self.yaw_angles)
            .field("yaw_angles_sorted", &self.yaw_angles_sorted)
            .field("tilt_angles", &self.tilt_angles)
            .field("tilt_angles_sorted", &self.tilt_angles_sorted)
            .field("power_setpoints", &self.power_setpoints)
            .field("power_setpoints_sorted", &self.power_setpoints_sorted)
            .field("hub_heights", &self.hub_heights)
            .field("hub_heights_sorted", &self.hub_heights_sorted)
            .field("turbine_map", &self.turbine_map)
            .field("rotor_diameters", &self.rotor_diameters)
            .field("tsrs", &self.tsrs)
            .field("ref_tilts", &self.ref_tilts)
            .field("state", &self.state)
            .finish()
    }
}

impl Farm {
    pub fn new(
        layout_x: Array1,
        layout_y: Array1,
        turbine_types: Vec<String>,
    ) -> crate::Result<Self> {
        let n_turbines = layout_x.len();

        if layout_x.len() != layout_y.len() {
            anyhow::bail!("layout_x and layout_y must have the same number of entries");
        }

        if turbine_types.len() != 1 && turbine_types.len() != n_turbines {
            anyhow::bail!(
                "turbine_type must have the same number of entries as layout_x/layout_y or have \
                a single turbine_type value"
            );
        }

        // 确保TurbineLibrary已经被初始化
        if TurbineLibrary::get_loaded_turbines().is_empty() {
            // 如果还没有初始化，这里应该预先加载一些默认类型
            // 或者您可以确保在创建Farm之前调用TurbineLibrary::initialize()之类的函数
        }

        // 获取TurbineType引用
        let mut tts = Vec::new();
        for t in &turbine_types {
            match TurbineLibrary::get_turbine(t) {
                Some(turbine_type) => tts.push(turbine_type),
                None => {
                    anyhow::bail!(
                        "Turbine type '{}' not found in turbine library. Available types: {:?}",
                        t,
                        TurbineLibrary::get_loaded_turbines()
                    );
                }
            }
        }

        // 如果只有一个类型，扩展到所有涡轮机
        let turbine_types = if tts.len() == 1 && n_turbines > 1 {
            vec![tts[0].clone(); n_turbines]
        } else {
            tts
        };

        let mut farm = Self {
            layout_x,
            layout_y,
            turbine_types: turbine_types,
            yaw_angles: Array2::zeros((1, n_turbines)),
            yaw_angles_sorted: Array2::zeros((1, n_turbines)),
            tilt_angles: Array2::zeros((1, n_turbines)),
            tilt_angles_sorted: Array2::zeros((1, n_turbines)),
            power_setpoints: Array2::zeros((1, n_turbines)),
            power_setpoints_sorted: Array2::zeros((1, n_turbines)),
            awc_modes: NdArray2::from_elem((1, n_turbines), "baseline".to_string()),
            awc_modes_sorted: NdArray2::from_elem((1, n_turbines), "baseline".to_string()),
            awc_amplitudes: Array2::zeros((1, n_turbines)),
            awc_amplitudes_sorted: Array2::zeros((1, n_turbines)),
            awc_frequencies: Array2::zeros((1, n_turbines)),
            awc_frequencies_sorted: Array2::zeros((1, n_turbines)),
            hub_heights: Array1::zeros(n_turbines),
            hub_heights_sorted: Array2::zeros((1, n_turbines)),
            turbine_map: vec![],
            turbine_type_map: NdArray2::from_elem((1, n_turbines), String::new()),
            turbine_type_map_sorted: NdArray2::from_elem((1, n_turbines), String::new()),

            rotor_diameters: Array1::zeros(n_turbines),
            rotor_diameters_sorted: Array2::zeros((1, n_turbines)),
            tsrs: Array1::zeros(n_turbines),
            tsrs_sorted: Array2::zeros((1, n_turbines)),
            ref_tilts: Array1::zeros(n_turbines),
            ref_tilts_sorted: Array2::zeros((1, n_turbines)),
            correct_cp_ct_for_tilt: vec![false; n_turbines],
            correct_cp_ct_for_tilt_sorted: Array2::zeros((1, n_turbines)),
            state: State::default(),
        };

        farm.map_turbine_types()?;

        Ok(farm)
    }

    fn map_turbine_types(&mut self) -> crate::Result<()> {
        // 构建涡轮机映射
        self.construct_hub_heights();
        self.construct_rotor_diameters();
        self.construct_turbine_tsrs();
        self.construct_turbine_ref_tilts();
        self.construct_turbine_correct_cp_ct_for_tilt();
        Ok(())
    }

    pub fn initialize(&mut self, _sorted_indices: &NdArray3<usize>) {
        // 根据排序索引对偏航角进行排序
        // 这里使用简单的实现，实际需要使用更复杂的索引操作
        self.yaw_angles_sorted = self.yaw_angles.clone();
        self.tilt_angles_sorted = self.tilt_angles.clone();
        self.power_setpoints_sorted = self.power_setpoints.clone();
        self.awc_modes_sorted = self.awc_modes.clone();
        self.awc_amplitudes_sorted = self.awc_amplitudes.clone();
        self.awc_frequencies_sorted = self.awc_frequencies.clone();

        self.state.initialized = true;
    }

    pub fn construct_hub_heights(&mut self) {
        self.hub_heights = Array1::from_vec(
            self.turbine_types
                .iter()
                .map(|turb| turb.hub_height)
                .collect(),
        );
    }

    pub fn construct_rotor_diameters(&mut self) {
        self.rotor_diameters = Array1::from_vec(
            self.turbine_types
                .iter()
                .map(|turb| turb.rotor_diameter)
                .collect(),
        );
    }

    pub fn construct_turbine_tsrs(&mut self) {
        self.tsrs = Array1::from_vec(
            self.turbine_types
                .iter()
                .map(|turb| turb.tsr.unwrap_or(8.0))
                .collect(),
        );
    }

    pub fn construct_turbine_ref_tilts(&mut self) {
        self.ref_tilts = Array1::from_vec(
            self.turbine_types
                .iter()
                .map(|turb| turb.power_thrust_table.ref_tilt.unwrap_or(5.0))
                .collect(),
        );
    }

    pub fn construct_turbine_correct_cp_ct_for_tilt(&mut self) {
        self.correct_cp_ct_for_tilt = self
            .turbine_types
            .iter()
            .map(|turb| turb.correct_cp_ct_for_tilt)
            .collect();
    }

    // ... 其他方法保持不变
    pub fn expand_farm_properties(&mut self, n_findex: usize, sorted_coord_indices: &Array2) {
        let n_turbines = self.n_turbines();

        // Helper function to broadcast array2 from (1, n) to (n_findex, n)
        let broadcast_array2 = |arr: &Array2| -> Array2 {
            if arr.shape()[0] == n_findex {
                arr.clone()
            } else if arr.shape()[0] == 1 {
                let mut expanded = Array2::zeros((n_findex, n_turbines));
                for fi in 0..n_findex {
                    for ti in 0..n_turbines {
                        expanded[[fi, ti]] = arr[[0, ti]];
                    }
                }
                expanded
            } else {
                Array2::zeros((n_findex, n_turbines))
            }
        };

        // Helper function to broadcast NdArray2<String> from (1, n) to (n_findex, n)
        let broadcast_string_array2 = |arr: &NdArray2<String>| -> NdArray2<String> {
            if arr.shape()[0] == n_findex {
                arr.clone()
            } else if arr.shape()[0] == 1 {
                let mut expanded = NdArray2::from_elem((n_findex, n_turbines), String::new());
                for fi in 0..n_findex {
                    for ti in 0..n_turbines {
                        expanded[[fi, ti]] = arr[[0, ti]].clone();
                    }
                }
                expanded
            } else {
                NdArray2::from_elem((n_findex, n_turbines), String::new())
            }
        };

        // Broadcast arrays from first findex to all findex
        let yaw_angles_expanded = broadcast_array2(&self.yaw_angles);
        let tilt_angles_expanded = broadcast_array2(&self.tilt_angles);
        let power_setpoints_expanded = broadcast_array2(&self.power_setpoints);
        let awc_modes_expanded = broadcast_string_array2(&self.awc_modes);
        let awc_amplitudes_expanded = broadcast_array2(&self.awc_amplitudes);
        let awc_frequencies_expanded = broadcast_array2(&self.awc_frequencies);
        let turbine_type_map_expanded = broadcast_string_array2(&self.turbine_type_map);

        // Helper function to sort array1 according to sorted_indices for each findex
        let sort_array1_for_findex = |arr: &Array1, fi: usize| -> Array1 {
            let mut sorted = Array1::zeros(n_turbines);
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                sorted[new_i] = arr[old_i];
            }
            sorted
        };
        let sort_bool_array1_for_findex = |arr: &Vec<bool>, fi: usize| -> Array1 {
            let mut sorted = Array1::zeros(n_turbines);
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                sorted[new_i] = if arr[old_i] { 1.0 } else { 0.0 }; // 将bool转换为f64
            }
            sorted
        };

        // Expand and sort hub_heights
        let _hub_heights_expanded = self
            .hub_heights
            .to_owned()
            .insert_axis(ndarray::Axis(0))
            .broadcast((n_findex, n_turbines))
            .unwrap()
            .to_owned();
        self.hub_heights_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            let sorted = sort_array1_for_findex(&self.hub_heights, fi);
            for ti in 0..n_turbines {
                self.hub_heights_sorted[[fi, ti]] = sorted[ti];
            }
        }

        // Expand and sort rotor_diameters
        self.rotor_diameters_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            let sorted = sort_array1_for_findex(&self.rotor_diameters, fi);
            for ti in 0..n_turbines {
                self.rotor_diameters_sorted[[fi, ti]] = sorted[ti];
            }
        }

        // Expand and sort tsrs
        self.tsrs_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            let sorted = sort_array1_for_findex(&self.tsrs, fi);
            for ti in 0..n_turbines {
                self.tsrs_sorted[[fi, ti]] = sorted[ti];
            }
        }

        // Expand and sort ref_tilts
        self.ref_tilts_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            let sorted = sort_array1_for_findex(&self.ref_tilts, fi);
            for ti in 0..n_turbines {
                self.ref_tilts_sorted[[fi, ti]] = sorted[ti];
            }
        }

        // Expand and sort correct_cp_ct_for_tilt
        self.correct_cp_ct_for_tilt_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            let sorted = sort_bool_array1_for_findex(&self.correct_cp_ct_for_tilt, fi);
            for ti in 0..n_turbines {
                self.correct_cp_ct_for_tilt_sorted[[fi, ti]] = sorted[ti];
            }
        }

        // Sort yaw_angles according to sorted_indices for each findex
        self.yaw_angles_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                self.yaw_angles_sorted[[fi, new_i]] = yaw_angles_expanded[[fi, old_i]];
            }
        }

        // Sort tilt_angles according to sorted_indices for each findex
        self.tilt_angles_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                self.tilt_angles_sorted[[fi, new_i]] = tilt_angles_expanded[[fi, old_i]];
            }
        }

        // Sort power_setpoints according to sorted_indices for each findex
        self.power_setpoints_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                self.power_setpoints_sorted[[fi, new_i]] = power_setpoints_expanded[[fi, old_i]];
            }
        }

        // Sort awc_modes according to sorted_indices for each findex
        self.awc_modes_sorted = NdArray2::from_elem((n_findex, n_turbines), "baseline".to_string());
        for fi in 0..n_findex {
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                self.awc_modes_sorted[[fi, new_i]] = awc_modes_expanded[[fi, old_i]].clone();
            }
        }

        // Sort awc_amplitudes according to sorted_indices for each findex
        self.awc_amplitudes_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                self.awc_amplitudes_sorted[[fi, new_i]] = awc_amplitudes_expanded[[fi, old_i]];
            }
        }

        // Sort awc_frequencies according to sorted_indices for each findex
        self.awc_frequencies_sorted = Array2::zeros((n_findex, n_turbines));
        for fi in 0..n_findex {
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                self.awc_frequencies_sorted[[fi, new_i]] = awc_frequencies_expanded[[fi, old_i]];
            }
        }

        // Sort turbine_type_map according to sorted_indices for each findex
        self.turbine_type_map_sorted = NdArray2::from_elem((n_findex, n_turbines), String::new());
        for fi in 0..n_findex {
            for new_i in 0..n_turbines {
                let old_i = sorted_coord_indices[[fi, new_i]] as usize;
                self.turbine_type_map_sorted[[fi, new_i]] =
                    turbine_type_map_expanded[[fi, old_i]].clone();
            }
        }
    }

    pub fn set_yaw_angles(&mut self, yaw_angles: Array2) {
        self.yaw_angles = yaw_angles;
    }

    pub fn set_yaw_angles_to_ref_yaw(&mut self, n_findex: usize) {
        let n_turbines = self.n_turbines();
        let yaw_angles = Array2::zeros((n_findex, n_turbines));
        self.set_yaw_angles(yaw_angles);
        self.yaw_angles_sorted = Array2::zeros((n_findex, n_turbines));
    }

    pub fn set_tilt_to_ref_tilt(&mut self, n_findex: usize) {
        let n_turbines = self.n_turbines();
        let tilt_angles = Array2::ones((n_findex, n_turbines))
            * self
                .ref_tilts
                .to_owned()
                .insert_axis(ndarray::Axis(0))
                .broadcast((n_findex, n_turbines))
                .unwrap()
                .to_owned();
        self.tilt_angles = tilt_angles;
        self.tilt_angles_sorted = Array2::ones((n_findex, n_turbines))
            * self
                .ref_tilts
                .to_owned()
                .insert_axis(ndarray::Axis(0))
                .broadcast((n_findex, n_turbines))
                .unwrap()
                .to_owned();
    }

    pub fn set_power_setpoints(&mut self, power_setpoints: Array2) {
        self.power_setpoints = power_setpoints;
    }

    pub fn set_power_setpoints_to_ref_power(&mut self, n_findex: usize) {
        let n_turbines = self.n_turbines();
        let power_setpoints = Array2::from_elem((n_findex, n_turbines), POWER_SETPOINT_DEFAULT);
        self.set_power_setpoints(power_setpoints);
        self.power_setpoints_sorted =
            Array2::from_elem((n_findex, n_turbines), POWER_SETPOINT_DEFAULT);
    }

    pub fn set_awc_modes(&mut self, awc_modes: NdArray2<String>) {
        self.awc_modes = awc_modes;
    }

    pub fn set_awc_modes_to_ref_mode(&mut self, n_findex: usize) {
        let n_turbines = self.n_turbines();
        let awc_modes = NdArray2::from_elem((n_findex, n_turbines), "baseline".to_string());
        self.set_awc_modes(awc_modes);
        self.awc_modes_sorted = NdArray2::from_elem((n_findex, n_turbines), "baseline".to_string());
    }

    pub fn set_awc_amplitudes(&mut self, awc_amplitudes: Array2) {
        self.awc_amplitudes = awc_amplitudes;
    }

    pub fn set_awc_amplitudes_to_ref_amp(&mut self, n_findex: usize) {
        let n_turbines = self.n_turbines();
        let awc_amplitudes = Array2::zeros((n_findex, n_turbines));
        self.set_awc_amplitudes(awc_amplitudes);
        self.awc_amplitudes_sorted = Array2::zeros((n_findex, n_turbines));
    }

    pub fn set_awc_frequencies(&mut self, awc_frequencies: Array2) {
        self.awc_frequencies = awc_frequencies;
    }

    pub fn set_awc_frequencies_to_ref_freq(&mut self, n_findex: usize) {
        let n_turbines = self.n_turbines();
        let awc_frequencies = Array2::zeros((n_findex, n_turbines));
        self.set_awc_frequencies(awc_frequencies);
        self.awc_frequencies_sorted = Array2::zeros((n_findex, n_turbines));
    }

    pub fn calculate_tilt_for_eff_velocities(
        &self,
        _rotor_effective_velocities: &Array2,
    ) -> Array2 {
        // 简单实现，实际需要根据有效速度计算倾斜角度
        Array2::zeros((1, self.n_turbines()))
    }

    pub fn finalize(&mut self, _unsorted_indices: &NdArray3<usize>) {
        // 恢复原始顺序
        // 这里简化实现，实际需要根据unsorted_indices重新排序
        self.state.converged = true;
    }

    pub fn coordinates(&self) -> Array2 {
        let n_turbines = self.n_turbines();
        let mut coords = Array2::zeros((n_turbines, 3));

        for i in 0..n_turbines {
            coords[[i, 0]] = self.layout_x[i];
            coords[[i, 1]] = self.layout_y[i];
            coords[[i, 2]] = if self.hub_heights.len() == n_turbines {
                self.hub_heights[i]
            } else {
                self.hub_heights[0]
            };
        }

        coords
    }

    pub fn n_turbines(&self) -> usize {
        self.layout_x.len()
    }

    /// Get yaw angles reference
    pub fn yaw_angles(&self) -> &Array2 {
        &self.yaw_angles
    }

    /// Get tilt angles reference
    pub fn tilt_angles(&self) -> &Array2 {
        &self.tilt_angles
    }

    /// Get hub heights reference
    pub fn hub_heights(&self) -> &Array1 {
        &self.hub_heights
    }

    /// Get rotor diameters reference
    pub fn rotor_diameters(&self) -> &Array1 {
        &self.rotor_diameters
    }

    /// Get turbine map reference
    pub fn turbine_map(&self) -> &[Turbine] {
        &self.turbine_map
    }

    /// Create grid from farm configuration
    pub fn create_grid(&self) -> crate::Result<TurbineGrid> {
        let coords = self.coordinates();
        TurbineGrid::new(
            coords,
            self.rotor_diameters.clone(),
            Array1::from_vec(vec![270.0]), // Default wind direction
            3,
        )
    }

    pub fn initialize_control_arrays(&mut self, n_findex: usize) {
        let n_turbines = self.n_turbines();

        self.yaw_angles = Array2::zeros((n_findex, n_turbines));
        self.yaw_angles_sorted = Array2::zeros((n_findex, n_turbines));

        self.tilt_angles = Array2::zeros((n_findex, n_turbines));
        self.tilt_angles_sorted = Array2::zeros((n_findex, n_turbines));

        self.power_setpoints = Array2::zeros((n_findex, n_turbines));
        self.power_setpoints_sorted = Array2::zeros((n_findex, n_turbines));

        self.awc_modes = NdArray2::from_elem((n_findex, n_turbines), "baseline".to_string());
        self.awc_modes_sorted = NdArray2::from_elem((n_findex, n_turbines), "baseline".to_string());

        self.awc_amplitudes = Array2::zeros((n_findex, n_turbines));
        self.awc_amplitudes_sorted = Array2::zeros((n_findex, n_turbines));

        self.awc_frequencies = Array2::zeros((n_findex, n_turbines));
        self.awc_frequencies_sorted = Array2::zeros((n_findex, n_turbines));

        self.turbine_type_map = NdArray2::from_elem((n_findex, n_turbines), String::new());
        self.turbine_type_map_sorted = NdArray2::from_elem((n_findex, n_turbines), String::new());
    }

    /// Set turbine layout
    ///
    /// Updates the turbine positions and reinitializes related properties.
    pub fn set_layout(&mut self, layout_x: &Array1, layout_y: &Array1) -> crate::Result<()> {
        if layout_x.len() != layout_y.len() {
            anyhow::bail!("layout_x and layout_y must have the same number of entries");
        }

        let n_turbines = layout_x.len();
        self.layout_x = layout_x.clone();
        self.layout_y = layout_y.clone();

        // 更新turbine_types以适应新的布局
        if self.turbine_types.len() == 1 {
            let single_type = self.turbine_types[0].clone();
            self.turbine_types = vec![single_type; n_turbines];
        } else if self.turbine_types.len() > n_turbines {
            self.turbine_types = self.turbine_types[..n_turbines].to_vec();
        } else if self.turbine_types.len() < n_turbines {
            let last_type = self.turbine_types.last().clone().unwrap();
            self.turbine_types = vec![last_type.clone(); n_turbines];
        }

        // Reconstruct derived properties
        self.turbine_map.clear();
        self.construct_turbine_map();
        self.construct_hub_heights();
        self.construct_rotor_diameters();
        self.construct_turbine_tsrs();
        self.construct_turbine_ref_tilts();
        self.construct_turbine_correct_cp_ct_for_tilt();

        Ok(())
    }

    // Calculate power for array of turbines
    pub fn power(
        velocities: &Array4,
        turbines: &[Turbine],
        air_density: Float,
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        average_method: AveragingMethod,
    ) -> crate::Result<Array2> {
        let shape = velocities.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];

        let mut power_output = Array::zeros((n_findex, n_turbines));

        for ti in 0..n_turbines {
            if ti < turbines.len() {
                let turbine_power = turbines[ti].calculate_power(
                    velocities,
                    air_density,
                    yaw_angles,
                    tilt_angles,
                    average_method,
                )?;

                for fi in 0..n_findex {
                    power_output[[fi, ti]] = turbine_power[[fi, 0]];
                }
            }
        }

        Ok(power_output)
    }

    /// Calculate thrust coefficient for array of turbines
    pub fn thrust_coefficient(
        velocities: &Array4,
        turbines: &[Turbine],
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        average_method: AveragingMethod,
    ) -> crate::Result<Array2> {
        let shape = velocities.shape();
        let n_findex = shape[0];
        let n_turbines = shape[1];

        let mut ct_output = Array::zeros((n_findex, n_turbines));

        for ti in 0..n_turbines {
            if ti < turbines.len() {
                let turbine_ct = turbines[ti].calculate_thrust_coefficient(
                    velocities,
                    yaw_angles,
                    tilt_angles,
                    average_method,
                )?;

                for fi in 0..n_findex {
                    ct_output[[fi, ti]] = turbine_ct[[fi, 0]];
                }
            }
        }

        Ok(ct_output)
    }

    /// Calculate axial induction from thrust coefficient
    pub fn axial_induction(
        velocities: &Array4,
        turbines: &[Turbine],
        yaw_angles: Option<&Array2>,
        tilt_angles: Option<&Array2>,
        average_method: AveragingMethod,
    ) -> crate::Result<Array2> {
        let ct = Self::thrust_coefficient(
            velocities,
            turbines,
            yaw_angles,
            tilt_angles,
            average_method,
        )?;

        let mut ai = Array::zeros(ct.dim());

        for ((i, j), &ct_val) in ct.indexed_iter() {
            if j < turbines.len() {
                ai[[i, j]] = turbines[j].calculate_axial_induction(ct_val);
            }
        }

        Ok(ai)
    }
}
