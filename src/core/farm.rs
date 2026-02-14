use crate::core::{State, TurbineGrid, Turbine};
use crate::types::{Array1, Array2, Float, NumericDict};
use crate::utilities::load_yaml;
use serde_yaml::Value;
use ndarray::{Array2 as NdArray2, Array3 as NdArray3};
use std::collections::HashMap;
use std::fmt;
use std::path::Path;

const POWER_SETPOINT_DEFAULT: Float = 5000.0; // 假设默认功率为5MW

#[derive(Clone)]
pub struct Farm {
    pub layout_x: Array1,
    pub layout_y: Array1,
    pub turbine_type: Vec<String>,
    pub turbine_library_path: std::path::PathBuf,
    pub turbine_definitions: Vec<NumericDict>,
    pub turbine_thrust_coefficient_functions: HashMap<String, String>,
    pub turbine_axial_induction_functions: HashMap<String, String>,
    pub turbine_tilt_interps: HashMap<String, String>,
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
    pub turbine_power_functions: HashMap<String, String>,
    pub turbine_power_thrust_tables: HashMap<String, NumericDict>,
    pub rotor_diameters: Array1,
    pub rotor_diameters_sorted: Array2,
    pub tsrs: Array1,
    pub tsrs_sorted: Array2,
    pub ref_tilts: Array1,
    pub ref_tilts_sorted: Array2,
    pub correct_cp_ct_for_tilt: Array1,
    pub correct_cp_ct_for_tilt_sorted: Array2,
    pub state: State,
    // 私有属性
    _turbine_types: Vec<String>,
    _turbine_definition_cache: HashMap<String, NumericDict>,
}

impl fmt::Debug for Farm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Farm")
            .field("layout_x", &self.layout_x)
            .field("layout_y", &self.layout_y)
            .field("turbine_type", &self.turbine_type)
            .field("turbine_library_path", &self.turbine_library_path)
            .field("turbine_definitions", &self.turbine_definitions)
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
            .field("_turbine_types", &self._turbine_types)
            .finish()
    }
}

impl Farm {
    pub fn new(
        layout_x: Array1,
        layout_y: Array1,
        turbine_type: Vec<String>,
    ) -> crate::Result<Self> {
        let n_turbines = layout_x.len();

        if layout_x.len() != layout_y.len() {
            anyhow::bail!("layout_x and layout_y must have the same number of entries");
        }

        if turbine_type.len() != 1 && turbine_type.len() != n_turbines {
            anyhow::bail!(
                "turbine_type must have the same number of entries as layout_x/layout_y or have \
                a single turbine_type value"
            );
        }

        let mut farm = Self {
            layout_x,
            layout_y,
            turbine_type,
            turbine_library_path: Path::new("./turbine_library").to_path_buf(),
            turbine_definitions: vec![],
            turbine_thrust_coefficient_functions: HashMap::new(),
            turbine_axial_induction_functions: HashMap::new(),
            turbine_tilt_interps: HashMap::new(),
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
            turbine_power_functions: HashMap::new(),
            turbine_power_thrust_tables: HashMap::new(),
            rotor_diameters: Array1::zeros(n_turbines),
            rotor_diameters_sorted: Array2::zeros((1, n_turbines)),
            tsrs: Array1::zeros(n_turbines),
            tsrs_sorted: Array2::zeros((1, n_turbines)),
            ref_tilts: Array1::zeros(n_turbines),
            ref_tilts_sorted: Array2::zeros((1, n_turbines)),
            correct_cp_ct_for_tilt: Array1::zeros(n_turbines),
            correct_cp_ct_for_tilt_sorted: Array2::zeros((1, n_turbines)),
            state: State::default(), // 使用 State::new() 替代 State::UNINITIALIZED
            _turbine_types: vec![],
            _turbine_definition_cache: HashMap::new(),
        };

        farm.initialize_turbine_cache()?;
        farm.map_turbine_types()?;

        Ok(farm)
    }

    fn initialize_turbine_cache(&mut self) -> crate::Result<()> {
        // 检查 turbine_type 是否为文件名或预定义类型
        for t in &self.turbine_type {
            if self._turbine_definition_cache.contains_key(t) {
                continue; // 如果已经加载，跳过
            }

            // 尝试从文件加载
            let internal_fn = Path::new("turbine_library").join(&format!("{}.yaml", t));
            let external_fn = self.turbine_library_path.join(&format!("{}.yaml", t));

            let yaml_path = if internal_fn.exists() {
                internal_fn
            } else if external_fn.exists() {
                external_fn
            } else {
                anyhow::bail!("The turbine type: {} does not exist in either the internal or external turbine library.", t);
            };

            let value: Value = load_yaml(yaml_path)?;
            let turbine_def: NumericDict = serde_yaml::from_value(value)
                .map_err(|e| anyhow::anyhow!("Failed to parse turbine definition: {}", e))?;
            self._turbine_definition_cache
                .insert(t.clone(), turbine_def);
        }

        // 确保_turbine_types与输入类型一致
        self._turbine_types = self.turbine_type.clone();

        // 如果只有一个涡轮机定义，扩展到N个涡轮机
        if self._turbine_types.len() == 1 {
            self._turbine_types = vec![self._turbine_types[0].clone(); self.n_turbines()];
        }

        Ok(())
    }

    fn map_turbine_types(&mut self) -> crate::Result<()> {
        self.turbine_definitions = self
            ._turbine_types
            .iter()
            .map(|t| {
                self._turbine_definition_cache
                    .get(t)
                    .cloned()
                    .unwrap_or_else(|| panic!("Turbine definition not found for type: {}", t))
            })
            .collect();

        // 构建涡轮机映射
        self.construct_turbine_map();
        self.construct_hub_heights();
        self.construct_rotor_diameters();
        self.construct_turbine_tsrs();
        self.construct_turbine_ref_tilts();
        self.construct_turbine_correct_cp_ct_for_tilt();
        self.construct_turbine_thrust_coefficient_functions();
        self.construct_turbine_axial_induction_functions();
        self.construct_turbine_tilt_interps();
        self.construct_turbine_power_functions();
        self.construct_turbine_power_thrust_tables();

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

        self.state.initialized = true; // 使用 state.initialized = true 替代 state = State::INITIALIZED
    }

    pub fn construct_hub_heights(&mut self) {
        self.hub_heights = Array1::from_vec(
            self.turbine_definitions
                .iter()
                .map(|turb| turb.get_scalar("hub_height").unwrap_or(0.0))
                .collect(),
        );
    }

    pub fn construct_rotor_diameters(&mut self) {
        self.rotor_diameters = Array1::from_vec(
            self.turbine_definitions
                .iter()
                .map(|turb| turb.get_scalar("rotor_diameter").unwrap_or(0.0))
                .collect(),
        );
    }

    pub fn construct_turbine_tsrs(&mut self) {
        self.tsrs = Array1::from_vec(
            self.turbine_definitions
                .iter()
                .map(|turb| turb.get_scalar("TSR").unwrap_or(0.0))
                .collect(),
        );
    }

    pub fn construct_turbine_ref_tilts(&mut self) {
        self.ref_tilts = Array1::from_vec(
            self.turbine_definitions
                .iter()
                .map(|turb| {
                    turb.get_scalar("power_thrust_table.ref_tilt")
                        .unwrap_or(0.0)
                })
                .collect(),
        );
    }

    pub fn construct_turbine_correct_cp_ct_for_tilt(&mut self) {
        // 简化实现，因为Turbine结构体中没有correct_cp_ct_for_tilt字段
        self.correct_cp_ct_for_tilt = Array1::zeros(self.turbine_map.len());
    }

    pub fn construct_turbine_map(&mut self) {
        // 从缓存创建涡轮机映射
        use crate::core::turbine::turbine_type::{TurbineType, PowerThrustTable};

        let mut turbine_map_unique: HashMap<String, Turbine> = HashMap::new();

        // 这里需要创建实际的Turbine实例，因为我们没有from_dict实现
        for (k, v) in &self._turbine_definition_cache {
            // 使用从定义中提取的参数创建TurbineType，然后创建Turbine实例
            let operation_model_str = v
                .get_array("operation_model")
                .and_then(|arr| arr.first())
                .map(|s| s.to_string())
                .unwrap_or_else(|| "cosine-loss".to_string());

            // Try to get power_thrust_table data
            let power_thrust_table = PowerThrustTable {
                wind_speed: v.get_array("power_thrust_table.wind_speed")
                    .or_else(|| v.get_array("wind_speed"))
                    .unwrap_or(&[])
                    .to_vec(),
                power: v.get_array("power_thrust_table.power")
                    .or_else(|| v.get_array("power"))
                    .unwrap_or(&[])
                    .to_vec(),
                thrust_coefficient: v.get_array("power_thrust_table.thrust_coefficient")
                    .or_else(|| v.get_array("thrust_coefficient"))
                    .unwrap_or(&[])
                    .to_vec(),
                ref_air_density: v.get_scalar("power_thrust_table.ref_air_density"),
                ref_tilt: v.get_scalar("power_thrust_table.ref_tilt"),
                cosine_loss_exponent_yaw: v.get_scalar("power_thrust_table.cosine_loss_exponent_yaw"),
                cosine_loss_exponent_tilt: v.get_scalar("power_thrust_table.cosine_loss_exponent_tilt"),
            };

            let turbine_type = TurbineType {
                name: k.clone(),
                rotor_diameter: v.get_scalar("rotor_diameter").unwrap_or(126.0),
                hub_height: v.get_scalar("hub_height").unwrap_or(90.0),
                tsr: v.get_scalar("TSR").unwrap_or(8.0),
                operation_model: operation_model_str.clone(),
                ref_tilt: v.get_scalar("power_thrust_table.ref_tilt"),
                correct_cp_ct_for_tilt: v.get_bool("correct_cp_ct_for_tilt"),
                power_thrust_table: Some(power_thrust_table),
                // Legacy fields (may be empty if using nested format)
                power_curve_wind_speeds: v
                    .get_array("power_curve_wind_speeds")
                    .unwrap_or(&[])
                    .to_vec(),
                power_curve_powers: v.get_array("power_curve_powers").unwrap_or(&[]).to_vec(),
                thrust_coefficient_wind_speeds: v
                    .get_array("thrust_coefficient_wind_speeds")
                    .unwrap_or(&[])
                    .to_vec(),
                thrust_coefficient_values: v
                    .get_array("thrust_coefficient_values")
                    .unwrap_or(&[])
                    .to_vec(),
                controller_dependent_turbine_parameters: None,
            };

            let turbine = Turbine {
                turbine_type,
                operation_model: operation_model_str,
            };

            turbine_map_unique.insert(k.clone(), turbine);
        }

        self.turbine_map = self
            ._turbine_types
            .iter()
            .map(|k| turbine_map_unique.get(k).unwrap().clone())
            .collect();
    }

    pub fn construct_turbine_thrust_coefficient_functions(&mut self) {
        for (i, turbine) in self.turbine_map.iter().enumerate() {
            let turbine_type = &self._turbine_types[i];
            // 存储函数的标识符而不是函数本身
            self.turbine_thrust_coefficient_functions.insert(
                turbine_type.clone(),
                format!("thrust_coefficient_{}", turbine.turbine_type),
            );
        }
    }

    pub fn construct_turbine_axial_induction_functions(&mut self) {
        for (i, turbine) in self.turbine_map.iter().enumerate() {
            let turbine_type = &self._turbine_types[i];
            // 存储函数的标识符而不是函数本身
            self.turbine_axial_induction_functions.insert(
                turbine_type.clone(),
                format!("axial_induction_{}", turbine.turbine_type),
            );
        }
    }

    pub fn construct_turbine_tilt_interps(&mut self) {
        for (i, turbine) in self.turbine_map.iter().enumerate() {
            let turbine_type = &self._turbine_types[i];
            // 存储函数的标识符而不是函数本身
            self.turbine_tilt_interps.insert(
                turbine_type.clone(),
                format!("tilt_interp_{}", turbine.turbine_type),
            );
        }
    }

    pub fn construct_turbine_power_functions(&mut self) {
        for (i, turbine) in self.turbine_map.iter().enumerate() {
            let turbine_type = &self._turbine_types[i];
            // 存储函数的标识符而不是函数本身
            self.turbine_power_functions.insert(
                turbine_type.clone(),
                format!("power_function_{}", turbine.turbine_type),
            );
        }
    }

    pub fn construct_turbine_power_thrust_tables(&mut self) {
        for turbine_type in &self._turbine_types {
            if let Some(turbine_def) = self._turbine_definition_cache.get(turbine_type) {
                self.turbine_power_thrust_tables
                    .insert(turbine_type.clone(), turbine_def.clone());
            }
        }
    }

    pub fn expand_farm_properties(
        &mut self,
        n_findex: usize,
        sorted_coord_indices: &Array2,
    ) {
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
            let sorted = sort_array1_for_findex(&self.correct_cp_ct_for_tilt, fi);
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
                self.turbine_type_map_sorted[[fi, new_i]] = turbine_type_map_expanded[[fi, old_i]].clone();
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
        self.state.converged = true; // 使用 state.converged = true 替代 State::USED
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

        // Update _turbine_types to match the new number of turbines
        if self._turbine_types.len() == 1 {
            self._turbine_types = vec![self._turbine_types[0].clone(); n_turbines];
        } else if self._turbine_types.len() > n_turbines {
            self._turbine_types = self._turbine_types[..n_turbines].to_vec();
        } else if self._turbine_types.len() < n_turbines {
            let last_type = self._turbine_types.last().unwrap().clone();
            self._turbine_types.extend(vec![last_type; n_turbines - self._turbine_types.len()]);
        }

        // Reconstruct turbine type definitions for the new layout
        self.turbine_definitions = self
            ._turbine_types
            .iter()
            .map(|t| {
                self._turbine_definition_cache
                    .get(t)
                    .cloned()
                    .unwrap_or_else(|| panic!("Turbine definition not found for type: {}", t))
            })
            .collect();

        // Reconstruct derived properties
        self.construct_turbine_map();
        self.construct_hub_heights();
        self.construct_rotor_diameters();
        self.construct_turbine_tsrs();
        self.construct_turbine_ref_tilts();
        self.construct_turbine_correct_cp_ct_for_tilt();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Array1;

    #[test]
    fn test_farm_creation() {
        let layout_x = Array1::from_vec(vec![0.0, 630.0, 1260.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 3];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types);
        
        assert!(farm.is_ok());
        let farm = farm.unwrap();
        assert_eq!(farm.n_turbines(), 3);
    }

    #[test]
    fn test_farm_creation_single_turbine_type() {
        let layout_x = Array1::from_vec(vec![0.0, 630.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string()];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types);
        
        assert!(farm.is_ok());
        let farm = farm.unwrap();
        assert_eq!(farm.n_turbines(), 2);
    }

    #[test]
    fn test_farm_layout_mismatch() {
        let layout_x = Array1::from_vec(vec![0.0, 630.0]);
        let layout_y = Array1::from_vec(vec![0.0]);  // Mismatch: 1 element vs 2
        let turbine_types = vec!["nrel_5MW".to_string(); 2];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types);
        
        assert!(farm.is_err());
    }

    #[test]
    fn test_farm_turbine_type_mismatch() {
        let layout_x = Array1::from_vec(vec![0.0, 630.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(), "iea_10MW".to_string(), "extra".to_string()];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types);
        
        assert!(farm.is_err());
    }

    #[test]
    fn test_farm_coordinates() {
        let layout_x = Array1::from_vec(vec![0.0, 500.0]);
        let layout_y = Array1::from_vec(vec![100.0, -100.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 2];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();
        let coords = farm.coordinates();
        
        assert_eq!(coords.shape()[0], 2);
        assert_eq!(coords.shape()[1], 3);
        assert_eq!(coords[[0, 0]], 0.0);
        assert_eq!(coords[[0, 1]], 100.0);
        assert_eq!(coords[[1, 0]], 500.0);
        assert_eq!(coords[[1, 1]], -100.0);
    }

    #[test]
    fn test_farm_hub_heights() {
        let layout_x = Array1::from_vec(vec![0.0]);
        let layout_y = Array1::from_vec(vec![0.0]);
        let turbine_types = vec!["nrel_5MW".to_string()];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();
        let heights = farm.hub_heights();
        
        assert_eq!(heights.len(), 1);
        assert!(heights[0] > 0.0);  // Should be around 90m for nrel_5MW
    }

    #[test]
    fn test_farm_rotor_diameters() {
        let layout_x = Array1::from_vec(vec![0.0]);
        let layout_y = Array1::from_vec(vec![0.0]);
        let turbine_types = vec!["nrel_5MW".to_string()];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();
        let diameters = farm.rotor_diameters();
        
        assert_eq!(diameters.len(), 1);
        assert!(diameters[0] > 0.0);  // Should be around 126m for nrel_5MW
    }

    #[test]
    fn test_farm_yaw_angles() {
        let layout_x = Array1::from_vec(vec![0.0]);
        let layout_y = Array1::from_vec(vec![0.0]);
        let turbine_types = vec!["nrel_5MW".to_string()];
        
        let farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();
        let yaw = farm.yaw_angles();
        
        assert_eq!(yaw.shape()[0], 1);
        assert_eq!(yaw.shape()[1], 1);
    }

    #[test]
    fn test_farm_initialize_control_arrays() {
        let layout_x = Array1::from_vec(vec![0.0, 630.0]);
        let layout_y = Array1::from_vec(vec![0.0, 0.0]);
        let turbine_types = vec!["nrel_5MW".to_string(); 2];
        
        let mut farm = Farm::new(layout_x, layout_y, turbine_types).unwrap();
        farm.initialize_control_arrays(3);
        
        assert_eq!(farm.yaw_angles.shape()[0], 3);
        assert_eq!(farm.yaw_angles.shape()[1], 2);
        assert_eq!(farm.tilt_angles.shape()[0], 3);
        assert_eq!(farm.tilt_angles.shape()[1], 2);
        assert_eq!(farm.power_setpoints.shape()[0], 3);
        assert_eq!(farm.power_setpoints.shape()[1], 2);
    }
}
