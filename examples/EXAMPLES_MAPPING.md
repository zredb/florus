# FLORUS 示例文件映射表

本文档记录了 Rust FLORUS 示例文件与 Python FLORIS v4.6.4 示例文件的对应关系。

## 📁 目录结构说明

所有示例文件现在都位于 `examples/` 主目录下，采用扁平化结构以兼容 Cargo 的示例发现机制。

**命名规则**: `{category}_{number}_{description}.rs`

---

## 📊 基础示例 (001-010)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `001_opening_floris_computing_power.rs` | `001_opening_floris_computing_power.py` | ✅ | FLORIS 入门和功率计算 |
| `002_visualization.rs` | `002_visualizations.py` | ✅ | 基础可视化 |
| `003_wind_data_objects.rs` | `003_wind_data_objects.py` | ✅ | WindData 对象 |
| `004_set.rs` | `004_set.py` | ✅ | 设置参数 |
| `005_getting_power.rs` | `005_getting_power.py` | ✅ | 获取功率 |
| `006_get_farm_aep.rs` | `006_get_farm_aep.py` | ✅ | 计算农场 AEP |
| `006_get_farm_aep_simple.rs` | - | ✅ | 简化版 AEP 计算（Rust 特有） |
| `007_sweeping_variables.rs` | `007_sweeping_variables.py` | ✅ | 变量扫描 |
| `008_uncertain_models.rs` | `008_uncertain_models.py` | ✅ | 不确定性模型 |
| `009_parallel_models.rs` | `009_parallel_models.py` | ✅ | 并行模型 |
| `010_compare_farm_power_with_neighbor.rs` | `010_compare_farm_power_with_neighbor.py` | ✅ | 与邻近农场比较 |

---

## 🎨 可视化示例 (visualizations_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `visualizations_001_layout_visualizations.rs` | `examples_visualizations/001_layout_visualizations.py` | ✅ | 布局可视化（生成8个PNG） |
| `visualizations_002_visualize_y_cut_plane.rs` | `examples_visualizations/002_visualize_y_cut_plane.py` | ⚠️ | Y 平面可视化（占位符） |
| `visualizations_003_visualize_cross_plane.rs` | `examples_visualizations/003_visualize_cross_plane.py` | ⚠️ | 横截面可视化（占位符） |
| `visualizations_004_visualize_rotor_values.rs` | `examples_visualizations/004_visualize_rotor_values.py` | ✅ | 转子值可视化 |
| `visualizations_005_visualize_flow_by_sweeping_turbines.rs` | `examples_visualizations/005_visualize_flow_by_sweeping_turbines.py` | ⚠️ | 涡轮机扫描流场（占位符） |

---

## 🎯 控制优化示例 (control_optimization_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `control_optimization_001_opt_yaw_single_ws.rs` | `examples_control_optimization/001_opt_yaw_single_ws.py` | ✅ | 单风速偏航优化 |
| `control_optimization_002_opt_yaw_single_ws_uncertain.rs` | `examples_control_optimization/002_opt_yaw_single_ws_uncertain.py` | ⚠️ | 带不确定性的偏航优化（占位符） |
| `control_optimization_003_opt_yaw_multiple_ws.rs` | `examples_control_optimization/003_opt_yaw_multiple_ws.py` | ⚠️ | 多风速偏航优化（占位符） |
| `control_optimization_004_optimize_yaw_aep.rs` | `examples_control_optimization/004_optimize_yaw_aep.py` | ⚠️ | AEP 偏航优化（占位符） |
| `control_optimization_005_optimize_yaw_aep_parallel.rs` | `examples_control_optimization/005_optimize_yaw_aep_parallel.py` | ⚠️ | 并行 AEP 优化（占位符） |
| `control_optimization_006_compare_yaw_optimizers.rs` | `examples_control_optimization/006_compare_yaw_optimizers.py` | ⚠️ | 比较偏航优化器（占位符） |
| `control_optimization_007_optimize_yaw_with_neighbor_farms.rs` | `examples_control_optimization/007_optimize_yaw_with_neighbor_farms.py` | ⚠️ | 与邻近农场优化（占位符） |
| `control_optimization_008_optimize_yaw_with_disabled_turbines.rs` | `examples_control_optimization/008_optimize_yaw_with_disabled_turbines.py` | ⚠️ | 禁用涡轮机优化（占位符） |

---

## 🎮 控制类型示例 (control_types_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `control_types_001_control_types_overview.rs` | `examples_control_types/001_control_types_overview.py` | ⚠️ | 控制类型概览（占位符） |
| `control_types_002_yaw_addition.rs` | `examples_control_types/002_yaw_addition.py` | ⚠️ | 偏航附加（占位符） |
| `control_types_003_tilt_addition.rs` | `examples_control_types/003_tilt_addition.py` | ⚠️ | 倾斜附加（占位符） |
| `control_types_004_power_setpoints.rs` | `examples_control_types/004_power_setpoints.py` | ⚠️ | 功率设定点（占位符） |
| `control_types_005_multiple_control_types.rs` | `examples_control_types/005_multiple_control_types.py` | ⚠️ | 多种控制类型（占位符） |

---

## 🌊 流场提取示例 (get_flow_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `get_flow_001_get_flow.rs` | `examples_get_flow/001_extract_wind_speed_at_turbines.py` | ✅ | 提取涡轮机风速 |
| `get_flow_002_get_flow_at_turbines.rs` | `examples_get_flow/002_extract_wind_speed_at_points.py` | ⚠️ | 在点提取风速（占位符） |
| `get_flow_003_get_flow_on_grid.rs` | `examples_get_flow/003_extract_turbulence_intensity_at_points.py` | ⚠️ | 网格上提取流场（占位符） |
| `get_flow_004_get_flow_with_custom_planes.rs` | `examples_get_flow/004_plot_velocity_deficit_profiles.py` | ⚠️ | 自定义平面流场（占位符） |

---

## 🌪️ 异构图示例 (heterogeneous_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `heterogeneous_001_heterogeneous_inflow.rs` | `examples_heterogeneous/001_heterogeneous_inflow.py` | ⚠️ | 非均匀入流（占位符） |
| `heterogeneous_002_het_map_from_file.rs` | `examples_heterogeneous/002_het_map_from_file.py` | ⚠️ | 从文件加载异构图（占位符） |
| `heterogeneous_003_het_multi_turbine.rs` | `examples_heterogeneous/003_het_multi_turbine.py` | ⚠️ | 多涡轮机异构（占位符） |
| `heterogeneous_004_het_with_wind_rose.rs` | `examples_heterogeneous/004_het_with_wind_rose.py` | ⚠️ | 带风玫瑰图的异构（占位符） |

---

## 🎈 浮式风机示例 (floating_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `floating_001_floating_turbine_models.rs` | `examples_floating/001_floating_turbine_models.py` | ⚠️ | 浮式风机模型（占位符） |
| `floating_002_floating_vs_fixedbottom_farm.rs` | `examples_floating/002_floating_vs_fixedbottom_farm.py` | ⚠️ | 浮式与固定式对比（占位符） |
| `floating_003_tilt_driven_vertical_wake_deflection.rs` | `examples_floating/003_tilt_driven_vertical_wake_deflection.py` | ⚠️ | 倾斜驱动垂直尾流（占位符） |

---

## 📐 布局优化示例 (layout_optimization_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `layout_optimization_001_layout_optimization.rs` | `examples_layout_optimization/001_layout_optimization.py` | ⚠️ | 布局优化（占位符） |
| `layout_optimization_002_layout_optimization_with_wind_rose.rs` | `examples_layout_optimization/002_layout_optimization_with_wind_rose.py` | ⚠️ | 带风玫瑰图的布局优化（占位符） |
| `layout_optimization_003_layout_optimization_gridded.rs` | `examples_layout_optimization/003_layout_optimization_gridded.py` | ⚠️ | 网格化布局优化（占位符） |
| `layout_optimization_004_layout_optimization_random_search.rs` | `examples_layout_optimization/004_layout_optimization_random_search.py` | ⚠️ | 随机搜索布局优化（占位符） |
| `layout_optimization_005_layout_optimization_boundary_grid.rs` | `examples_layout_optimization/005_layout_optimization_boundary_grid.py` | ⚠️ | 边界网格布局优化（占位符） |

---

## ⚖️ 载荷优化示例 (load_optimization_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `load_optimization_001_load_optimization.rs` | `examples_load_optimization/001_load_optimization.py` | ⚠️ | 载荷优化（占位符） |
| `load_optimization_002_load_optimization_with_wind_rose.rs` | `examples_load_optimization/002_load_optimization_with_wind_rose.py` | ⚠️ | 带风玫瑰图的载荷优化（占位符） |

---

## 🔢 多维示例 (multidim_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `multidim_001_multi_dimensional_cp_ct.rs` | `examples_multidim/001_multi_dimensional_cp_ct.py` | ⚠️ | 多维 CP/CT 表（占位符） |
| `multidim_002_multi_dimensional_cp_ct_2Hs.rs` | `examples_multidim/002_multi_dimensional_cp_ct_2Hs.py` | ⚠️ | 双 Hs 多维 CP/CT（占位符） |
| `multidim_003_multi_dimensional_Tp_Hs.rs` | `examples_multidim/003_multi_dimensional_Tp_Hs.py` | ⚠️ | Tp-Hs 多维表（占位符） |

---

## 🔧 运行模型示例 (operation_models_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `operation_models_001_operation_model.rs` | `examples_operation_models/001_operation_model.py` | ⚠️ | 运行模型（占位符） |

---

## 🌀 经验高斯模型示例 (emgauss_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `emgauss_001_empirical_gauss_velocity_deficit_parameters.rs` | `examples_emgauss/001_empirical_gauss_velocity_deficit_parameters.py` | ⚠️ | 经验高斯速度亏损参数（占位符） |
| `emgauss_002_empirical_gauss_helix.rs` | `examples_emgauss/002_empirical_gauss_helix.py` | ⚠️ | 经验高斯螺旋（占位符） |

---

## 🎲 不确定性示例 (uncertain_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `uncertain_001_uncertain_floris_model.rs` | `examples_uncertain/001_uncertain_floris_model.py` | ⚠️ | 不确定性 FLORIS 模型（占位符） |
| `uncertain_002_approx_floris_model.rs` | `examples_uncertain/002_approx_floris_model.py` | ⚠️ | 近似 FLORIS 模型（占位符） |
| `uncertain_003_parallel_uncertain.rs` | `examples_uncertain/003_parallel_uncertain.py` | ⚠️ | 并行不确定性（占位符） |

---

## 🌬️ 风资源网格示例 (wind_resource_grid_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `wind_resource_grid_001_wind_resource_grid.rs` | `examples_wind_resource_grid/001_wind_resource_grid.py` | ⚠️ | 风资源网格（占位符） |
| `wind_resource_grid_002_wrg_from_file.rs` | `examples_wind_resource_grid/002_wrg_from_file.py` | ⚠️ | 从文件加载 WRG（占位符） |
| `wind_resource_grid_003_wrg_with_heterogeneous.rs` | `examples_wind_resource_grid/003_wrg_with_heterogeneous.py` | ⚠️ | 带异构的 WRG（占位符） |
| `wind_resource_grid_004_wrg_visualization.rs` | `examples_wind_resource_grid/004_wrg_visualization.py` | ⚠️ | WRG 可视化（占位符） |
| `wind_resource_grid_005_wrg_multi_site.rs` | `examples_wind_resource_grid/005_wrg_multi_site.py` | ⚠️ | 多站点 WRG（占位符） |

---

## 💨 风数据示例 (wind_data_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `wind_data_001_wind_data_comparisons.rs` | `examples_wind_data/001_wind_data_comparisons.py` | ⚠️ | 风数据比较（占位符） |
| `wind_data_002_generate_ti.rs` | `examples_wind_data/002_generate_ti.py` | ⚠️ | 生成 TI 表（占位符） |
| `wind_data_003_generate_value.rs` | `examples_wind_data/003_generate_value.py` | ⚠️ | 生成值表（占位符） |

---

## 🏭 涡轮机示例 (turbine_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `turbine_001_turbine_library.rs` | `examples_turbine/001_turbine_library.py` | ⚠️ | 涡轮机库（占位符） |
| `turbine_002_turbine_interaction.rs` | `examples_turbine/002_turbine_interaction.py` | ⚠️ | 涡轮机交互（占位符） |
| `turbine_003_custom_turbine.rs` | `examples_turbine/003_custom_turbine.py` | ⚠️ | 自定义涡轮机（占位符） |

---

## 🅿️ TurboPark 示例 (turbopark_*)

| Rust 文件 | Python 文件 | 状态 | 说明 |
|-----------|-------------|------|------|
| `turbopark_001_turbopark_model.rs` | `examples_turbopark/001_turbopark.py` | ⚠️ | TurboPark 模型（占位符） |

---

## 📊 统计摘要

### 按类别统计

| 类别 | 文件数 | 已完成 | 占位符 | 完成率 |
|------|--------|--------|--------|--------|
| 基础示例 | 11 | 11 | 0 | 100% |
| visualizations | 5 | 2 | 3 | 40% |
| control_optimization | 8 | 1 | 7 | 12.5% |
| get_flow | 4 | 1 | 3 | 25% |
| control_types | 5 | 0 | 5 | 0% |
| heterogeneous | 4 | 0 | 4 | 0% |
| floating | 3 | 0 | 3 | 0% |
| layout_optimization | 5 | 0 | 5 | 0% |
| load_optimization | 2 | 0 | 2 | 0% |
| multidim | 3 | 0 | 3 | 0% |
| operation_models | 1 | 0 | 1 | 0% |
| emgauss | 2 | 0 | 2 | 0% |
| uncertain | 3 | 0 | 3 | 0% |
| wind_resource_grid | 5 | 0 | 5 | 0% |
| wind_data | 3 | 0 | 3 | 0% |
| turbine | 3 | 0 | 3 | 0% |
| turbopark | 1 | 0 | 1 | 0% |
| **总计** | **69** | **15** | **54** | **21.7%** |

---

## 🚀 运行示例

### 查看所有可用示例
```bash
cargo build --examples
ls examples/*.rs
```

### 运行特定示例
```bash
# 基础示例
cargo run --release --example 001_opening_floris_computing_power

# 可视化示例
cargo run --release --example visualizations_001_layout_visualizations

# 控制优化示例
cargo run --release --example control_optimization_001_opt_yaw_single_ws

# 流场提取示例
cargo run --release --example get_flow_001_get_flow
```

### 编译所有示例
```bash
cargo build --release --examples
```

---

## 📝 命名约定说明

### 文件名格式
```
{category}_{number}_{description}.rs
```

**示例**:
- `control_optimization_001_opt_yaw_single_ws.rs`
  - Category: `control_optimization`
  - Number: `001`
  - Description: `opt_yaw_single_ws`

### 类别前缀对照表

| 前缀 | 对应 Python 文件夹 | 说明 |
|------|-------------------|------|
| (无前缀) | examples/ | 基础示例 |
| `visualizations_` | examples_visualizations/ | 可视化 |
| `control_optimization_` | examples_control_optimization/ | 控制优化 |
| `control_types_` | examples_control_types/ | 控制类型 |
| `get_flow_` | examples_get_flow/ | 流场提取 |
| `heterogeneous_` | examples_heterogeneous/ | 异构图 |
| `floating_` | examples_floating/ | 浮式风机 |
| `layout_optimization_` | examples_layout_optimization/ | 布局优化 |
| `load_optimization_` | examples_load_optimization/ | 载荷优化 |
| `multidim_` | examples_multidim/ | 多维 |
| `operation_models_` | examples_operation_models/ | 运行模型 |
| `emgauss_` | examples_emgauss/ | 经验高斯 |
| `uncertain_` | examples_uncertain/ | 不确定性 |
| `wind_resource_grid_` | examples_wind_resource_grid/ | 风资源网格 |
| `wind_data_` | examples_wind_data/ | 风数据 |
| `turbine_` | examples_turbine/ | 涡轮机 |
| `turbopark_` | examples_turbopark/ | TurboPark |

---

## ✅ 优势

1. **Cargo 兼容**: 所有示例都在主目录，Cargo 可以自动发现
2. **清晰分类**: 通过前缀快速识别示例类别
3. **保持对应**: 与 Python FLORIS 一一对应，便于参考
4. **易于搜索**: 可以使用通配符搜索特定类别的示例
5. **简洁命名**: 去掉了冗余的 "examples_" 前缀

---

**最后更新**: 2026-04-20  
**FLORUS 版本**: v0.1.0  
**对应 Python FLORIS**: v4.6.4  
**总示例数**: 69个
