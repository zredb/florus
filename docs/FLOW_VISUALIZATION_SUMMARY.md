# FLORUS Flow Visualization - Implementation Summary

## Overview

成功实现了FLORUS中完整的尾流可视化功能，包括：
1. 水平面可视化 (Horizontal Plane)
2. Y平面可视化 (Y-Plane) 
3. 交叉面可视化 (Cross-Plane)

这些功能完全对应Python FLORIS的可视化能力。

## 实现的功能

### 1. FlorisModel API方法

#### `calculate_horizontal_plane(height, x_resolution, y_resolution)`
- 计算指定高度的水平切面（x-y平面）
- 参数：
  - `height`: 切面高度（米）
  - `x_resolution`: x方向分辨率
  - `y_resolution`: y方向分辨率
- 返回：`CutPlane`对象

#### `calculate_y_plane(x_resolution, z_resolution, crossstream_dist)`
- 计算沿风向的垂直切面（x-z平面）
- 参数：
  - `x_resolution`: x方向分辨率
  - `z_resolution`: z方向分辨率
  - `crossstream_dist`: 横向距离（y坐标）
- 返回：`CutPlane`对象

#### `calculate_cross_plane(y_resolution, z_resolution, downstream_dist)`
- 计算垂直于风向的切面（y-z平面）
- 参数：
  - `y_resolution`: y方向分辨率
  - `z_resolution`: z方向分辨率
  - `downstream_dist`: 下游距离（x坐标）
- 返回：`CutPlane`对象

### 2. FlowFieldPlanarGrid增强

修复并增强了`FlowFieldPlanarGrid`以支持所有三种切面类型：

#### 添加的功能
- 添加`turbine_hub_heights`字段存储轮毂高度
- 修复z方向范围计算（使用hub heights而不是turbine coordinates）
- 重写`hub_heights()`方法返回存储的值

#### 支持的切面类型
- `"z"` - 水平面（x-y平面，固定z高度）
- `"y"` - Y平面（x-z平面，固定y位置）
- `"x"` - 交叉面（y-z平面，固定x位置）

### 3. Core Solver修复

修复了`sequential_solver`以支持FlowFieldPlanarGrid：
- 使用`grid_second_dim`而不是`n_turbines`初始化wake_field
- 解决了形状不匹配导致的运行时错误

## 示例程序

### 006_flow_visualization.rs
演示水平面可视化的三个场景：
1. 基本水平流场可视化
2. 偏航涡轮机流场可视化
3. 带转子和标签的流场可视化

### 007_y_and_cross_plane_visualization.rs
演示Y平面和交叉面可视化：
1. Y平面可视化（沿风向的垂直切面）
2. 交叉面可视化（垂直于风向的切面）

## 生成的输出文件

所有可视化图像保存在`examples/outputs/visualization/`目录：

1. `01_horizontal_flow.png` - 基本水平流场（200x200分辨率）
2. `02_yawed_flow.png` - 偏航涡轮机流场（第二台涡轮机偏航30°）
3. `03_flow_with_rotors_labels.png` - 带转子和标签的流场
4. `04_y_plane.png` - Y平面可视化（200x100分辨率）
5. `05_cross_plane.png` - 交叉面可视化（100x100分辨率）

## 修复的关键Bug

### 1. FlowFieldPlanarGrid hub_heights问题
**问题**: Grid::hub_heights()访问不存在的索引
**修复**: 
- 添加`turbine_hub_heights: Array1`字段到FlowFieldPlanarGrid
- 更新构造函数接收hub heights参数
- 重写Grid trait的`hub_heights()`方法

### 2. wake_field形状不匹配
**问题**: combination model的function中形状不兼容
**修复**: 使用`grid_second_dim`而不是`n_turbines`初始化wake_field
```rust
// 修复前
let mut wake_field: Array4 = Array::zeros((n_findex, n_turbines, grid_y_dim, grid_z_dim));

// 修复后
let grid_second_dim = shape[1];
let mut wake_field: Array4 = Array::zeros((n_findex, grid_second_dim, grid_y_dim, grid_z_dim));
```

### 3. z方向范围计算错误
**问题**: Y平面和交叉面的z范围使用turbine_coordinates的z值（都是0）
**修复**: 使用`turbine_hub_heights`计算z方向范围
```rust
// 修复前
let max_z = z.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
(0.001, 6.0 * max_z)

// 修复后
let min_height = turbine_hub_heights.iter().cloned().fold(Float::INFINITY, Float::min);
let max_height = turbine_hub_heights.iter().cloned().fold(Float::NEG_INFINITY, Float::max);
(min_height - 2.0 * max_diameter, max_height + 2.0 * max_diameter)
```

### 4. Farm初始化turbine_types问题
**问题**: 单个turbine type时只创建一个Turbine对象
**修复**: 检测单一类型并复制给所有涡轮机位置

### 5. 坐标旋转索引越界
**问题**: 
- `rotate_coordinates_rel_west`访问不存在的z坐标
- `reverse_rotate_coordinates_rel_west`索引越界
**修复**: 添加边界检查

## 使用方法

```rust
use florus::{FlorisModel, Result};
use florus::visualization::flow_visualization;
use ndarray::Array;

fn main() -> Result<()> {
    // 初始化模型
    let mut fmodel = FlorisModel::from_file("examples/inputs/gch.yaml")?;
    
    // 设置涡轮机布局
    let layout_x = Array::from_vec(vec![0.0, 500.0, 1000.0]);
    let layout_y = Array::from_vec(vec![0.0, 0.0, 0.0]);
    fmodel.set_layout(&layout_x, &layout_y)?;
    
    // 设置风况
    let wind_speeds = Array::from_vec(vec![8.0]);
    let wind_directions = Array::from_vec(vec![270.0]);
    fmodel.set(
        Some(wind_speeds),
        Some(wind_directions),
        None, None, None,
        None, None, None, None,
        None, None, None, None, None,
        None,
    )?;
    
    // 方法1: 水平面可视化
    let horizontal_plane = fmodel.calculate_horizontal_plane(90.0, 200, 200)?;
    flow_visualization::visualize_cut_plane(
        &horizontal_plane,
        "horizontal_plane.png",
        Some(1.0),
        Some(8.0),
        "coolwarm",
        false,
        "Horizontal Plane",
    )?;
    
    // 方法2: Y平面可视化
    let y_plane = fmodel.calculate_y_plane(200, 100, 0.0)?;
    flow_visualization::visualize_cut_plane(
        &y_plane,
        "y_plane.png",
        Some(3.0),
        Some(9.0),
        "coolwarm",
        false,
        "Y Cut Plane",
    )?;
    
    // 方法3: 交叉面可视化
    let cross_plane = fmodel.calculate_cross_plane(100, 100, 500.0)?;
    flow_visualization::visualize_cut_plane(
        &cross_plane,
        "cross_plane.png",
        Some(3.0),
        Some(9.0),
        "coolwarm",
        false,
        "Cross Plane",
    )?;
    
    Ok(())
}
```

## 技术细节

### CutPlane数据结构
```rust
pub struct CutPlane {
    pub data: CutPlaneData,
    pub normal_vector: String,
    pub resolution: (usize, usize),
}

pub struct CutPlaneData {
    pub x1: Array1,  // 第一坐标轴
    pub x2: Array1,  // 第二坐标轴
    pub x3: Array1,  // 第三坐标轴（法线方向）
    pub u: Array1,   // x方向速度
    pub v: Array1,   // y方向速度
    pub w: Array1,   // z方向速度
}
```

### 坐标系统
- 使用旋转坐标系进行尾流计算
- 转换回惯性坐标系用于可视化
- 正确处理不同数量的grid点和turbine点

### 可视化流程
1. 创建FlowFieldPlanarGrid
2. 设置到Core的grid字段
3. 重新初始化流场
4. 调用solve_for_viz求解
5. 提取指定切面的数据
6. 过滤到指定位置（容差0.1）
7. 创建CutPlane对象
8. 使用plotters渲染图像

## 与Python FLORIS的对比

| 功能 | Python FLORIS | Rust FLORUS | 状态 |
|------|--------------|-------------|------|
| calculate_horizontal_plane() | ✓ | ✓ | ✓ 完成 |
| calculate_y_plane() | ✓ | ✓ | ✓ 完成 |
| calculate_cross_plane() | ✓ | ✓ | ✓ 完成 |
| visualize_cut_plane() | ✓ | ✓ | ✓ 完成 |
| visualize_cut_plane_with_rotors() | ✓ | ✓ | ✓ 完成 |
| visualize_cut_plane_with_rotors_and_labels() | ✓ | ✓ | ✓ 完成 |
| 多turbine类型支持 | ✓ | ✓ | ✓ 完成 |
| yaw角度支持 | ✓ | ✓ | ✓ 完成 |
| CutPlane数据结构 | ✓ | ✓ | ✓ 完成 |

## 性能优化

1. 使用ndarray的高效多维数组操作
2. 优化的形状广播和扩展
3. release模式编译以获得最佳性能
4. 使用clone而不是不必要的重新分配

## 未来改进方向

1. 添加3D可视化支持
2. 添加动画生成能力
3. 优化大规模网格的内存使用
4. 添加更多颜色映射选项
5. 支持自定义后处理功能
6. 添加等高线标签支持
7. 添加涡轮机转子详细可视化
8. 支持交互式可视化（WebGL等）

## 结论

FLORUS现在具备完整的尾流可视化功能，完全对应Python FLORIS的能力。用户可以：
- 可视化水平流场分布
- 分析垂直切面的速度分布
- 研究偏航对尾流的影响
- 生成高质量的可视化图像

所有功能都经过测试验证，生成的图像清晰准确地展示了流场特征。
