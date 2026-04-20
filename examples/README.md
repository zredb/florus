# FLORIS Rust 示例说明

本文档描述了 FLORIS Rust 实现 (florus) 中的示例程序。

## 已实现的示例

### 基础示例

- **001_opening_floris_computing_power.rs** - 打开 FLORIS 并计算功率
  - 演示基本的 FlorisModel 初始化和功率计算
  - 使用 GCH 尾流模型

- **004_set.rs** - 设置方法演示
  - 设置风况（风速、风向、湍流强度）
  - 设置风电场布局
  - 设置控制参数（偏航角）
  - 重置操作

- **005_getting_power.rs** - 获取功率
  - 演示如何获取涡轮机和风电场功率
  - TimeSeries 条件：扫描风向（250°-289°）
  - WindRose 条件：多风向×多风速组合
  - 计算尾流损失

- **006_get_farm_aep_simple.rs** - 计算风电场 AEP
  - 演示年发电量（AEP）计算
  - 均匀和非均匀频率情况
  - 尾流损失分析

### 高级示例

- **007_sweeping_variables.rs** - 变量扫描
  - 扫描风速（5.0-9.9 m/s）
  - 扫描风向（250°-289°）
  - 扫描湍流强度（0.03-0.19）
  - 扫描偏航角（-30° to 30°）
  - 展示如何在单次运行中处理多个工况

- **008_uncertain_models.rs** - 不确定性模型（简化版）
  - 手动模拟风向不确定性
  - 高斯加权平均计算
  - 与标称情况对比
  - 注：完整 UncertainFlorisModel 将在未来版本实现

- **009_parallel_models.rs** - 并行模型（概念演示）
  - 展示多模型实例的使用
  - 验证结果一致性
  - 提供并行化建议（Rayon、tokio）
  - 注：完整 ParFlorisModel 将在未来版本实现

- **010_compare_farm_power_with_neighbor.rs** - 与邻近风电场对比
  - 演示有/无邻近风机的功率对比
  - 8 涡轮机场景（4+4）
  - 分析不同风向下的功率差异
  - 注：完整 turbine_weights 功能将在未来版本实现

## 运行示例

```bash
# 运行单个示例
cargo run --release --example 001_opening_floris_computing_power
cargo run --release --example 004_set
cargo run --release --example 005_getting_power
cargo run --release --example 006_get_farm_aep_simple
cargo run --release --example 007_sweeping_variables
cargo run --release --example 008_uncertain_models
cargo run --release --example 009_parallel_models
cargo run --release --example 010_compare_farm_power_with_neighbor

# 编译所有示例
cargo build --release --examples
```

## 与 Python FLORIS 的对比

### 已实现的功能
✅ 基本 FlorisModel 功能  
✅ 风况设置  
✅ 布局设置  
✅ 功率计算  
✅ AEP 计算  
✅ 变量扫描  
✅ 尾流模型（Gauss, GCH）  

### 待实现的功能
⏳ UncertainFlorisModel（完整实现）  
⏳ ParFlorisModel（完整实现）  
⏳ turbine_weights 参数  
⏳ Power setpoints 和 derating  
⏳ Turbine disable 功能  

## 已知限制

1. **多次调用 set_wind_conditions**：✅ 已修复！现在可以安全地多次调用 `set_wind_conditions` 来切换不同的工况。

2. **并行化**：目前未实现完整的并行模型，但提供了概念演示和未来实现建议。

3. **不确定性建模**：当前仅提供了手动模拟方法，完整的 UncertainFlorisModel 类尚未实现。

## 数值一致性

所有示例的计算结果已与 Python FLORIS 4.6.4 版本进行对比验证，确保数值一致性（差异 < 2.5%）。

关键修复：
- Gauss 速度模型中涡轮机中心坐标计算（使用网格点平均值）
- 数组排序与未排序的正确使用
- 控制数组的正确初始化

## 下一步开发计划

1. 实现完整的 UncertainFlorisModel
2. 实现 ParFlorisModel 支持真正的并行计算
3. 添加 turbine_weights 参数支持
4. 实现 power setpoints 和 derating 功能
5. 添加更多可视化示例
