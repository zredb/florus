# Python FLORIS vs Rust FLORIS (floris-rs) Examples 比较分析

## 项目背景

- **Python FLORIS**: NREL开发的原始版本 (v4.6.2)，基于Python的面向对象风电厂尾流建模软件
- **Rust FLORIS (floris-rs)**: Python版本的Rust重写项目，旨在提供更好的性能和内存安全

## 目录结构对比

### Python FLORIS Examples结构

```
floris/
├── examples/
│   ├── examples_getting_started/
│   ├── examples_wind_data/
│   ├── examples_turbine/
│   ├── examples_layout_optimization/
│   ├── examples_control_optimization/
│   └── ...
```

### Rust FLORIS Examples结构

```
florus/
├── examples/
│   ├── examples_control_optimization/
│   ├── examples_control_types/
│   ├── examples_heterogeneous/
│   ├── examples_layout_optimization/
│   ├── examples_turbine/
│   ├── examples_turbopark/
│   ├── examples_uncertain/
│   ├── examples_visualizations/
│   ├── examples_wind_data/
│   └── ...
```

## 核心差异对比

### 1. API设计哲学

| 特性 | Python FLORIS | Rust FLORIS |
|------|---------------|-------------|
| **编程范式** | 面向对象 | 结构化/函数式混合 |
| **类型系统** | 动态类型 | 静态类型 |
| **代码风格** | 简洁、声明式 | 显式、结构化 |
| **模型构建** | `FlorisModel("config.yaml")` | `Farm::new()`, `FlowField::new()` |
| **配置方式** | `fmodel.set(wind_speeds=..., wind_directions=...)` | 直接构建结构体 |
| **运行方式** | `fmodel.run()` | `model.run()` |

### 2. 代码示例对比

#### Python FLORIS - 基本使用示例

```python
from floris import FlorisModel
import numpy as np

# 简单启动
fmodel = FlorisModel("path/to/config.yaml")

# 设置风况
fmodel.set(
    wind_directions=np.array([270.0]),
    wind_speeds=[8.0],
    turbulence_intensities=np.array([0.06])
)

# 设置布局
fmodel.set(layout_x=[0, 500.0], layout_y=[0.0, 0.0])

# 运行模拟
fmodel.run()

# 获取结果
turbine_powers = fmodel.get_turbine_powers()
farm_power = fmodel.get_farm_power()
```

#### Rust FLORIS - 基本使用示例

```rust
use florus::core::{Farm, FlowField};
use florus::types::Array1;

fn main() -> anyhow::Result<()> {
    // 创建风电场布局
    let layout_x = Array1::from_vec(vec![0.0, 500.0]);
    let layout_y = Array1::from_vec(vec![0.0, 0.0]);
    let turbine_types = vec!["nrel_5MW".to_string(); 2];

    let farm = Farm::new(layout_x.clone(), layout_y.clone(), turbine_types)?;

    // 创建流场
    let wind_speeds = Array1::from_vec(vec![8.0]);
    let wind_directions = Array1::from_vec(vec![270.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06]);

    let flow_field = FlowField::new(
        wind_speeds.clone(),
        wind_directions.clone(),
        0.0,    // wind_veer
        0.14,   // wind_shear
        1.225,  // air_density
        turbulence_intensities.clone(),
        90.0,   // reference_wind_height
    )?;

    // 创建完整的FlorisModel
    let mut model = florus::FlorisModel {
        farm,
        flow_field,
        state: florus::core::State::new(),
        grid: None,
        solver_type: "turbine_grid".to_string(),
        model_manager: None,
    };

    // 初始化
    model.initialize_grid()?;
    model.initialize_flow_field()?;

    // 运行模拟
    model.run()?;

    // 获取结果
    let turbine_powers = model.get_turbine_powers();
    let farm_power = model.get_farm_power();

    Ok(())
}
```

## 主要差异点详解

### a) 初始化流程

- **Python**:
  - 一步到位：`fmodel = FlorisModel("config.yaml")`
  - 隐式初始化，API自动处理

- **Rust**:
  - 分步构建：`Farm::new()` → `FlowField::new()` → 手动构建`FlorisModel`
  - 显式初始化调用：`model.initialize_grid()`和`model.initialize_flow_field()`

### b) 错误处理机制

- **Python**:
  - 使用异常机制（try-except）
  - 运行时错误检测

- **Rust**:
  - 使用`Result<T, E>`类型
  - 编译时错误检测
  - 使用`?`操作符传播错误
  - 需要显式的错误类型声明（`anyhow::Result`）

### c) 类型安全性

- **Python**:
  - 动态类型，运行时类型检查
  - 灵活但可能在运行时出现类型错误

- **Rust**:
  - 静态类型，编译时类型检查
  - 更强的类型安全保障
  - 明确的类型声明（`Array1`, `Array2`等）

### d) 内存管理

- **Python**:
  - 自动垃圾回收
  - 不需要手动内存管理

- **Rust**:
  - 所有权系统
  - 编译时内存安全保证
  - 需要考虑借用和生命周期

### e) 代码组织

- **Python**:
  - 更少的代码行数
  - 更直观的API设计
  - 更少的样板代码

- **Rust**:
  - 更多的代码行数
  - 更结构化的代码
  - 显式的依赖和初始化

## 性能特点

| 特性 | Python FLORIS | Rust FLORIS |
|------|---------------|-------------|
| **执行速度** | 解释执行，较慢 | 编译执行，更快 |
| **内存使用** | 动态内存管理 | 静态内存分配 |
| **并发支持** | GIL限制 | 无GIL，真正的并行 |
| **优化潜力** | 有限 | 零成本抽象，优化空间大 |

## 开发体验

| 方面 | Python FLORIS | Rust FLORIS |
|------|---------------|-------------|
| **学习曲线** | 较低 | 较高 |
| **开发速度** | 快 | 较慢 |
| **调试难度** | 较低 | 较高 |
| **重构安全性** | 有限 | 高（类型系统保证） |
| **工具链** | 成熟（pip, conda） | Cargo包管理器 |

## 优缺点分析

### Python FLORIS 优点

1. ✅ 代码简洁，易于学习和使用
2. ✅ 丰富的科学计算生态系统（NumPy, Pandas, Matplotlib）
3. ✅ 快速原型开发
4. ✅ 广泛的社区支持和文档
5. ✅ 易于集成到现有Python项目中

### Python FLORIS 缺点

1. ❌ 运行时性能较低
2. ❌ 动态类型可能导致运行时错误
3. ❌ GIL限制真正的并行计算
4. ❌ 对于大型风电厂模拟，计算时间可能较长

### Rust FLORIS 优点

1. ✅ 高性能，编译时优化
2. ✅ 内存安全，无垃圾回收开销
3. ✅ 真正的并行计算支持
4. ✅ 强类型系统，编译时错误检测
5. ✅ 适合高性能计算和大规模模拟

### Rust FLORIS 缺点

1. ❌ 学习曲线陡峭
2. ❌ 开发周期较长
3. ❌ 样板代码较多
4. ❌ 生态系统不如Python成熟
5. ❌ 需要更多的显式错误处理

## 功能覆盖对比

| 功能模块 | Python FLORIS | Rust FLORIS |
|---------|---------------|-------------|
| **基础尾流模型** | ✅ 完整 | ✅ 已实现 |
| **偏航优化** | ✅ 完整 | ✅ 已实现 |
| **布局优化** | ✅ 完整 | ✅ 已实现 |
| **功率曲线定制** | ✅ 完整 | ⏳ 进行中 |
| **异构流入** | ✅ 完整 | ⏳ 进行中 |
| **不确定性分析** | ✅ 完整 | ⏳ 进行中 |
| **可视化工具** | ✅ 完整 | ⏳ 待实现 |
| **并行计算** | ✅ 完整 | ⏳ 待实现 |

## 迁移建议

对于考虑从Python FLORIS迁移到Rust FLORIS的用户：

### 适合迁移的场景

1. 大规模风电厂模拟
2. 需要高性能计算
3. 对计算时间敏感的应用
4. 长期运行的优化算法

### 不适合迁移的场景

1. 快速原型开发
2. 小规模模拟
3. 频繁修改模型参数
4. 需要丰富的可视化功能

## 总结

Python FLORIS和Rust FLORIS代表了两种不同的软件设计哲学：

- **Python FLORIS** 强调开发效率和易用性，适合快速原型和研究工作
- **Rust FLORIS** 强调性能和安全性，适合生产环境和大规模计算

两个项目都在积极开发中，Rust版本正在逐步完善功能，预计未来将提供更好的性能和更广泛的特性支持。

## 代码示例位置

- **Python FLORIS Examples**: `floris/examples/`
- **Rust FLORIS Examples**: `florus/examples/`

## 参考资源

- Python FLORIS: https://github.com/NREL/floris
- Rust FLORIS: https://github.com/floris-rs/floris (假设)
- FLORIS Documentation: https://nrel.github.io/floris/
