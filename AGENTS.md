# FLORUS Agent Guidelines

## Build/Test Commands
- **Build**: `cargo build --release`
- **Single test**: `cargo test <test_name>` (e.g., `cargo test solver::tests::test_axial_induction_low_ct`)
- **Run all tests**: `cargo test`
- **Check formatting**: `cargo fmt -- --check`
- **Lint**: `cargo clippy -- -D warnings`

## Architecture

- **Crate type**: Library + binary. `lib.rs` exports the public API; `main.rs` is a thin CLI.
- **Python reference**: `floris/` contains FLORIS v4.6 (Python). This Rust crate ports it. When adding features, check the Python source first.
- **4D array convention**: All spatial fields use `[findex, turbine, y, z]` indexing (ndarray `Array4`).
- **Wake model dispatch**: Wake models are `Box<dyn Trait>` dispatched via `WakeModelManager`. Implementations live in `src/core/wake/{velocity,deflection,turbulence,combination}/`.
- **Turbine types**: Loaded from YAML in `turbine_library/`. Cp/Ct surfaces can come from `.npz` (Python) or `.csv` files.

## Build quirks

- `.cargo/config.toml` sets `rustflags = ["-A", "warnings"]` — **warnings are suppressed in `cargo build`**. Always use `cargo clippy -- -D warnings` for linting, not build output.
- `build.rs` copies `default_inputs.yaml` to `OUT_DIR` and optionally runs Python+scipy to convert `src/core/wake/wake_velocity/turbopark_lookup_table.mat` to a Rust lookup table. Falls back gracefully if Python is unavailable.

## Code Style

### Imports & Organization
- Use absolute crate paths: `crate::core::wake::WakeModelManager`
- Group imports: std → crates → local (`use std::...; use crate::...;`)
- Doc comments: `///` for public API, `//!` for module-level

### Naming Conventions
- Types: `PascalCase` (structs, enums, traits)
- Functions/variables: `snake_case`
- Constants: `SCREAMING_SNAKE_CASE`
- Array types: `Array1`, `Array2`, `Array3`, `Array4` (see `types.rs`)
- Float type alias: `Float = f64`

### Error Handling
- Return type: `crate::Result<T>` = `anyhow::Result<T>`
- Use `anyhow::bail!("...")` for early returns
- Avoid `unwrap()` in production code (use `?` operator)
- Never suppress errors with `as any`

### Array Indexing
- 4D arrays: `[findex, turbine, y, z]`
- Use `ndarray::s![]` macro for slicing
- Clone arrays before mutation to avoid borrow checker issues

### Testing
- Tests go in `#[cfg(test)]` module at end of file
- Use `approx::assert_relative_eq!` for floating-point comparisons
- Regression tests in `src/core/wake/floris_regression_tests.rs` validate against Python FLORIS values

### Wake Models
- Traits: `VelocityModel`, `DeflectionModel`, `TurbulenceModel`, `CombinationModel`
- Implement in `src/core/wake/{submodule}/`
- Use `Box<dyn Trait>` for dynamic dispatch

### Operation Models
- `set_operation_model(name)` MUST be called AFTER `set_layout()` — `set_layout()` recreates the Farm and overwrites the operation model
- Supported: "simple", "cosine-loss", "simple-derating", "peak-shaving", "mixed", "awc", "unified-momentum", "controller-dependent"
- Power setpoints: Use `POWER_SETPOINT_DEFAULT` (1e12 from `operation_models`) as sentinel for "no limit". Zero is NOT unlimited — it means disabled.
- `farm.rs` has its own `RATED_POWER_DEFAULT = 5000.0` — do NOT confuse with `POWER_SETPOINT_DEFAULT`

### Known Limitations
- `get_turbine_powers()` bypasses the operation model for power calculation (uses inline Cp lookup). Power setpoint clamping is applied post-hoc.
- `get_turbine_thrust_coefficients()` uses `wind_speeds[fi]` (freestream) not rotor-averaged velocity, and doesn't go through the operation model.
- Disabled turbines via `set_operation(disable_turbines=...)` set power to 0 but don't trigger wake re-computation.
- `sample_flow_at_points()` only samples at findex=0 — loop through wind directions individually for multi-condition sampling.
- `sample_ti_at_points()` returns constant inflow TI — the turbulence field isn't populated after wake solving.

#### 实现

florus项目是根据python库floris(在floris文件夹下)实现的, 除了基本功能组织存在小差异外, 其大部分功能与floris一致, 尤其是提供给用户的接口,计算结果也要一致, 目前已经实现了基本功能, 请根据floris中的examples代码, 实现florus对应的examples代码, 要求双侧运行结果一致, 如果不一致就修改florus实现, 知道结果一致, 实现过程中不能使用简化实现,也不能使用占位符实现,  如果不知道如何实现, 可以参考floris中实现;o