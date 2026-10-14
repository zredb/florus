# FLORUS Agent Guidelines

## Build / Test / Lint

```bash
cargo build --release          # build
cargo test                     # all tests
cargo test <test_name>         # single test (e.g. solver::tests::test_axial_induction_low_ct)
cargo clippy -- -D warnings    # lint — ALWAYS use this, not build output
cargo fmt -- --check           # formatting check
```

### Build quirks

- `.cargo/config.toml` sets `rustflags = ["-A", "warnings"]` — **all warnings suppressed in `cargo build`**. Use `cargo clippy -- -D warnings` for actual lint results.
- `build.rs` copies `default_inputs.yaml` to `OUT_DIR` and optionally runs Python+scipy to convert `src/core/wake/wake_velocity/turbopark_lookup_table.mat` to a Rust lookup table. Falls back gracefully if Python unavailable.

## Architecture

- **Crate type**: Library + binary. `lib.rs` exports the public API; `main.rs` is a thin CLI.
- **Python reference**: `floris/` contains FLORIS v4.6 (Python). **This is the source of truth for behavior.** When adding features or fixing logic, check the Python source first and port accurately.
- **4D array convention**: All spatial fields use `[findex, turbine, y, z]` indexing (`ndarray::Array4`).
- **Wake model dispatch**: `Box<dyn Trait>` via `WakeModelManager`. Implementations in `src/core/wake/{velocity,deflection,turbulence,combination}/`.
- **Turbine types**: Loaded from YAML in `turbine_library/`. Cp/Ct surfaces from `.npz` (Python) or `.csv`.

### Key source layout

```
src/
├── floris_model.rs       # FlorisModel — main user-facing API (mirrors Python FlorisModel)
├── floris_config.rs      # Config parsing from YAML
├── types.rs              # Float = f64, Array1-4 type aliases, NumericDict
├── core/
│   ├── farm.rs           # Farm struct, turbine management
│   ├── flow_field.rs     # Flow field computation
│   ├── solver.rs         # Wake solver loop
│   ├── state.rs          # Mutable state for iterations
│   ├── grid/             # Turbine grid, flow grid
│   ├── turbines/         # Turbine types, Cp/Ct tables, operation models
│   └── wake/             # Wake models (velocity, deflection, turbulence, combination)
├── wind_data/            # TimeSeries, WindRose, WindRoseWRG
├── optimization/         # Yaw, layout, load optimization
├── visualization/        # Cut plane plotting
└── output_formatter.rs   # Output formatting
```

## Code Style

- **Imports**: Absolute crate paths (`crate::core::wake::WakeModelManager`), grouped: std → crates → local.
- **Doc comments**: `///` for public API, `//!` for module-level.
- **Types**: `PascalCase` structs/enums/traits, `snake_case` functions/variables, `SCREAMING_SNAKE_CASE` constants.
- **Float**: Use `Float` alias (= `f64`), not raw `f64`.
- **Arrays**: `Array1`, `Array2`, `Array3`, `Array4` from `types.rs`. Use `ndarray::s![]` for slicing. Clone before mutation to avoid borrow checker issues.
- **Error handling**: Return `crate::Result<T>` (= `anyhow::Result<T>`). Use `anyhow::bail!("...")` for early returns. Avoid `unwrap()` in production code.

## Testing

- Tests go in `#[cfg(test)]` module at end of file.
- Use `approx::assert_relative_eq!` for floating-point comparisons.
- **Regression tests** in `src/core/wake/floris_regression_tests.rs` validate against Python FLORIS values.
- **Consistency test design**: `tests/FLORUS_FLORIS_CONSISTENCY_TEST_DESIGN.md` documents the golden-value cross-validation strategy and test categories.
- **Golden values** live in `tests/golden_values/` (JSON files generated from Python FLORIS).
- Tolerances: 1% relative for power/AEP, 0.5% for Ct, 0.1% for velocity.

## Operation Models

- `set_operation_model(name)` MUST be called **AFTER** `set_layout()` — `set_layout()` recreates the Farm and overwrites the operation model.
- Supported: `"simple"`, `"cosine-loss"`, `"simple-derating"`, `"peak-shaving"`, `"mixed"`, `"awc"`, `"unified-momentum"`, `"controller-dependent"`.
- Power setpoints: Use `POWER_SETPOINT_DEFAULT` (1e12, from `operation_models`) as sentinel for "no limit". **Zero means disabled, not unlimited.**
- `farm.rs` previously carried a bogus `RATED_POWER_DEFAULT = 5000.0` used to seed power setpoints; it was removed in Bug #30 — `set_power_setpoints_to_ref_power` now uses `POWER_SETPOINT_DEFAULT` (1e12), matching Python.

## Wake Model Traits

- `VelocityModel`, `DeflectionModel`, `TurbulenceModel`, `CombinationModel`
- Implement in `src/core/wake/{submodule}/`.
- Use `Box<dyn Trait>` for dynamic dispatch.

## Examples

- Python examples: `floris/examples/`
- Rust examples: `examples/` (mirrors Python naming: `001_opening_floris_computing_power.rs`, etc.)
- Input YAMLs: `examples/inputs/` (gch.yaml, jensen.yaml, turbopark.yaml, etc.)
- **When porting examples from Python**: match the Rust output to Python output numerically. If results differ, fix the Rust implementation — do not use simplified or placeholder implementations.

## Known Limitations

- `get_turbine_powers()` bypasses the operation model for power calculation (uses inline Cp lookup). Power setpoint clamping is applied post-hoc.
- `get_turbine_thrust_coefficients()` uses `wind_speeds[fi]` (freestream) not rotor-averaged velocity, and doesn't go through the operation model.
- Floating-tilt correction (`compute_tilt_angles_for_floating_turbines` + `correct_cp_ct_for_tilt`) is not wired into operation models: Python applies it inside each model's `thrust_coefficient()` (e.g. controller_dependent_operation_model.py:262-269), Rust passes `correct_cp_ct_for_tilt: None`. Only affects floating turbines with tilt interpolation; fixed-bottom turbines are unaffected (`np.where(correct_cp_ct_for_tilt, ...)` keeps the original tilt).
- Disabled turbines via `set_operation(disable_turbines=...)` set power to 0; the solver DOES see the (near-zero) setpoint in Ct/AI since the Bug #35 wiring, but disable does not itself trigger wake re-computation (call `run()` after `set_operation`).
- `sample_flow_at_points()` only samples at findex=0 — loop through wind directions individually for multi-condition sampling.
- `sample_ti_at_points()` returns constant inflow TI — the turbulence field isn't populated after wake solving.
- `TurbineCubatureGrid`'s `average_method` (`"simple-cubature"`) + `cubature_weights` feed **only** the solver's thrust/axial-induction path (Python `core/solver.py`). `get_turbine_powers()` still averages with the cubic-mean default — Python `_get_turbine_powers` does not forward the grid's averaging to `power()`; do not wire the cubature weights into the power path (verified against Python: gauss_cubature_grid upstream power 1745305.95 W from cubic-mean vs 1735364.46 W from simple-cubature).

参考的python 代码在:D:\code\florus\floris下
