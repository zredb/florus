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

- **Multi-dimensional Cp/Ct tables** (e.g. `iea_15MW_multi_dim_cp_ct.yaml`): the operation-model dispatch (`Turbine::power/thrust_coefficient/axial_induction`) supports only 1-D `power`/`thrust_coefficient` curves. `get_turbine_powers()` falls back to the legacy inline estimate (with a warning) for such fleets; `get_turbine_ais()` returns an error. Python routes multi-dim tables through a separate `turbine.py` branch that is not ported.
- Disabled turbines via `set_operation(disable_turbines=...)` set power to 0; the solver DOES see the (near-zero) setpoint in Ct/AI since the Bug #35 wiring, but disable does not itself trigger wake re-computation (call `run()` after `set_operation` — Python requires the same). Query APIs now refuse stale state: `get_turbine_powers()` errors with Python's "without first running" message after any state-invalidating change.
- `solve_for_points`/`sample_flow_at_points`/`sample_ti_at_points` support the sequential-model family only; `cc` and `empirical_gauss` full-flow solver variants are not ported (sampling errors for those models; `turbopark` is refused exactly like Python).
- `get_operation_models()` is a placeholder using hardcoded Ct thresholds, not the operation models' state.
- **Many Rust examples are narrative placeholders or use their own scenarios instead of faithful ports** of the Python examples (control_optimization_003–008, layout_optimization_001–005, heterogeneous_001–004, multidim_001–003, wind_resource_grid_001–003, load_optimization, floating_002, uncertain_002). Aligned 1:1 so far: `001_opening_floris_computing_power` (values bit-identical), `get_flow_001`/`011`, `006_get_farm_aep` (AEP within 0.03%, wake losses identical), `control_optimization_003`–`008` (SR/Geometric optima bit-identical; AEP within 0.1%), `layout_optimization_004` (gridded layout bit-identical). See `output/examples_comparison.txt` for the full side-by-side inventory.
- **Engine gaps that block the remaining example ports** (not example-only work):
  - *Heterogeneous inflow is not wired into the solver* — `HeterogeneousMap` exists as a data structure and is plumbed through the wind-data layer, but `FlowField::initialize_flow_field` never applies the speed multipliers, so `heterogeneous_001`–`004` and `layout_optimization_002` cannot reproduce Python.
  - *Multi-dimensional Cp/Ct tables are not usable through the operation models* — `Turbine::power`/`thrust_coefficient`/`axial_induction` bail for `CpCtTable::MultiDimensional` (only the legacy inline power path reads them), so `multidim_001`–`003` and `floating_002` (iea_15MW_floating is a multi-dim turbine) cannot reproduce Python.
  - *Layout optimizers differ in algorithm*: FLORUS's `LayoutOptimizationScipy`/`RandomSearch` are coordinate-descent / heuristic-spread, while Python uses SLSQP and a genetic algorithm; only the deterministic `LayoutOptimizationGridded` was ported exactly.
  - *The WRG parser does not reproduce Python's per-turbine roses* — the grid metadata, the Weibull frequency construction and the layout application are implemented (and the bilinear interpolation was fixed to index the `[nx, ny]` data correctly, which removed an out-of-bounds panic at grid edges), but the file is parsed with a different grid-point layout, so `wind_resource_grid_002` yields different expected powers than Python.
- *Yaw optimizers were re-implemented to Python parity* (`YawOptimizationSR` = Ny_passes [5,4] two-pass search, `YawOptimizationGeometric` = analytic trapezoid, `YawOptimizationScipy` = per-condition block-coordinate search); `set_yaw_angles` now invalidates the solved state so the solver re-sorts the yaw angles (Python `set()` resets `core.state`).
- `TurbineCubatureGrid`'s `average_method` (`"simple-cubature"`) + `cubature_weights` feed **only** the solver's thrust/axial-induction path (Python `core/solver.py`). `get_turbine_powers()` still averages with the cubic-mean default — Python `_get_turbine_powers` does not forward the grid's averaging to `power()`; do not wire the cubature weights into the power path (verified against Python: gauss_cubature_grid upstream power 1745305.95 W from cubic-mean vs 1735364.46 W from simple-cubature).

## Query API semantics (fixed 2026-10, cross-validated in `tests/cross_validation_query_apis.rs`)

- `get_turbine_powers()` routes through the operation model over the SOLVED field with the cubic-mean default (Python `floris_model.py:554-570`), returns layout-order `(n_findex, n_turbines)` in Watts and a `Result` (state check mirrors Python's RuntimeError).
- `get_turbine_thrust_coefficients()` / `get_turbine_ais()` route through `farm.thrust_coefficient` / `farm.axial_induction` with the grid's averaging method + cubature weights, then unsort to layout order (Python `floris_model.py:1017-1057`).
- `turbine_average_velocities()` (Python property parity, `floris_model.py:1798-1804`): rotor average over the solved field with the grid's method/weights, layout order — previously returned a single rotor point in sorted order. `get_turbine_velocities()` remains a layout-order FIRST-POINT sample (no Python equivalent).
- `get_farm_aep(wind_data, hours)` uses the `unpack_freq` contract (Python `get_farm_AEP`): per-condition frequency column `(n_findex, 1)` weights each condition's farm power; a `(n_findex, n_turbines)` table (WindRoseByTurbine) weights each turbine's power. `WindRose::unpack` returns the flattened, zero-freq-filtered frequency aligned with `n_findex`; `TimeSeries::unpack` returns uniform `1/N`.
- `Farm::set_layout` resets the operation arrays (yaw/tilt/setpoints/AWC) when the turbine COUNT changes — previously they kept the old count and `expand_farm_properties` panicked (index out of bounds).
- Floating tilt: `correct_cp_ct_for_tilt` + `floating_tilt_table` are wired into cosine-loss power (interp at the post-density/yaw effective velocity), cosine-loss & controller-dependent thrust (interp at the uncorrected rotor average), and axial induction. Only `cosine-loss` and `controller-dependent` use the interp — matching Python. Test turbine: `nrel_5MW_floating`.
- `sample_flow_at_points` / `sample_ti_at_points` solve the wake ON the points for ALL findices via `Core::solve_for_points` (PointsGrid + `full_flow_sequential_solver`, Python `core.py:230-266`) and return `(n_findex, n_points)`. The solve runs on a scratch flow field; the model state stays intact (Python clobbers `self.flow_field` instead).

参考的python 代码在:D:\code\florus\floris下
