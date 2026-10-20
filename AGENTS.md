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

- **Multi-dimensional Cp/Ct tables are supported** (2026-10): the multi-dimensional data file is loaded as a *set of one-dimensional tables keyed by the condition tuple* (`MultiDimConditionTables`, Python `turbine.py:_initialize_multidim_power_thrust_table`), and for every wind condition the nearest specified value of each condition dimension selects the table (Python `turbine.py:_select_multidim_condition`). The operation models are called once per selected table and the results are masked per condition, exactly as Python does. New APIs: `FlorisModel::set_multidim_conditions(...)` (scalar or per-condition values), `flow_field.multidim_conditions` in the input file, and `TimeSeries::multidim_conditions` (the wind-data route overrides the input-file values, as in Python).
- Disabled turbines via `set_operation(disable_turbines=...)` set power to 0; the solver DOES see the (near-zero) setpoint in Ct/AI since the Bug #35 wiring, but disable does not itself trigger wake re-computation (call `run()` after `set_operation` — Python requires the same). Query APIs now refuse stale state: `get_turbine_powers()` errors with Python's "without first running" message after any state-invalidating change.
- `solve_for_points`/`sample_flow_at_points`/`sample_ti_at_points` support the sequential-model family only; `cc` and `empirical_gauss` full-flow solver variants are not ported (sampling errors for those models; `turbopark` is refused exactly like Python).
- `get_operation_models()` is a placeholder using hardcoded Ct thresholds, not the operation models' state.
- **Many Rust examples are narrative placeholders or use their own scenarios instead of faithful ports** of the Python examples (layout_optimization_001–005, wind_resource_grid_001–003, load_optimization, floating_002, uncertain_002). Aligned 1:1 so far: `001_opening_floris_computing_power` (values bit-identical), `get_flow_001`/`011`, `006_get_farm_aep` (AEP within 0.03%, wake losses identical), `control_optimization_003`–`008` (SR/Geometric optima bit-identical; AEP within 0.1%), `layout_optimization_004` (gridded layout bit-identical), `heterogeneous_001`–`004` (powers bit-identical to Python; the 2D/3D flow-field statistics agree within 0.02 m/s), `multidim_001`–`003` (powers bit-identical to Python). See `output/examples_comparison.txt` for the full side-by-side inventory.
- **Engine gaps that block the remaining example ports** (not example-only work):
  - *Layout optimizers differ in algorithm*: FLORUS's `LayoutOptimizationScipy`/`RandomSearch` are coordinate-descent / heuristic-spread, while Python uses SLSQP and a genetic algorithm; only the deterministic `LayoutOptimizationGridded` was ported exactly.
  - *The WRG parser does not reproduce Python's per-turbine roses* — the grid metadata, the Weibull frequency construction and the layout application are implemented (and the bilinear interpolation was fixed to index the `[nx, ny]` data correctly, which removed an out-of-bounds panic at grid edges), but the file is parsed with a different grid-point layout, so `wind_resource_grid_002` yields different expected powers than Python.
  - *No `ApproxFlorisModel`* — Python's `uncertain_002` example uses the surrogate model (`floris/simulation/approx.py`), which has no Rust equivalent.
- *Heterogeneous inflow is wired into the solver* (2026-10): the speed multipliers are interpolated onto the grid points in the **inertial frame** and multiply the sheared initial field, exactly as Python's `initialize_velocity_field` does (`floris/core/flow_field.py`). The application lives in `FlowField::initialize_flow_field_with_coordinates`, so every flow-field initialization (turbine grid, visualization grid, points grid) applies the multipliers — Python's `copy.deepcopy(flow_field)` behaviour. New API: `FlorisModel::set_heterogeneous_inflow_config(HeterogeneousInflowConfig)`; the configuration can also be given in the input file under `flow_field.heterogeneous_inflow_config` (2D and 3D, with `speed_multipliers`, `x`, `y`, optional `z`). Interpolation: `nearest` matches scipy exactly; `linear` uses a Bowyer-Watson Delaunay triangulation with barycentric coordinates (scipy `LinearNDInterpolator`), which can pick a different diagonal than qhull for co-circular points; the 3D case uses trilinear interpolation on structured grids (identical to scipy for affinely varying data, e.g. the FLORIS examples) and falls back to nearest neighbours for scattered locations.
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
