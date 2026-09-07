# FLORUS-FLORIS Consistency Test Design

## Goal

Ensure FLORUS (Rust) and FLORIS (Python) produce numerically consistent results across all public API functions. Tests run identical configurations and compare outputs within defined tolerances.

## Test Architecture

### Strategy: Golden Value Cross-Validation

1. **Generate golden values** from Python FLORIS for each test case
2. **Run same configuration** in FLORUS Rust
3. **Compare outputs** within tolerances (relative 1%, absolute fallback)
4. **Automate**: Script generates Python expected values, Rust tests consume them

### Tolerance Levels

| Metric | Tolerance | Justification |
|--------|-----------|---------------|
| Turbine power (W) | 1% relative OR 10W absolute | Floating point rounding, interpolation differences |
| Thrust coefficient | 0.5% relative | Ct lookup table interpolation |
| Axial induction | 1% relative | Derived from Ct |
| Rotor-averaged velocity | 0.1% relative | Grid resolution dependent |
| AEP (Wh) | 1% relative | Aggregate of per-condition errors |
| Farm power (W) | 1% relative | Sum of per-turbine errors |

---

## Test Categories

### Category 1: Wake Model Regression (Existing + Extended)

**Current state**: `floris_regression_tests.rs` covers TurbOPark, Gauss, Jensen, None with 4 wind speeds.

**Extend to match FLORIS regression test matrix**:

| Velocity Model | Deflection Model | Combination | Turbulence | Yaw | Test Name |
|----------------|-----------------|-------------|------------|-----|-----------|
| turbopark | gauss | fls | none | No | turbopark_regression ✅ exists |
| turbopark | gauss | fls | none | Yes | turbopark_yawed ✅ exists |
| gauss | gauss | fls | crespo_hernandez | No | gauss_regression ✅ exists |
| gauss | gauss | fls | crespo_hernandez | Yes | gauss_yawed ✅ exists |
| gauss | gauss | sosfs | crespo_hernandez | No | **gch_regression** (NEW) |
| gauss | gauss | sosfs | crespo_hernandez | Yes | **gch_yawed** (NEW) |
| jensen | jimenez | sosfs | crespo_hernandez | No | jensen_regression ✅ exists |
| jensen | jimenez | sosfs | crespo_hernandez | Yes | **jensen_yawed** (NEW) |
| turboparkgauss | gauss | sosfs | crespo_hernandez | No | **turboparkgauss_regression** (NEW) |
| none | none | fls | none | No | none_regression ✅ exists |
| cumulative_gauss_curl | gauss | sosfs | crespo_hernandez | No | **cc_regression** (NEW) |
| empirical_gauss | empirical_gauss | sosfs | wake_induced_mixing | No | **empirical_gauss_regression** (NEW) |

**Golden value source**: FLORIS `tests/reg_tests/` baseline arrays.

**Configuration** (matches FLORIS conftest.py):
- 3 turbines at [0, 630, 1260] m (5D, 10D spacing)
- 16 wind conditions: 4 directions × 4 speeds
- Directions: [270, 270, 270, 270, 360, 360, 360, 360, 285, 285, 285, 285, 315, 315, 315, 315]
- Speeds: [8, 9, 10, 11] repeated
- TI: 0.1 (constant)
- Solver: turbine_grid with 3 grid points
- Hub height: 90m, rotor diameter: 126m

---

### Category 2: FlorisModel High-Level API

**Source**: FLORIS `floris_model_integration_test.py`

#### 2.1 Basic Power Computation

| Test | Config | Expected (from FLORIS) |
|------|--------|----------------------|
| `test_single_condition_power` | 2 turbines, 270°, 8 m/s, TI=0.06 | farm_power = 3,507,908.92 W |
| `test_multi_condition_power` | 3 turbines, 4 conditions | Per-turbine powers match baseline |
| `test_no_wake_power` | run_no_wake vs run | No-wake > wake power |
| `test_get_farm_power` | Sum of turbine powers | matches sum |
| `test_get_farm_AEP` | TimeSeries, uniform freq | AEP = sum(powers * freq) * 8760 |

#### 2.2 Set/Reset Operations

| Test | Action | Expected |
|------|--------|----------|
| `test_set_layout` | Change from 3 to 2 turbines | n_turbines=2, power recomputed |
| `test_set_wind_conditions` | Change speed from 8 to 12 m/s | Higher power |
| `test_set_yaw_angles` | T0 yaw=10° | T0 power reduced, T1 power increased |
| `test_set_power_setpoints` | Limit to 3MW | Power capped at 3MW |
| `test_set_disable_turbines` | Disable T0 | T0 power=0 |
| `test_reset_operation` | After set_operation | Returns to defaults |
| `test_set_operation_model` | Switch to "simple" | Model changes |
| `test_set_wind_shear` | Change shear exponent | Velocity profile changes |
| `test_set_air_density` | Change from 1.225 to 1.0 | Power scales with density^(1/3) |
| `test_set_reference_wind_height` | Change ref height | Wind speed at ref height changes |

#### 2.3 Flow Sampling

| Test | Config | Expected |
|------|--------|----------|
| `test_sample_flow_at_points` | Met mast at known position | Velocity matches grid interpolation |
| `test_sample_ti_at_points` | Same position | TI matches inflow TI (known limitation) |
| `test_calculate_horizontal_plane` | Hub height slice | 2D plane shape correct, values in range |
| `test_calculate_y_plane` | Streamwise vertical slice | Shape and values correct |
| `test_calculate_cross_plane` | Perpendicular slice at 500m | Shape and values correct |

#### 2.4 Turbine Data

| Test | Expected |
|------|----------|
| `test_get_turbine_powers` | Powers match baseline |
| `test_get_turbine_thrust_coefficients` | Ct values from lookup table |
| `test_get_turbine_ais` | Axial induction derived from Ct |
| `test_turbine_average_velocities` | Rotor-averaged velocities |
| `test_n_findex` | Returns correct count |
| `test_n_turbines` | Returns correct count |
| `test_layout_x_layout_y` | Returns correct coordinates |

---

### Category 3: Wind Data Types

#### 3.1 TimeSeries

| Test | Expected |
|------|----------|
| `test_timeseries_new` | Correct dimensions, default values |
| `test_timeseries_to_wind_rose` | Conversion produces correct bins |
| `test_timeseries_assign_ti_iec` | TI follows IEC 61400-1 formula |
| `test_timeseries_assign_value_piecewise` | Value follows piecewise linear |
| `test_timeseries_set_model` | fmodel.set_wind_data(ts) works |

#### 3.2 WindRose

| Test | Expected |
|------|----------|
| `test_windrose_new` | Correct shape validation |
| `test_windrose_downsample` | Reduces direction resolution |
| `test_windrose_upsample` | Increases direction resolution |
| `test_windrose_assign_ti_iec` | TI follows IEC formula |
| `test_windrose_assign_value_piecewise` | Value follows piecewise linear |
| `test_windrose_set_model` | fmodel.set_wind_data(wr) works |
| `test_windrose_frequencies` | Frequencies sum to ~1.0 |

#### 3.3 WindRoseWRG

| Test | Expected |
|------|----------|
| `test_wrg_read` | Parses WRG file correctly |
| `test_wrg_get_wind_rose_at_point` | Interpolated values reasonable |
| `test_wrg_set_layout` | Per-turbine wind roses created |
| `test_wrg_set_model` | fmodel.set_wind_data(wrg) works |

---

### Category 4: Operation Models

| Test | Model | Config | Expected |
|------|-------|--------|----------|
| `test_simple_power` | simple | 10 m/s, no yaw | Power from Cp table |
| `test_cosine_loss_power` | cosine-loss | 10 m/s, yaw=10° | Power reduced by cos^p |
| `test_simple_derating` | simple-derating | SP=3MW, available=5MW | Power capped at 3MW |
| `test_peak_shaving` | peak-shaving | 11 m/s, TI=0.1 | Ct reduced near rated |
| `test_mixed_disabled` | mixed | disable_turbine=true | Power=0 |
| `test_mixed_derating` | mixed | SP=2MW on T0 | T0 capped, T1 free |
| `test_awc_baseline` | awc | baseline mode | Same as simple |
| `test_controller_dependent` | controller-dependent | 10 m/s | Matches controller model |
| `test_unified_momentum` | unified-momentum | 10 m/s | Matches UMM model |

---

### Category 5: Optimization

#### 5.1 Yaw Optimization

| Test | Method | Expected |
|------|--------|----------|
| `test_sr_optimize_single` | SerialRefine, 3 turbines, 270° | Optimized power > baseline |
| `test_sr_optimize_multiple_wd` | SerialRefine, WindRose | Per-WD yaw angles |
| `test_geometric_optimize` | Geometric, 3 turbines | Fast optimization |
| `test_optimize_yaw_angles` | Top-level function | Returns yaw angles |

#### 5.2 Layout Optimization

| Test | Method | Expected |
|------|--------|----------|
| `test_random_search` | RandomSearch, 4 turbines | Improved layout |
| `test_boundary_grid` | BoundaryGrid, square | Points inside boundary |
| `test_min_distance` | All points > min_dist | No violations |

---

### Category 6: Visualization (Shape/Existence Only)

| Test | Expected |
|------|----------|
| `test_visualize_cut_plane` | PNG file created, size > 0 |
| `test_visualize_horizontal_plane` | PNG file created |
| `test_layout_visualization` | PNG files created |

---

## Implementation Plan

### Phase 1: Golden Value Extraction (Python)

Create `tests/generate_golden_values.py` that:
1. Runs each test configuration in FLORIS
2. Saves expected outputs as JSON/CSV files in `tests/golden_values/`
3. Files: `wake_regression.json`, `floris_model_api.json`, `wind_data.json`, etc.

### Phase 2: Rust Test Implementation

Add tests to existing files or create new test modules:

1. **Extend `floris_regression_tests.rs`** with 6 new wake model configurations
2. **Create `floris_model_api_tests.rs`** for high-level API tests
3. **Create `wind_data_consistency_tests.rs`** for wind data type tests
4. **Create `operation_model_tests.rs`** for operation model tests
5. **Create `optimization_consistency_tests.rs`** for optimization tests

### Phase 3: Automated Comparison

Add `tests/consistency_test.rs` integration test that:
1. Reads golden values from JSON
2. Runs FLORUS with same configurations
3. Compares within tolerances
4. Reports pass/fail with details

---

## Test Data Files

```
tests/
├── golden_values/
│   ├── wake_regression/
│   │   ├── gauss_baseline.json
│   │   ├── jensen_baseline.json
│   │   ├── turbopark_baseline.json
│   │   ├── turboparkgauss_baseline.json
│   │   ├── cc_baseline.json
│   │   ├── empirical_gauss_baseline.json
│   │   └── none_baseline.json
│   ├── floris_model_api/
│   │   ├── power_computation.json
│   │   ├── set_operations.json
│   │   └── flow_sampling.json
│   ├── wind_data/
│   │   ├── timeseries.json
│   │   ├── windrose.json
│   │   └── wrg.json
│   └── optimization/
│       ├── yaw_optimization.json
│       └── layout_optimization.json
├── generate_golden_values.py
└── consistency_test.rs
```

---

## Critical File Paths

### FLORIS Python (golden value source)
- `floris/tests/conftest.py` - fixtures, constants
- `floris/tests/reg_tests/gauss_regression_test.py` - Gauss baselines
- `floris/tests/reg_tests/jensen_jimenez_regression_test.py` - Jensen baselines
- `floris/tests/reg_tests/turbopark_regression_test.py` - TurbOPark baselines
- `floris/tests/reg_tests/turboparkgauss_regression_test.py` - TurboParkGauss baselines
- `floris/tests/reg_tests/cumulative_curl_regression_test.py` - CC baselines
- `floris/tests/reg_tests/empirical_gauss_regression_test.py` - Empirical Gauss baselines
- `floris/tests/reg_tests/none_regression_test.py` - No-wake baselines
- `floris/tests/floris_model_integration_test.py` - API integration tests
- `floris/tests/turbine_operation_models_unit_test.py` - Operation model tests

### FLORUS Rust (test implementation)
- `src/core/wake/floris_regression_tests.rs` - existing regression tests
- `src/core/wake/wake_tests.rs` - wake physics tests
- `src/floris_model.rs` (tests module) - model API tests
- `src/wind_data/` (tests) - wind data tests
- `src/optimization/` (tests) - optimization tests
