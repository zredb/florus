#!/usr/bin/env python3
"""
FLORUS-FLORIS Golden Value Generator

This script generates golden values from Python FLORIS for cross-validation
with the Rust FLORUS implementation. It produces JSON files that can be
consumed by Rust tests.

Usage:
    python tests/generate_golden_values.py                    # Generate all
    python tests/generate_golden_values.py --category wake    # Generate wake regression only
    python tests/generate_golden_values.py --output /path     # Custom output directory
"""

import argparse
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

# Add floris to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "floris"))

from floris.floris_model import FlorisModel
from floris.core import Core
from floris.wind_data import TimeSeries, WindRose


# =============================================================================
# Test Configuration (must match Rust tests)
# =============================================================================

WIND_DIRECTIONS = [
    270.0, 270.0, 270.0, 270.0,
    360.0, 360.0, 360.0, 360.0,
    285.0, 285.0, 285.0, 285.0,
    315.0, 315.0, 315.0, 315.0,
]

WIND_SPEEDS = [
    8.0, 9.0, 10.0, 11.0,
    8.0, 9.0, 10.0, 11.0,
    8.0, 9.0, 10.0, 11.0,
    8.0, 9.0, 10.0, 11.0,
]

TURBULENCE_INTENSITIES = [0.1] * 16

N_FINDEX = 16
N_TURBINES = 3
ROTOR_DIAMETER = 126.0
TURBINE_GRID_RESOLUTION = 2

X_COORDS = [0.0, 5.0 * ROTOR_DIAMETER, 10.0 * ROTOR_DIAMETER]
Y_COORDS = [0.0, 0.0, 0.0]
Z_COORDS = [90.0, 90.0, 90.0]


# =============================================================================
# Helper Functions
# =============================================================================

def create_floris_model(
    velocity_model: str = "gauss",
    deflection_model: str = "gauss",
    combination_model: str = "fls",
    turbulence_model: str = "crespo_hernandez",
    wind_speeds: Optional[List[float]] = None,
    wind_directions: Optional[List[float]] = None,
    turbulence_intensities: Optional[List[float]] = None,
    layout_x: Optional[List[float]] = None,
    layout_y: Optional[List[float]] = None,
    turbine_type: str = "nrel_5MW",
    operation_model: str = "cosine-loss",
    yaw_angles: Optional[List[List[float]]] = None,
    turbine_grid_points: Optional[int] = None,
) -> FlorisModel:
    """Create a FlorisModel with specified configuration."""
    
    if wind_speeds is None:
        wind_speeds = WIND_SPEEDS
    if wind_directions is None:
        wind_directions = WIND_DIRECTIONS
    if turbulence_intensities is None:
        turbulence_intensities = TURBULENCE_INTENSITIES
    if layout_x is None:
        layout_x = X_COORDS
    if layout_y is None:
        layout_y = Y_COORDS
    
    n_findex = len(wind_speeds)
    
    config = {
        "name": "golden_value_test",
        "description": "Test configuration for golden value generation",
        "floris_version": "v4",
        "logging": {
            "console": {
                "enable": False,
                "level": "WARNING",
            },
            "file": {
                "enable": False,
                "level": "WARNING",
            },
        },
        "solver": {
            "type": "turbine_grid",
            "turbine_grid_points": (
                turbine_grid_points
                if turbine_grid_points is not None
                else TURBINE_GRID_RESOLUTION
            ),
        },
        "farm": {
            "layout_x": layout_x,
            "layout_y": layout_y,
            "turbine_type": [turbine_type] * len(layout_x),
        },
        "flow_field": {
            "wind_speeds": wind_speeds,
            "wind_directions": wind_directions,
            "turbulence_intensities": turbulence_intensities,
            "wind_shear": 0.12,
            "wind_veer": 0.0,
            "air_density": 1.225,
            "reference_wind_height": 90.0,
        },
        "wake": {
            "model_strings": {
                "velocity_model": velocity_model,
                "deflection_model": deflection_model,
                "combination_model": combination_model,
                "turbulence_model": turbulence_model,
            },
            "enable_secondary_steering": True,
            "enable_yaw_added_recovery": True,
            "enable_transverse_velocities": True,
            "enable_active_wake_mixing": False,
            "wake_deflection_parameters": {
                "gauss": {
                    "ad": 0.0,
                    "alpha": 0.58,
                    "bd": 0.0,
                    "beta": 0.077,
                    "dm": 1.0,
                    "ka": 0.38,
                    "kb": 0.004,
                },
                "jimenez": {
                    "ad": 0.0,
                    "bd": 0.0,
                    "kd": 0.05,
                },
                "empirical_gauss": {
                    "horizontal_deflection_gain_D": 3.0,
                    "vertical_deflection_gain_D": -1,
                    "deflection_rate": 22,
                    "mixing_gain_deflection": 0.0,
                    "yaw_added_mixing_gain": 0.0,
                },
            },
            "wake_velocity_parameters": {
                "gauss": {
                    "alpha": 0.58,
                    "beta": 0.077,
                    "ka": 0.38,
                    "kb": 0.004,
                },
                "jensen": {
                    "we": 0.05,
                },
                "turbopark": {
                    "A": 0.04,
                    "sigma_max_rel": 4.0,
                },
                "turboparkgauss": {
                    "A": 0.04,
                    "include_mirror_wake": True,
                },
                "cc": {
                    "a_s": 0.179367259,
                    "b_s": 0.0118889215,
                    "c_s1": 0.0563691592,
                    "c_s2": 0.13290157,
                    "a_f": 3.11,
                    "b_f": -0.68,
                    "c_f": 2.41,
                    "alpha_mod": 1.0,
                },
                "empirical_gauss": {
                    "wake_expansion_rates": [0.023, 0.008],
                    "breakpoints_D": [10],
                    "sigma_0_D": 0.28,
                    "smoothing_length_D": 2.0,
                    "mixing_gain_velocity": 2.0,
                    "awc_wake_exp": 1.2,
                    "awc_wake_denominator": 400,
                },
            },
            "wake_turbulence_parameters": {
                "crespo_hernandez": {
                    "initial": 0.1,
                    "constant": 0.5,
                    "ai": 0.8,
                    "downstream": -0.32,
                },
                "wake_induced_mixing": {
                    "atmospheric_ti_gain": 0.0,
                },
            },
        },
    }
    
    fmodel = FlorisModel(config)
    
    if operation_model != "cosine-loss":
        fmodel.set_operation_model(operation_model)
    
    if yaw_angles is not None:
        fmodel.set(yaw_angles=yaw_angles)
    
    if velocity_model in ["turbopark", "empirical_gauss"]:
        fmodel.core.wake.enable_secondary_steering = False
        fmodel.core.wake.enable_transverse_velocities = False
        fmodel.core.wake.enable_yaw_added_recovery = False
    
    return fmodel


def extract_results(fmodel: FlorisModel) -> Dict[str, Any]:
    """Extract results from a FlorisModel after running."""
    fmodel.run()
    
    results = {
        "turbine_powers": fmodel.get_turbine_powers().tolist(),
        "farm_power": fmodel.get_farm_power().tolist(),
        "n_findex": fmodel.core.flow_field.n_findex,
        "n_turbines": fmodel.core.farm.n_turbines,
    }
    
    # Extract velocities if available
    if hasattr(fmodel.core.flow_field, 'u'):
        results["velocities_shape"] = list(fmodel.core.flow_field.u.shape)
    
    return results


def save_json(data: Any, filepath: Path) -> None:
    """Save data to JSON file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"  Saved: {filepath}")


# =============================================================================
# Wake Model Regression Tests
# =============================================================================

def generate_wake_regression_values(output_dir: Path) -> None:
    """Generate golden values for wake model regression tests."""
    print("\n=== Generating Wake Model Regression Values ===")
    
    wake_configs = [
        # (velocity_model, deflection_model, combination_model, turbulence_model, yawed)
        ("turbopark", "gauss", "fls", "none", False),
        ("turbopark", "gauss", "fls", "none", True),
        ("gauss", "gauss", "fls", "crespo_hernandez", False),
        ("gauss", "gauss", "fls", "crespo_hernandez", True),
        ("gauss", "gauss", "sosfs", "crespo_hernandez", False),
        ("gauss", "gauss", "sosfs", "crespo_hernandez", True),
        ("jensen", "jimenez", "sosfs", "none", False),
        ("jensen", "jimenez", "sosfs", "none", True),
        ("turboparkgauss", "gauss", "sosfs", "crespo_hernandez", False),
        ("cc", "gauss", "sosfs", "crespo_hernandez", False),
        ("empirical_gauss", "empirical_gauss", "sosfs", "wake_induced_mixing", False),
    ]
    
    for vel_model, def_model, comb_model, turb_model, yawed in wake_configs:
        test_name = f"{vel_model}_{'yawed' if yawed else 'regression'}"
        print(f"\n  Generating {test_name}...")
        
        yaw_angles = None
        if yawed:
            # Upstream turbine yawed 5 degrees
            yaw_angles = [[5.0, 0.0, 0.0]] * N_FINDEX
        
        try:
            fmodel = create_floris_model(
                velocity_model=vel_model,
                deflection_model=def_model,
                combination_model=comb_model,
                turbulence_model=turb_model,
                yaw_angles=yaw_angles,
            )
            
            results = extract_results(fmodel)
            results["config"] = {
                "velocity_model": vel_model,
                "deflection_model": def_model,
                "combination_model": comb_model,
                "turbulence_model": turb_model,
                "yawed": yawed,
            }
            
            save_json(results, output_dir / "wake_regression" / f"{test_name}.json")
            
        except Exception as e:
            print(f"    ERROR: {e}")


def generate_wake_ti4d_value(output_dir: Path) -> None:
    """Generate the staggered-layout golden value exercising per-point 4D TI.

    Turbine 1 is offset laterally (3 D) and turbine 0 is yawed 5 deg so the
    wake-added-turbulence field evolves per point and feeds back into the
    velocity models through the evolving 4D TI field (Front 2).
    """
    print("\n=== Generating Staggered-TI Wake Golden Value ===")

    d = 126.0
    fmodel = create_floris_model(
        velocity_model="gauss",
        deflection_model="gauss",
        combination_model="sosfs",
        turbulence_model="crespo_hernandez",
        wind_speeds=[8.0, 10.0],
        wind_directions=[270.0, 270.0],
        turbulence_intensities=[0.1, 0.1],
        layout_x=[0.0, 5 * d, 10 * d],
        layout_y=[0.0, 3 * d, 0.0],
        yaw_angles=[[5.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
        turbine_grid_points=3,
    )

    results = extract_results(fmodel)
    results["config"] = {
        "velocity_model": "gauss",
        "deflection_model": "gauss",
        "combination_model": "sosfs",
        "turbulence_model": "crespo_hernandez",
        "yawed": True,
    }
    save_json(results, output_dir / "wake_regression" / "gauss_staggered_ti.json")


def generate_wake_hetero_sort_value(output_dir: Path) -> None:
    """Generate the heterogeneous-fleet golden value with a per-findex sort
    order change (Bug #34).

    Turbines 0 and 2 are nrel_5MW, turbine 1 is iea_10MW with a large lateral
    offset. At 270 deg the sort order is 0, 1, 2; at 290 deg the rotated
    coordinates reorder to 1, 0, 2, so each findex maps different layout
    turbines to the sorted positions — exercising per-findex rotor diameter /
    hub height and Cp/Ct curve dispatch.
    """
    print("\n=== Generating Heterogeneous-Sort Wake Golden Value ===")

    fmodel = FlorisModel(
        {
            "name": "hetero_sort",
            "description": "heterogeneous fleet with per-findex sort order change",
            "floris_version": "v4",
            "logging": {
                "console": {"enable": False, "level": "WARNING"},
                "file": {"enable": False, "level": "WARNING"},
            },
            "solver": {"type": "turbine_grid", "turbine_grid_points": 2},
            "farm": {
                "layout_x": [0.0, 252.0, 504.0],
                "layout_y": [0.0, 700.0, 0.0],
                "turbine_type": ["nrel_5MW", "iea_10MW", "nrel_5MW"],
            },
            "flow_field": {
                "wind_speeds": [8.0, 8.0],
                "wind_directions": [270.0, 290.0],
                "turbulence_intensities": [0.1, 0.1],
                "wind_shear": 0.12,
                "wind_veer": 0.0,
                "air_density": 1.225,
                "reference_wind_height": 90.0,
            },
            "wake": {
                "model_strings": {
                    "velocity_model": "gauss",
                    "deflection_model": "gauss",
                    "combination_model": "sosfs",
                    "turbulence_model": "crespo_hernandez",
                },
                "enable_secondary_steering": True,
                "enable_yaw_added_recovery": True,
                "enable_transverse_velocities": True,
                "enable_active_wake_mixing": False,
                "wake_deflection_parameters": {
                    "gauss": {
                        "ad": 0.0,
                        "alpha": 0.58,
                        "bd": 0.0,
                        "beta": 0.077,
                        "dm": 1.0,
                        "ka": 0.38,
                        "kb": 0.004,
                    },
                },
                "wake_velocity_parameters": {
                    "gauss": {"ka": 0.38, "kb": 0.004, "alpha": 0.58, "beta": 0.077},
                },
                "wake_turbulence_parameters": {
                    "crespo_hernandez": {
                        "initial": 0.1,
                        "constant": 0.5,
                        "ai": 0.8,
                        "downstream": -0.32,
                    },
                },
            },
        }
    )
    fmodel.set(yaw_angles=[[5.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    fmodel.run()

    powers = fmodel.get_turbine_powers()
    results = {
        "turbine_powers": [[float(p) for p in row] for row in powers.tolist()],
        "farm_power": [float(p) for p in fmodel.get_farm_power().tolist()],
        "n_findex": 2,
        "n_turbines": 3,
        "velocities_shape": list(fmodel.core.flow_field.u.shape),
        "config": {
            "velocity_model": "gauss",
            "deflection_model": "gauss",
            "combination_model": "sosfs",
            "turbulence_model": "crespo_hernandez",
            "yawed": True,
        },
    }
    save_json(results, output_dir / "wake_regression" / "gauss_hetero_sort.json")


# =============================================================================
# FlorisModel API Tests
# =============================================================================

def generate_floris_model_api_values(output_dir: Path) -> None:
    """Generate golden values for FlorisModel API tests."""
    print("\n=== Generating FlorisModel API Values ===")
    
    # 1. Single condition power
    print("\n  Generating single_condition_power...")
    try:
        fmodel = create_floris_model(
            wind_speeds=[8.0],
            wind_directions=[270.0],
            turbulence_intensities=[0.06],
            layout_x=[0.0, 630.0],
            layout_y=[0.0, 0.0],
        )
        results = extract_results(fmodel)
        save_json(results, output_dir / "floris_model_api" / "single_condition_power.json")
    except Exception as e:
        print(f"    ERROR: {e}")
    
    # 2. No wake power
    print("\n  Generating no_wake_power...")
    try:
        fmodel = create_floris_model(
            wind_speeds=[8.0],
            wind_directions=[270.0],
            turbulence_intensities=[0.06],
            layout_x=[0.0, 630.0],
            layout_y=[0.0, 0.0],
        )
        fmodel.run_no_wake()
        results = {
            "turbine_powers": fmodel.get_turbine_powers().tolist(),
            "farm_power": fmodel.get_farm_power().tolist(),
        }
        save_json(results, output_dir / "floris_model_api" / "no_wake_power.json")
    except Exception as e:
        print(f"    ERROR: {e}")
    
    # 3. Yaw effect
    print("\n  Generating yaw_effect...")
    try:
        fmodel = create_floris_model(
            wind_speeds=[8.0],
            wind_directions=[270.0],
            turbulence_intensities=[0.06],
            layout_x=[0.0, 630.0],
            layout_y=[0.0, 0.0],
            yaw_angles=[[10.0, 0.0]],
        )
        results = extract_results(fmodel)
        save_json(results, output_dir / "floris_model_api" / "yaw_effect.json")
    except Exception as e:
        print(f"    ERROR: {e}")
    
    # 4. Disable turbines
    print("\n  Generating disable_turbines...")
    try:
        fmodel = create_floris_model(
            wind_speeds=[8.0, 8.0],
            wind_directions=[270.0, 270.0],
            turbulence_intensities=[0.06, 0.06],
            layout_x=[0.0, 630.0, 1260.0],
            layout_y=[0.0, 0.0, 0.0],
        )
        fmodel.run()
        powers_before = fmodel.get_turbine_powers().tolist()
        
        # Disable first turbine
        fmodel.set(disable_turbines=[[True, False, False], [True, False, False]])
        fmodel.run()
        powers_after = fmodel.get_turbine_powers().tolist()
        
        results = {
            "powers_before": powers_before,
            "powers_after": powers_after,
        }
        save_json(results, output_dir / "floris_model_api" / "disable_turbines.json")
    except Exception as e:
        print(f"    ERROR: {e}")
    
    # 5. Air density effect
    print("\n  Generating air_density_effect...")
    try:
        fmodel_1225 = create_floris_model(
            wind_speeds=[10.0],
            wind_directions=[270.0],
            turbulence_intensities=[0.06],
        )
        fmodel_1225.core.flow_field.air_density = 1.225
        fmodel_1225.run()
        power_1225 = fmodel_1225.get_turbine_powers().tolist()
        
        fmodel_10 = create_floris_model(
            wind_speeds=[10.0],
            wind_directions=[270.0],
            turbulence_intensities=[0.06],
        )
        fmodel_10.core.flow_field.air_density = 1.0
        fmodel_10.run()
        power_10 = fmodel_10.get_turbine_powers().tolist()
        
        results = {
            "power_air_density_1225": power_1225,
            "power_air_density_10": power_10,
            "ratio": power_10[0][0] / power_1225[0][0] if power_1225[0][0] > 0 else None,
        }
        save_json(results, output_dir / "floris_model_api" / "air_density_effect.json")
    except Exception as e:
        print(f"    ERROR: {e}")


# =============================================================================
# Operation Model Tests
# =============================================================================

def generate_operation_model_values(output_dir: Path) -> None:
    """Generate golden values for operation model tests."""
    print("\n=== Generating Operation Model Values ===")
    
    models = [
        "simple",
        "cosine-loss",
        "simple-derating",
        "peak-shaving",
        "mixed",
        "awc",
        "controller-dependent",
        "unified-momentum",
    ]
    
    for model_name in models:
        print(f"\n  Generating {model_name}...")
        try:
            fmodel = create_floris_model(
                wind_speeds=[10.0],
                wind_directions=[270.0],
                turbulence_intensities=[0.06],
                layout_x=[0.0],
                layout_y=[0.0],
                operation_model=model_name,
            )
            
            # Add yaw for cosine-loss model
            if model_name == "cosine-loss":
                fmodel.set(yaw_angles=[[10.0]])
            
            results = extract_results(fmodel)
            results["operation_model"] = model_name
            
            save_json(results, output_dir / "operation_models" / f"{model_name}.json")
            
        except Exception as e:
            print(f"    ERROR: {e}")
    
    # Derating test
    print("\n  Generating derating_3MW...")
    try:
        fmodel = create_floris_model(
            wind_speeds=[10.0],
            wind_directions=[270.0],
            turbulence_intensities=[0.06],
            layout_x=[0.0],
            layout_y=[0.0],
            operation_model="simple-derating",
        )
        fmodel.set(power_setpoints=[[3e6]])  # 3MW limit
        results = extract_results(fmodel)
        save_json(results, output_dir / "operation_models" / "derating_3MW.json")
    except Exception as e:
        print(f"    ERROR: {e}")

    # Peak-shaving with TI above the nrel_5MW peak_shaving_TI_threshold (0.1)
    # so the thrust limit is actually applied — the default TI=0.06 case never
    # shaves and cannot discriminate a broken implementation (Bug #35).
    print("\n  Generating peak_shaving_high_ti...")
    try:
        fmodel = create_floris_model(
            wind_speeds=[10.0],
            wind_directions=[270.0],
            turbulence_intensities=[0.15],
            layout_x=[0.0],
            layout_y=[0.0],
            operation_model="peak-shaving",
        )
        results = extract_results(fmodel)
        results["operation_model"] = "peak-shaving"
        save_json(results, output_dir / "operation_models" / "peak_shaving_high_ti.json")
    except Exception as e:
        print(f"    ERROR: {e}")


# =============================================================================
# Wind Data Tests
# =============================================================================

def generate_wind_data_values(output_dir: Path) -> None:
    """Generate golden values for wind data tests."""
    print("\n=== Generating Wind Data Values ===")
    
    # TimeSeries
    print("\n  Generating timeseries...")
    try:
        ts = TimeSeries(
            wind_speeds=np.array([8.0, 10.0, 12.0]),
            wind_directions=np.array([270.0, 270.0, 270.0]),
            turbulence_intensities=np.array([0.06, 0.06, 0.06]),
        )
        
        results = {
            "timeseries": {
                "n_findex": ts.n_findex,
                "wind_speeds": ts.wind_speeds.tolist(),
                "wind_directions": ts.wind_directions.tolist(),
                "turbulence_intensities": ts.turbulence_intensities.tolist(),
            }
        }
        save_json(results, output_dir / "wind_data" / "timeseries.json")
    except Exception as e:
        print(f"    ERROR: {e}")
    
    # WindRose
    print("\n  Generating windrose...")
    try:
        wr = WindRose(
            wind_directions=np.array([0, 90, 180, 270]),
            wind_speeds=np.array([6, 8, 10, 12, 14]),
            ti_table=np.array([
                [0.06, 0.08, 0.10, 0.12, 0.14],
                [0.06, 0.08, 0.10, 0.12, 0.14],
                [0.06, 0.08, 0.10, 0.12, 0.14],
                [0.06, 0.08, 0.10, 0.12, 0.14],
            ]),
            freq_table=np.array([
                [0.1, 0.1, 0.1, 0.1, 0.05],
                [0.1, 0.1, 0.1, 0.1, 0.05],
                [0.1, 0.1, 0.1, 0.1, 0.05],
                [0.1, 0.1, 0.1, 0.1, 0.05],
            ]),
        )
        
        results = {
            "windrose": {
                "n_wind_directions": len(wr.wind_directions),
                "n_wind_speeds": len(wr.wind_speeds),
                "wind_directions": wr.wind_directions.tolist(),
                "wind_speeds": wr.wind_speeds.tolist(),
                "freq_table_shape": list(wr.freq_table.shape),
                "freq_table_sum": float(wr.freq_table.sum()),
            }
        }
        save_json(results, output_dir / "wind_data" / "windrose.json")
    except Exception as e:
        print(f"    ERROR: {e}")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate golden values for FLORUS-FLORIS cross-validation"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=str(Path(__file__).parent / "golden_values"),
        help="Output directory for golden values (default: tests/golden_values/)"
    )
    parser.add_argument(
        "--category", "-c",
        type=str,
        choices=["wake", "api", "operation", "wind_data", "all"],
        default="all",
        help="Category of tests to generate (default: all)"
    )
    
    args = parser.parse_args()
    output_dir = Path(args.output)
    
    print("=" * 60)
    print("FLORUS-FLORIS Golden Value Generator")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Category: {args.category}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate based on category
    if args.category in ["wake", "all"]:
        generate_wake_regression_values(output_dir)
    
    if args.category in ["api", "all"]:
        generate_floris_model_api_values(output_dir)
    
    if args.category in ["operation", "all"]:
        generate_operation_model_values(output_dir)
    
    if args.category in ["wind_data", "all"]:
        generate_wind_data_values(output_dir)
    
    print("\n" + "=" * 60)
    print("Golden value generation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
