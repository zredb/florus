#!/usr/bin/env python3
"""
Convert turbopark_lookup_table.mat to Rust lookup table data.
This script is called by build.rs during cargo build.
"""

import sys
import os
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Convert .mat file to Rust lookup table')
    parser.add_argument('--out-dir', type=str, help='Output directory for generated file')
    args = parser.parse_args()

    mat_file = 'src/core/wake/wake_velocity/turbopark_lookup_table.mat'

    if args.out_dir:
        output_file = os.path.join(args.out_dir, 'turbopark_lookup.rs')
    else:
        output_file = 'src/turbopark_lookup.rs'

    # Generate dist and radius_down arrays (always the same)
    dist = np.linspace(0, 10, 100)
    radius_down = np.linspace(0, 20, 100)

    # Try to load overlap_gauss from .mat file
    overlap_gauss = None
    if os.path.exists(mat_file):
        try:
            import scipy.io
            print(f"Loading {mat_file}...", file=sys.stderr)
            mat = scipy.io.loadmat(mat_file)
            keys = [k for k in mat.keys() if not k.startswith('__')]
            print(f"Keys found: {keys}", file=sys.stderr)

            # Find the overlap_gauss array
            for k in keys:
                data = mat[k]
                if hasattr(data, 'shape') and len(data.shape) == 2 and data.shape[0] == 100 and data.shape[1] == 100:
                    overlap_gauss = np.asarray(data)
                    print(f"Found 100x100 overlap table: {k}", file=sys.stderr)
                    break
        except Exception as e:
            print(f"Error loading .mat file: {e}", file=sys.stderr)
    else:
        print(f"Mat file not found: {mat_file}", file=sys.stderr)

    # If overlap_gauss not found, generate it (fallback)
    if overlap_gauss is None:
        print("Generating overlap_gauss from Simpson integral...", file=sys.stderr)
        overlap_gauss = np.zeros((100, 100))

        for i, d in enumerate(dist):
            for j, rd in enumerate(radius_down):
                if rd <= 0:
                    overlap_gauss[i, j] = np.exp(-d**2 / 2)
                elif d > 10:
                    overlap_gauss[i, j] = np.exp(-d**2 / 2)
                else:
                    n_points = 100
                    dr = rd / (n_points - 1)
                    integral = 0.0
                    for k in range(n_points):
                        r = k * dr
                        decay_val = -(r**2 + d**2 - 2 * d * r)
                        if decay_val < 0:
                            decay_val = 0
                        decay = np.exp(decay_val / 2)
                        integrand = r * decay
                        weight = 1 if k == 0 or k == n_points - 1 else (2 if k % 2 == 0 else 4)
                        integral += weight * integrand
                    integral *= dr / 3.0
                    area = np.pi * rd**2
                    overlap_gauss[i, j] = max(0.0, min(1.0, integral / area))

    # Convert to Python lists for Rust output
    dist_list = [float(v) for v in dist.tolist()]
    radius_down_list = [float(v) for v in radius_down.tolist()]
    overlap_list = [[float(v) for v in row] for row in overlap_gauss.tolist()]

    # Generate Rust code
    with open(output_file, 'w') as f:
        f.write("/// Auto-generated lookup table for TurbOPark wake model.\n")
        f.write("/// Generated from turbopark_lookup_table.mat by build.rs.\n\n")
        f.write("pub const OVERLAP_DIST: [f64; 100] = [\n    ")
        f.write(", ".join([f"{v:.6}" for v in dist_list]))
        f.write("\n];\n\n")
        f.write("pub const OVERLAP_RADIUS_DOWN: [f64; 100] = [\n    ")
        f.write(", ".join([f"{v:.6}" for v in radius_down_list]))
        f.write("\n];\n\n")
        f.write("pub const OVERLAP_GAUSS: [[f64; 100]; 100] = [\n")
        for row in overlap_list:
            f.write("    [")
            f.write(", ".join([f"{v:.10}" for v in row]))
            f.write("],\n")
        f.write("];\n")

    print(f"Generated {output_file} with {len(overlap_list)}x{len(overlap_list[0])} values", file=sys.stderr)

if __name__ == '__main__':
    main()
