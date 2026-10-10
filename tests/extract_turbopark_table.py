"""Extract TurbOPark overlap lookup table from MATLAB .mat file to Rust source.

Uses OnceLock + Vec for heap-allocated storage to avoid 16MB stack overflow.
"""
import scipy.io
import numpy as np
from pathlib import Path

mat_file = Path(__file__).parent.parent / "floris" / "floris" / "core" / "wake_velocity" / "turbopark_lookup_table.mat"
data = scipy.io.loadmat(str(mat_file))

dist = data['overlap_lookup_table'][0][0][0][0]
radius_down = data['overlap_lookup_table'][0][0][1][0]
overlap_gauss = data['overlap_lookup_table'][0][0][2]

print(f"dist shape: {dist.shape}, range: [{dist[0]}, {dist[-1]}]")
print(f"radius_down shape: {radius_down.shape}, range: [{radius_down[0]}, {radius_down[-1]}]")
print(f"overlap_gauss shape: {overlap_gauss.shape}")
print(f"overlap_gauss range: [{np.nanmin(overlap_gauss)}, {np.nanmax(overlap_gauss)}]")

lines = []
lines.append("// Auto-generated from turbopark_lookup_table.mat - DO NOT EDIT")
lines.append("use std::sync::OnceLock;")
lines.append("use crate::types::Float;")
lines.append("")
lines.append(f"pub const DIST_N: usize = {len(dist)};")
lines.append(f"pub const RADIUS_DOWN_N: usize = {len(radius_down)};")
lines.append("")
lines.append("static DIST: OnceLock<Vec<Float>> = OnceLock::new();")
lines.append("static RADIUS_DOWN: OnceLock<Vec<Float>> = OnceLock::new();")
lines.append("static OVERLAP_GAUSS: OnceLock<Vec<Float>> = OnceLock::new();")
lines.append("")
lines.append("pub fn dist() -> &'static [Float] {")
lines.append("    DIST.get_or_init(|| vec![")
for v in dist:
    lines.append(f"        {v:.15e},")
lines.append("    ])")
lines.append("}")
lines.append("")
lines.append("pub fn radius_down() -> &'static [Float] {")
lines.append("    RADIUS_DOWN.get_or_init(|| vec![")
for v in radius_down:
    lines.append(f"        {v:.15e},")
lines.append("    ])")
lines.append("}")
lines.append("")
lines.append("/// Returns flattened overlap_gauss data (DIST_N x RADIUS_DOWN_N, row-major)")
lines.append("pub fn overlap_gauss() -> &'static [Float] {")
lines.append("    OVERLAP_GAUSS.get_or_init(|| {")
lines.append("        let mut data = Vec::with_capacity(DIST_N * RADIUS_DOWN_N);")
for row in overlap_gauss:
    for v in row:
        if np.isnan(v):
            lines.append("        data.push(f64::NAN);")
        else:
            lines.append(f"        data.push({v:.15e});")
lines.append("        data")
lines.append("    })")
lines.append("}")
lines.append("")
lines.append("/// Get overlap value at (i, j) index")
lines.append("#[inline]")
lines.append("pub fn overlap_at(i: usize, j: usize) -> Float {")
lines.append("    overlap_gauss()[i * RADIUS_DOWN_N + j]")
lines.append("}")

output = Path(__file__).parent.parent / "src" / "core" / "wake" / "wake_velocity" / "turbopark_lookup_matlab.rs"
output.write_text("\n".join(lines))
print(f"\nGenerated: {output}")
print(f"Total lines: {len(lines)}")
print(f"Data points: {len(dist) * len(radius_down)} ({len(dist)} x {len(radius_down)})")
