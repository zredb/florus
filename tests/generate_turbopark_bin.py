"""Serialize TurbOPark overlap lookup table as a compact binary blob.

Format (all little-endian f64):
  header: 8 bytes DIST_N (u64), 8 bytes RADIUS_DOWN_N (u64)
  dist:    DIST_N * 8 bytes
  radius_down: RADIUS_DOWN_N * 8 bytes
  overlap_gauss: DIST_N * RADIUS_DOWN_N * 8 bytes (row-major)

Total: 16 + (1001 + 2001 + 1001*2001) * 8 = 16 + 16,032,008 ~16MB
"""
import struct
import numpy as np
from pathlib import Path

mat_file = Path(__file__).parent.parent / "floris" / "floris" / "core" / "wake_velocity" / "turbopark_lookup_table.mat"
import scipy.io
data = scipy.io.loadmat(str(mat_file))

dist = data['overlap_lookup_table'][0][0][0][0].astype(np.float64)
radius_down = data['overlap_lookup_table'][0][0][1][0].astype(np.float64)
overlap_gauss = data['overlap_lookup_table'][0][0][2].astype(np.float64)

out_dir = Path(__file__).parent.parent / "data"
out_dir.mkdir(exist_ok=True)
out_file = out_dir / "turbopark_overlap.bin"

with open(out_file, "wb") as f:
    f.write(struct.pack("<QQ", len(dist), len(radius_down)))
    f.write(dist.tobytes())
    f.write(radius_down.tobytes())
    f.write(overlap_gauss.tobytes())

print(f"Written: {out_file} ({out_file.stat().st_size} bytes)")
print(f"dist: {len(dist)} [{dist[0]}, {dist[-1]}]")
print(f"radius_down: {len(radius_down)} [{radius_down[0]}, {radius_down[-1]}]")
print(f"overlap: {overlap_gauss.shape} [{np.nanmin(overlap_gauss)}, {np.nanmax(overlap_gauss)}]")
