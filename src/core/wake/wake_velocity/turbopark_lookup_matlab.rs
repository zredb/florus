// Auto-generated from turbopark_lookup_table.mat via generate_turbopark_bin.py
// Uses include_bytes! + OnceLock for efficient heap-allocated storage.
use std::sync::OnceLock;

use crate::types::Float;

pub const DIST_N: usize = 1001;
pub const RADIUS_DOWN_N: usize = 2001;

static STATE: OnceLock<TableData> = OnceLock::new();

struct TableData {
    dist: Vec<Float>,
    radius_down: Vec<Float>,
    overlap_flat: Vec<Float>,
}

fn load_blob() -> &'static TableData {
    STATE.get_or_init(|| {
        const BLOB: &[u8] = include_bytes!("../../../../data/turbopark_overlap.bin");

        let hdr_len = 16; // 2 x u64
        let dist_n = u64::from_le_bytes(BLOB[0..8].try_into().unwrap()) as usize;
        let radius_down_n = u64::from_le_bytes(BLOB[8..16].try_into().unwrap()) as usize;
        assert_eq!(dist_n, DIST_N);
        assert_eq!(radius_down_n, RADIUS_DOWN_N);

        let dist_bytes = &BLOB[hdr_len..hdr_len + dist_n * 8];
        let dist: Vec<Float> = dist_bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| Float::from_le_bytes(*c))
            .collect();

        let rd_off = hdr_len + dist_n * 8;
        let radius_down_bytes = &BLOB[rd_off..rd_off + radius_down_n * 8];
        let radius_down: Vec<Float> = radius_down_bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| Float::from_le_bytes(*c))
            .collect();

        let ov_off = rd_off + radius_down_n * 8;
        let overlap_flat: Vec<Float> = BLOB[ov_off..]
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| Float::from_le_bytes(*c))
            .collect();
        assert_eq!(overlap_flat.len(), dist_n * radius_down_n);

        TableData {
            dist,
            radius_down,
            overlap_flat,
        }
    })
}

pub fn dist() -> &'static [Float] {
    &load_blob().dist
}

pub fn radius_down() -> &'static [Float] {
    &load_blob().radius_down
}

#[inline]
pub fn overlap_at(i: usize, j: usize) -> Float {
    load_blob().overlap_flat[i * RADIUS_DOWN_N + j]
}
