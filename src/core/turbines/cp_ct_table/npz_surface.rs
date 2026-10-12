//! Cp/Ct surface loader for the controller-dependent (TUM) operation model.
//!
//! Reads the `.npz` surface file referenced by a turbine YAML's
//! `controller_dependent_turbine_parameters.cp_ct_data_file`. The npz is a
//! zip archive of `.npy` members; Python loads it with `numpy.load` and uses
//! `cp_lut`/`ct_lut` (2D over `[tsr, pitch]`) plus the `tsr_lut`/`pitch_lut`
//! axes through `RegularGridInterpolator(..., method="cubic")`.

use crate::types::{Array1, Array2, Float};
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use zip::ZipArchive;

/// The four arrays required by Python
/// `ControllerDependentTurbine.thrust_coefficient` (line 362-365).
#[derive(Debug, Clone)]
pub struct CpCtSurface {
    /// Cp over [tsr, pitch], shape (n_tsr, n_pitch)
    pub cp_lut: Array2,
    /// Ct over [tsr, pitch], shape (n_tsr, n_pitch)
    pub ct_lut: Array2,
    /// Tip-speed-ratio axis (length n_tsr)
    pub tsr_lut: Array1,
    /// Pitch angle axis in degrees (length n_pitch)
    pub pitch_lut: Array1,
}

/// Load `cp_lut`, `ct_lut`, `tsr_lut`, `pitch_lut` from a numpy `.npz` file.
pub fn load_cp_ct_surface(path: &Path) -> anyhow::Result<CpCtSurface> {
    let file = File::open(path)
        .map_err(|e| anyhow::anyhow!("cannot open cp/ct surface {}: {}", path.display(), e))?;
    let mut archive = ZipArchive::new(file)
        .map_err(|e| anyhow::anyhow!("cannot parse npz {}: {}", path.display(), e))?;

    let mut arrays: HashMap<String, NpyArray> = HashMap::new();
    for i in 0..archive.len() {
        let mut member = archive
            .by_index(i)
            .map_err(|e| anyhow::anyhow!("cannot read npz member: {}", e))?;
        let name = member.name().to_string();
        // Strip the ".npy" suffix numpy appends to array names inside npz
        let key = name.strip_suffix(".npy").unwrap_or(&name).to_string();
        let mut contents = Vec::new();
        member
            .read_to_end(&mut contents)
            .map_err(|e| anyhow::anyhow!("cannot read {}: {}", name, e))?;
        let arr = parse_npy(&contents)
            .map_err(|e| anyhow::anyhow!("cannot parse npy member {}: {}", name, e))?;
        arrays.insert(key, arr);
    }

    let cp = arrays
        .remove("cp_lut")
        .ok_or_else(|| anyhow::anyhow!("npz {} is missing cp_lut", path.display()))?;
    let ct = arrays
        .remove("ct_lut")
        .ok_or_else(|| anyhow::anyhow!("npz {} is missing ct_lut", path.display()))?;
    let tsr = arrays
        .remove("tsr_lut")
        .ok_or_else(|| anyhow::anyhow!("npz {} is missing tsr_lut", path.display()))?;
    let pitch = arrays
        .remove("pitch_lut")
        .ok_or_else(|| anyhow::anyhow!("npz {} is missing pitch_lut", path.display()))?;

    let cp_lut = cp.into_array2()?;
    let ct_lut = ct.into_array2()?;
    let tsr_lut = tsr.into_array1()?;
    let pitch_lut = pitch.into_array1()?;

    // Sanity: the 2D luts are indexed [tsr, pitch] (numpy save order)
    if cp_lut.shape() != ct_lut.shape()
        || cp_lut.shape()[0] != tsr_lut.len()
        || cp_lut.shape()[1] != pitch_lut.len()
    {
        anyhow::bail!(
            "cp/ct surface shape mismatch: cp {:?}, ct {:?}, tsr {}, pitch {}",
            cp_lut.shape(),
            ct_lut.shape(),
            tsr_lut.len(),
            pitch_lut.len()
        );
    }

    Ok(CpCtSurface {
        cp_lut,
        ct_lut,
        tsr_lut,
        pitch_lut,
    })
}

/// A minimal parsed .npy payload: shape + float64 data (any dtype order).
struct NpyArray {
    shape: Vec<usize>,
    fortran_order: bool,
    data: Vec<f64>,
}

impl NpyArray {
    fn into_array1(self) -> anyhow::Result<Array1> {
        if self.shape.len() != 1 {
            anyhow::bail!("expected 1D array, got shape {:?}", self.shape);
        }
        Ok(Array1::from_vec(self.data))
    }

    fn into_array2(self) -> anyhow::Result<Array2> {
        if self.shape.len() != 2 {
            anyhow::bail!("expected 2D array, got shape {:?}", self.shape);
        }
        let (r, c) = (self.shape[0], self.shape[1]);
        if self.data.len() != r * c {
            anyhow::bail!("npy data length mismatch");
        }
        // np.save stores C-order by default; honor the fortran_order flag by
        // transposing the row-major interpretation when set.
        let mut out = Array2::zeros((r, c));
        for i in 0..r {
            for j in 0..c {
                let k = if self.fortran_order {
                    j * r + i
                } else {
                    i * c + j
                };
                out[[i, j]] = self.data[k];
            }
        }
        Ok(out)
    }
}

/// Parse a `.npy` v1/v2 file with float64 (or int64/int32/float32) contents.
fn parse_npy(data: &[u8]) -> anyhow::Result<NpyArray> {
    if data.len() < 10 || &data[0..6] != b"\x93NUMPY" {
        anyhow::bail!("not a valid NPY file");
    }
    let (major, _minor) = (data[6], data[7]);
    let header_len_pos = if major == 1 { 8 } else { 10 };
    let header_len = u16::from_le_bytes([data[header_len_pos], data[header_len_pos + 1]]) as usize;
    let header_start = header_len_pos + 2;
    let header = std::str::from_utf8(&data[header_start..header_start + header_len])
        .map_err(|e| anyhow::anyhow!("invalid npy header: {}", e))?;

    // Parse 'descr', 'fortran_order', 'shape' out of the python-dict header.
    let descr = header
        .split("'descr':")
        .nth(1)
        .and_then(|s| s.split('\'').nth(1))
        .ok_or_else(|| anyhow::anyhow!("npy header missing descr"))?
        .to_string();
    let fortran_order = header
        .split("'fortran_order':")
        .nth(1)
        .map(|s| s.trim_start().starts_with("True"))
        .unwrap_or(false);
    let shape_str = header
        .split("'shape':")
        .nth(1)
        .and_then(|s| s.split(')').next())
        .ok_or_else(|| anyhow::anyhow!("npy header missing shape"))?;
    let shape: Vec<usize> = shape_str
        .trim_start_matches(" (")
        .trim_end()
        .split(',')
        .filter_map(|t| t.trim().parse::<usize>().ok())
        .collect();

    let body = &data[header_start + header_len..];
    let n: usize = shape.iter().product();
    let mut out = Vec::with_capacity(n);
    match descr.as_str() {
        "<f8" | "|f8" => {
            if body.len() < n * 8 {
                anyhow::bail!("npy body too short");
            }
            for k in 0..n {
                let b: [u8; 8] = body[k * 8..k * 8 + 8]
                    .try_into()
                    .map_err(|_| anyhow::anyhow!("slice error"))?;
                out.push(f64::from_le_bytes(b));
            }
        }
        "<i8" | "|i8" => {
            if body.len() < n * 8 {
                anyhow::bail!("npy body too short");
            }
            for k in 0..n {
                let b: [u8; 8] = body[k * 8..k * 8 + 8]
                    .try_into()
                    .map_err(|_| anyhow::anyhow!("slice error"))?;
                out.push(i64::from_le_bytes(b) as Float);
            }
        }
        "<i4" | "|i4" => {
            if body.len() < n * 4 {
                anyhow::bail!("npy body too short");
            }
            for k in 0..n {
                let b: [u8; 4] = body[k * 4..k * 4 + 4]
                    .try_into()
                    .map_err(|_| anyhow::anyhow!("slice error"))?;
                out.push(i32::from_le_bytes(b) as Float);
            }
        }
        "<f4" | "|f4" => {
            if body.len() < n * 4 {
                anyhow::bail!("npy body too short");
            }
            for k in 0..n {
                let b: [u8; 4] = body[k * 4..k * 4 + 4]
                    .try_into()
                    .map_err(|_| anyhow::anyhow!("slice error"))?;
                out.push(f32::from_le_bytes(b) as Float);
            }
        }
        other => anyhow::bail!("unsupported npy dtype: {}", other),
    }

    Ok(NpyArray {
        shape,
        fortran_order,
        data: out,
    })
}
