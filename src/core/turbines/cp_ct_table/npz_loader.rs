pub struct NpzDataLoader;
use crate::types::{Array2, DynArray};

use super::{CpCtError, DimensionRange, MultiDimTable};
use ndarray::ArrayD;
use std::{fs::File, io::Read, path::Path};
use zip::ZipArchive;

impl NpzDataLoader {
    /// 加载.npz文件数据
    pub fn load_npz_file(file_path: &Path) -> Result<MultiDimTable, CpCtError> {
        if !file_path.exists() {
            return Err(CpCtError::FileNotFound(
                file_path.to_string_lossy().to_string(),
            ));
        }

        if file_path.extension().and_then(|ext| ext.to_str()) != Some("npz") {
            return Err(CpCtError::InvalidNpzFormat(
                "文件必须是.npz格式".to_string(),
            ));
        }

        let file = File::open(file_path)?;
        let mut archive = ZipArchive::new(file)
            .map_err(|e| CpCtError::FileFormat(format!("无法解析NPZ文件: {}", e)))?;

        let mut cp_data: Option<DynArray> = None;
        let mut ct_data: Option<DynArray> = None;
        let mut dimensions: Vec<DimensionRange> = vec![];

        for i in 0..archive.len() {
            let mut file = archive.by_index(i)
                .map_err(|e| CpCtError::FileFormat(format!("无法读取NPZ文件内容: {}", e)))?;
            
            let name = file.name().to_string();
            println!("Found array: {}", name);

            let mut contents = Vec::new();
            file.read_to_end(&mut contents)
                .map_err(|e| CpCtError::FileFormat(format!("无法读取 {}: {}", name, e)))?;

            // Parse .npy format
            let arr = parse_npy(&contents)
                .map_err(|e| CpCtError::FileFormat(format!("无法解析 {}: {}", name, e)))?;

            match name.as_str() {
                "cp_values.npy" => {
                    cp_data = Some(arr.into_dyn());
                    println!("Loaded cp_values with shape: {:?}", cp_data.as_ref().map(|a| a.shape()));
                }
                "ct_values.npy" => {
                    ct_data = Some(arr.into_dyn());
                    println!("Loaded ct_values with shape: {:?}", ct_data.as_ref().map(|a| a.shape()));
                }
                _ => {
                    println!("Found unused file: {}", name);
                }
            }
        }

        let cp_data = cp_data.ok_or_else(|| CpCtError::LoadFailed("Missing cp_values array".to_string()))?;
        let ct_data = ct_data.ok_or_else(|| CpCtError::LoadFailed("Missing ct_values array".to_string()))?;

        Ok(MultiDimTable::new(dimensions, cp_data, ct_data)?)
    }
}

fn parse_npy(data: &[u8]) -> Result<ndarray::Array2<f64>, CpCtError> {
    // Simple .npy parser for fortran-order 2D arrays
    if data.len() < 10 {
        return Err(CpCtError::FileFormat("NPY file too short".to_string()));
    }

    // Check magic bytes
    if &data[0..6] != b"\x93NUMPY" {
        return Err(CpCtError::FileFormat("Not a valid NPY file".to_string()));
    }

    let version = (data[6], data[7]);
    let header_start = if version == (1, 0) { 8 } else { 12 };

    // Parse header to get shape and dtype
    let mut pos = header_start;
    while pos < data.len() && data[pos] != b'\n' {
        pos += 1;
    }
    let header = std::str::from_utf8(&data[header_start..pos])
        .map_err(|e| CpCtError::FileFormat(format!("Invalid header: {}", e)))?;

    // Extract shape from header - look for '(', ')' and commas
    let shape_start = header.find('(').unwrap_or(header.len());
    let shape_end = header.find(')').unwrap_or(header.len());
    let shape_str = &header[shape_start+1..shape_end];
    
    let shape: Vec<usize> = shape_str
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    if shape.len() != 2 {
        return Err(CpCtError::FileFormat(format!("Expected 2D array, got {}D", shape.len())));
    }

    // Data starts at position pos + 1 (skip newline), aligned to 64 bytes
    let data_start = (pos + 1 + 63) & !63;

    if data.len() < data_start + shape[0] * shape[1] * 8 {
        return Err(CpCtError::FileFormat("NPY file data truncated".to_string()));
    }

    let data_slice = &data[data_start..];
    let mut arr = ndarray::Array2::<f64>::zeros((shape[0], shape[1]));

    // NPY stores in C order (row-major), but we want to keep as-is
    for i in 0..shape[0] {
        for j in 0..shape[1] {
            let bytes: [u8; 8] = data_slice[(i * shape[1] + j) * 8..(i * shape[1] + j + 1) * 8].try_into()
                .map_err(|_| CpCtError::FileFormat("Failed to read float".to_string()))?;
            arr[[i, j]] = f64::from_le_bytes(bytes);
        }
    }

    Ok(arr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_load_npy_file_success() {
        let mut test_file = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_file.push("turbine_library/demo_cp_ct_surfaces/iea_10MW_demo_cp_ct_surface.npz");

        let result = NpzDataLoader::load_npz_file(&test_file);

        assert!(
            result.is_ok(),
            "Expected successful loading of .npz file, but got error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_load_nonexistent_file() {
        let nonexistent_file = Path::new("nonexistent_file.npz");
        let result = NpzDataLoader::load_npz_file(nonexistent_file);

        assert!(result.is_err());
        match result.err().unwrap() {
            CpCtError::FileNotFound(_) => {}
            _ => panic!("Expected FileNotFound error"),
        }
    }

    #[test]
    fn test_load_invalid_extension() {
        let invalid_file = Path::new("invalid_file.txt");
        let result = NpzDataLoader::load_npz_file(invalid_file);

        assert!(result.is_err());
        match result.err().unwrap() {
            CpCtError::InvalidNpzFormat(_) => {}
            _ => panic!("Expected InvalidNpzFormat error"),
        }
    }
}
