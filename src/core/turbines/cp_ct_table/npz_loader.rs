pub struct NpzDataLoader;
use crate::types::{Array2, DynArray};

use super::{CpCtError, DimensionRange, MultiDimTable};
use ndarray::ArrayD;
use ndarray_npy::NpzReader;
use std::{fs::File, path::Path};

impl NpzDataLoader {
    /// 加载.npy文件数据
    pub fn load_npz_file(file_path: &Path) -> Result<MultiDimTable, CpCtError> {
        // 检查文件是否存在
        if !file_path.exists() {
            return Err(CpCtError::FileNotFound(
                file_path.to_string_lossy().to_string(),
            ));
        }

        // 检查文件扩展名
        if file_path.extension().and_then(|ext| ext.to_str()) != Some("npz") {
            return Err(CpCtError::InvalidNpzFormat(
                "文件必须是.npz格式".to_string(),
            ));
        }
        let file = File::open(file_path)?;
        let mut npz = NpzReader::new(file)
            .map_err(|e| CpCtError::FileFormat(format!("无法解析NPZ文件: {}", e)))?;

        // 列出npz文件中的所有数组名称
        for name in npz
            .names()
            .map_err(|e| CpCtError::FileFormat(format!("无法读取NPZ文件内容: {}", e)))?
        {
            println!("Found array: {}", name);
        }

        // 读取所有数组并合并成MultiDimTable
        let mut cp_data: Option<DynArray> = None;
        let mut ct_data: Option<DynArray> = None;
        let mut dimensions: Vec<DimensionRange> = vec![];

        for name in npz
            .names()
            .map_err(|e| CpCtError::FileFormat(format!("无法读取NPZ文件内容: {}", e)))?
        {
            match name.as_str() {
                "cp_values" => {
                    let arr = npz
                        .by_name::<f64>("cp_values")
                        .map_err(|e| CpCtError::FileFormat(format!("无法读取cp_values: {}", e)))?;
                    let shape = arr.dim();
                    cp_data = Some(arr.into_dyn());
                    println!("Loaded cp_values with shape: {:?}", shape);
                }
                "ct_values" => {
                    let arr = npz
                        .by_name::<f64>("ct_values")
                        .map_err(|e| CpCtError::FileFormat(format!("无法读取ct_values: {}", e)))?;
                    let shape = arr.dim();
                    ct_data = Some(arr.into_dyn());
                    println!("Loaded ct_values with shape: {:?}", shape);
                }
                _ => {
                    println!("Found unused array: {}", name);
                }
            }
        }

        // 如果已经读取到数据，则使用它们，否则返回错误
        let cp_data =
            cp_data.ok_or_else(|| CpCtError::LoadFailed("Missing cp_values array".to_string()))?;
        let ct_data =
            ct_data.ok_or_else(|| CpCtError::LoadFailed("Missing ct_values array".to_string()))?;

        return Ok(MultiDimTable::new(dimensions, cp_data, ct_data)?);

        Ok(MultiDimTable::new(dimensions, cp_data, ct_data)?)
    }
}

#[cfg(test)]
mod tests {

    use crate::core::turbines::cp_ct_table::TableConditions;

    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_load_npy_file_success() {
        // 构建测试文件路径
        let mut test_file = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        test_file.push("src/turbine_library/demo_cp_ct_surfaces/iea_10MW_demo_cp_ct_surface.npz");

        // 测试加载功能
        let result = NpzDataLoader::load_npz_file(&test_file);

        // 验证结果
        assert!(
            result.is_ok(),
            "Expected successful loading of .npz file, but got error: {:?}",
            result.err()
        );

        let multi_dim_table = result.unwrap();

        // 验证加载的数据是否有效
        let conditions = TableConditions {
            wind_speed: 10.0,
            ..Default::default()
        };

        // 尝试获取Cp和Ct值以验证数据完整性
        let cp_result = multi_dim_table.get_cp(&conditions);
        let ct_result = multi_dim_table.get_ct(&conditions);

        assert!(
            cp_result.is_ok(),
            "Should be able to get Cp value: {:?}",
            cp_result.err()
        );
        assert!(
            ct_result.is_ok(),
            "Should be able to get Ct value: {:?}",
            ct_result.err()
        );

        let cp_value = cp_result.unwrap();
        let ct_value = ct_result.unwrap();

        // 验证Cp和Ct值在合理范围内
        assert!(
            cp_value >= 0.0 && cp_value <= 1.0,
            "Cp value should be in [0, 1] range, got: {}",
            cp_value
        );
        assert!(
            ct_value >= 0.0 && ct_value <= 1.0,
            "Ct value should be in [0, 1] range, got: {}",
            ct_value
        );
    }

    #[test]
    fn test_load_nonexistent_file() {
        let nonexistent_file = Path::new("nonexistent_file.npz");
        let result = NpzDataLoader::load_npz_file(nonexistent_file);

        assert!(result.is_err());
        match result.err().unwrap() {
            CpCtError::FileNotFound(_) => {} // 期望的错误类型
            _ => panic!("Expected FileNotFound error"),
        }
    }

    #[test]
    fn test_load_invalid_extension() {
        let invalid_file = Path::new("invalid_file.txt");
        let result = NpzDataLoader::load_npz_file(invalid_file);

        assert!(result.is_err());
        match result.err().unwrap() {
            CpCtError::InvalidNpzFormat(_) => {} // 期望的错误类型
            _ => panic!("Expected InvalidNpzFormat error"),
        }
    }
}
// ... existing code ...
