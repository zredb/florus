use serde::{Deserialize, Serialize};
use std::fs::File;
use std::path::{self, Path};
use std::sync::Mutex;
use std::{collections::HashMap, sync::OnceLock};

use crate::core::turbines::TurbineTypeError;

use super::turbine_type::TurbineType;

#[derive(Debug, Clone)]
pub struct TurbineLibrary {
    pub turbines: HashMap<String, TurbineType>,
    pub external_turbine_paths: Vec<std::path::PathBuf>,
}

static TURBINE_LIBRARY: OnceLock<Mutex<TurbineLibrary>> = OnceLock::new();

impl TurbineLibrary {
    /// 创建一个新的空的 TurbineLibrary 实例
    fn new() -> Self {
        Self {
            turbines: HashMap::new(),
            external_turbine_paths: Vec::new(),
        }
    }
    /// 获取全局唯一的 TurbineLibrary 实例
    pub fn instance() -> std::sync::MutexGuard<'static, TurbineLibrary> {
        TURBINE_LIBRARY
            .get_or_init(|| Mutex::new(TurbineLibrary::new()))
            .lock()
            .expect("Failed to acquire lock on TurbineLibrary")
    }
    pub fn set_external_turbine_paths<P: AsRef<Path>>(&mut self, paths: Vec<P>) {
        self.external_turbine_paths = paths
            .into_iter()
            .map(|p| p.as_ref().to_path_buf())
            .collect();
    }
    pub fn add_external_turbine_path<P: AsRef<Path>>(&mut self, path: P) {
        self.external_turbine_paths
            .push(path.as_ref().to_path_buf());
    }

    /// 初始化并加载所有可用的风机（如果尚未初始化）
    pub fn init_if_needed() -> Result<(), TurbineTypeError> {
        let mut library = Self::instance();

        // 如果库为空，则加载所有风机
        if library.turbines.is_empty() {
            library.load_all_available_turbines()?;
        }

        Ok(())
    }

    /// 添加一个风机到库中
    fn add_turbine(&mut self, turbine_type: TurbineType) {
        self.turbines
            .insert(turbine_type.name.clone(), turbine_type);
    }

    /// 从文件路径加载风机并添加到库中
    fn load_turbine_from_file(&mut self, path: &str) -> Result<(), TurbineTypeError> {
        let turbine = TurbineType::load_turbine_type(path)?;
        self.add_turbine(turbine);
        Ok(())
    }
    fn load_internal_turbines(&mut self) -> Result<(), TurbineTypeError> {
        let internal_turbine_files = [
            "turbine_library/nrel_5MW.yaml",
            "turbine_library/iea_10MW.yaml",
            "turbine_library/iea_15MW.yaml",
            "turbine_library/iea_15MW_floating_multi_dim_cp_ct.yaml",
            "turbine_library/iea_15MW_multi_dim_cp_ct.yaml",
            "turbine_library/iea_22MW.yaml",
        ];

        for file in &internal_turbine_files {
            if Path::new(file).exists() {
                match self.load_turbine_from_file(file) {
                    Ok(_) => {},  // Silent success
                    Err(e) => eprintln!("Failed to load {}: {}", file, e),
                }
            } else {
                eprintln!("File does not exist: {}", file);
            }
        }
        Ok(())
    }
    fn load_external_turbines(&mut self) -> Result<(), TurbineTypeError> {
        let mut external_files_to_load = Vec::new();
        for path in &self.external_turbine_paths {
            if Path::new(path).exists() {
                for file in std::fs::read_dir(path)? {
                    let file = file?;
                    if file.path().extension().and_then(|s| s.to_str()) == Some("yaml") {
                        external_files_to_load.push(file.path().to_string_lossy().to_string());
                    }
                }
            } else {
                eprintln!("External path does not exist: {:?}", path);
            }
        }

        // 加载收集到的外部风机文件
        for file_path in external_files_to_load {
            match self.load_turbine_from_file(&file_path) {
                Ok(_) => {},  // Silent success
                Err(e) => eprintln!("Failed to load external turbine {}: {}", file_path, e),
            }
        }

        Ok(())
    }
    /// 从 turbine_library 目录加载所有可用的风机
    pub fn load_all_available_turbines(&mut self) -> Result<(), TurbineTypeError> {
        self.load_internal_turbines()?;
        self.load_external_turbines()?;
        Ok(())
    }

    /// 获取指定名称的风机
    pub fn get_turbine(name: &str) -> Option<TurbineType> {
        // 首先尝试获取风机
        let library = Self::instance();
        if let Some(turbine) = library.turbines.get(name) {
            return Some(turbine.clone());
        }
        None
    }

    /// 获取所有已加载的风机类型
    pub fn get_loaded_turbines() -> Vec<String> {
        let library = Self::instance();
        library.turbines.keys().cloned().collect()
    }

    /// 清除所有缓存的风机（用于测试或重新加载）
    pub fn clear_cache() {
        let mut library = Self::instance();
        library.turbines.clear();
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_new_creates_empty_library() {
        let library = TurbineLibrary::new();
        assert!(library.turbines.is_empty());
        assert!(library.external_turbine_paths.is_empty());
    }

    #[test]
    fn test_set_and_add_external_paths() {
        let mut library = TurbineLibrary::new();

        // 测试设置外部路径
        let paths = vec!["path1".to_string(), "path2".to_string()];
        library.set_external_turbine_paths(paths.clone());
        assert_eq!(library.external_turbine_paths, paths);

        // 测试添加单个外部路径
        library.add_external_turbine_path("path3".to_string());
        assert_eq!(
            library.external_turbine_paths,
            vec![
                "path1".to_string(),
                "path2".to_string(),
                "path3".to_string()
            ]
        );
    }

    #[test]
    fn test_load_turbine_from_nonexistent_file() {
        let mut library = TurbineLibrary::new();
        let result = library.load_turbine_from_file("nonexistent_file.yaml");

        // 验证加载失败的情况
        assert!(result.is_err());
    }

    #[test]
    fn test_load_internal_turbines() {
        let mut library = TurbineLibrary::new();
        let initial_count = library.turbines.len();

        // 尝试加载内部风机（即使某些文件不存在，也应该不崩溃）
        let result = library.load_internal_turbines();

        // 即使某些文件不存在，函数也应成功返回
        assert!(result.is_ok());

        // 不管怎样，至少应该尝试处理内部文件列表
        // 验证没有出现panic
    }

    #[test]
    fn test_get_loaded_turbines() {
        let loaded_turbines = TurbineLibrary::get_loaded_turbines();

        assert!(loaded_turbines.contains(&"turbine1".to_string()));
        assert!(loaded_turbines.contains(&"turbine2".to_string()));
        assert_eq!(loaded_turbines.len(), 2);
    }

    #[test]
    fn test_clear_cache() {
        // 使用静态方法清除缓存
        TurbineLibrary::clear_cache();
        // 重新获取实例验证缓存已被清除
        let library_after_clear = TurbineLibrary::instance();
        assert!(library_after_clear.turbines.is_empty());
    }

    #[test]
    fn test_init_if_needed_first_time() {
        // 先清空缓存以模拟首次初始化
        TurbineLibrary::clear_cache();

        let result = TurbineLibrary::init_if_needed();

        assert!(result.is_ok());

        // 检查是否执行了加载操作
        let library = TurbineLibrary::instance();
        // 这里会尝试加载内部风机文件，即使它们不存在也不会出错
    }

    #[test]
    fn test_init_if_needed_subsequent_calls() {
        // 先初始化一次
        TurbineLibrary::clear_cache(); // 确保开始时是空的
        TurbineLibrary::init_if_needed().unwrap();

        // 再次调用，应该不会重复加载
        let result = TurbineLibrary::init_if_needed();

        assert!(result.is_ok());
    }
}
