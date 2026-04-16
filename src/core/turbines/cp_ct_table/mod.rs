// cp_ct_table.rs
// 统一的Cp/Ct表枚举实现

mod csv_loader;
mod npz_loader;

pub use csv_loader::CsvDataLoader;
use std::collections::HashMap;

use crate::types::{Array1, DynArray, Float};
use derive_builder::Builder;
use serde::{Deserialize, Serialize};
use thiserror::Error;
// ============================================================================
// 错误类型
// ============================================================================

#[derive(Error, Debug)]
pub enum CpCtError {
    #[error("维度不匹配: {0}")]
    DimensionMismatch(String),

    #[error("插值失败: {0}")]
    InterpolationFailed(String),

    #[error("数据验证失败: {0}")]
    ValidationFailed(String),

    #[error("文件格式错误: {0}")]
    FileFormat(String),

    #[error("文件未找到: {0}")]
    FileNotFound(String),

    #[error("不支持的操作: {0}")]
    UnsupportedOperation(String),

    #[error("CSV格式错误: {0}")]
    InvalidCsvFormat(String),

    #[error("Npz格式错误: {0}")]
    InvalidNpzFormat(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("IO error: {0}")]
    Csv(#[from] csv::Error),

    #[error("加载失败: {0}")]
    LoadFailed(String),
}

// ============================================================================
// 统一的Cp/Ct表枚举
// ============================================================================

/// 统一的Cp/Ct表枚举，支持一维和多维表
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum CpCtTable {
    /// 一维表：风速 -> (Cp, Ct)
    OneDimensional(OneDimTable),

    /// 多维表：多个维度 -> (Cp, Ct)
    MultiDimensional(MultiDimTable),
}

impl CpCtTable {
    /// 获取Cp值
    pub fn get_cp(&self, conditions: &TableConditions) -> Result<f64, CpCtError> {
        match self {
            CpCtTable::OneDimensional(table) => table.get_cp(conditions.wind_speed),
            CpCtTable::MultiDimensional(table) => table.get_cp(conditions),
        }
    }

    /// 获取Ct值
    pub fn get_ct(&self, conditions: &TableConditions) -> Result<f64, CpCtError> {
        match self {
            CpCtTable::OneDimensional(table) => table.get_ct(conditions.wind_speed),
            CpCtTable::MultiDimensional(table) => table.get_ct(conditions),
        }
    }

    /// 批量获取Cp值
    pub fn batch_get_cp(&self, conditions_list: &[TableConditions]) -> Result<Vec<f64>, CpCtError> {
        conditions_list
            .iter()
            .map(|cond| self.get_cp(cond))
            .collect()
    }

    /// 批量获取Ct值
    pub fn batch_get_ct(&self, conditions_list: &[TableConditions]) -> Result<Vec<f64>, CpCtError> {
        conditions_list
            .iter()
            .map(|cond| self.get_ct(cond))
            .collect()
    }

    /// 获取表类型
    pub fn table_type(&self) -> TableType {
        match self {
            CpCtTable::OneDimensional(_) => TableType::OneDimensional,
            CpCtTable::MultiDimensional(_) => TableType::MultiDimensional,
        }
    }

    /// 获取支持的维度
    pub fn supported_dimensions(&self) -> Vec<Dimension> {
        match self {
            CpCtTable::OneDimensional(_) => vec![Dimension::WindSpeed],
            CpCtTable::MultiDimensional(table) => table
                .dimensions
                .iter()
                .map(|d| d.dimension.clone())
                .collect(),
        }
    }
    pub fn wind_speeds(&self) -> &[f64] {
        match self {
            CpCtTable::OneDimensional(table) => table.wind_speeds.as_slice().unwrap_or(&[]),
            CpCtTable::MultiDimensional(table) => table
                .dimensions
                .iter()
                .find(|d| d.dimension == Dimension::WindSpeed)
                .map(|d| d.values.as_slice())
                .unwrap_or(&[]),
        }
    }
}

// ============================================================================
// 表类型定义
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableType {
    OneDimensional,
    MultiDimensional,
    Hybrid,
}

// ============================================================================
// 一维表实现
// ============================================================================

/// 一维Cp/Ct表
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct OneDimTable {
    pub wind_speeds: Array1,
    pub cp_values: Array1,
    pub ct_values: Array1,
    pub interpolation: InterpolationMethod,
}

impl OneDimTable {
    // pub fn new(config: &OneDimTableConfig) -> Result<Self, CpCtError> {
    //     // 验证数据
    //     if config.wind_speeds.len() != config.cp_values.len()
    //         || config.wind_speeds.len() != config.ct_values.len()
    //     {
    //         return Err(CpCtError::DimensionMismatch(
    //             "风速、Cp、Ct数组长度必须相同".to_string(),
    //         ));
    //     }

    //     // 验证风速单调递增
    //     for i in 1..config.wind_speeds.len() {
    //         if config.wind_speeds[i] <= config.wind_speeds[i - 1] {
    //             return Err(CpCtError::ValidationFailed(
    //                 "风速值必须单调递增".to_string(),
    //             ));
    //         }
    //     }

    //     // 验证Cp值范围
    //     for &cp in &config.cp_values {
    //         if cp < 0.0 || cp > 0.593 {
    //             return Err(CpCtError::ValidationFailed(format!(
    //                 "Cp值超出合理范围: {}",
    //                 cp
    //             )));
    //         }
    //     }

    //     // 验证Ct值范围
    //     for &ct in &config.ct_values {
    //         if ct < 0.0 || ct > 1.0 {
    //             return Err(CpCtError::ValidationFailed(format!(
    //                 "Ct值超出合理范围: {}",
    //                 ct
    //             )));
    //         }
    //     }

    //     Ok(Self {
    //         wind_speeds: config.wind_speeds.clone(),
    //         cp_values: config.cp_values.clone(),
    //         ct_values: config.ct_values.clone(),
    //         interpolation: config.interpolation,
    //     })
    // }

    pub fn get_cp(&self, wind_speed: f64) -> Result<f64, CpCtError> {
        let slice = self
            .cp_values
            .as_slice()
            .ok_or_else(|| CpCtError::ValidationFailed("无法获取cp_values的切片".to_string()))?;
        let cp = self.interpolate(slice, wind_speed)?;
        Ok(cp)
    }

    pub fn get_ct(&self, wind_speed: f64) -> Result<f64, CpCtError> {
        let slice = self
            .ct_values
            .as_slice()
            .ok_or_else(|| CpCtError::ValidationFailed("无法获取ct_values的切片".to_string()))?;
        let ct = self.interpolate(slice, wind_speed)?;
        Ok(ct)
    }

    fn interpolate(&self, values: &[f64], wind_speed: f64) -> Result<f64, CpCtError> {
        // 边界检查
        if wind_speed <= self.wind_speeds[0] {
            return Ok(values[0]);
        }

        if wind_speed >= self.wind_speeds[self.wind_speeds.len() - 1] {
            return Ok(values[values.len() - 1]);
        }

        // 查找插值区间
        for i in 0..self.wind_speeds.len() - 1 {
            if wind_speed >= self.wind_speeds[i] && wind_speed <= self.wind_speeds[i + 1] {
                return match self.interpolation {
                    InterpolationMethod::Linear => {
                        let t = (wind_speed - self.wind_speeds[i])
                            / (self.wind_speeds[i + 1] - self.wind_speeds[i]);
                        Ok(values[i] * (1.0 - t) + values[i + 1] * t)
                    }
                    InterpolationMethod::CubicSpline => {
                        self.cubic_spline_interpolate(i, wind_speed, values)
                    }
                    InterpolationMethod::NearestNeighbor => {
                        let dist1 = (wind_speed - self.wind_speeds[i]).abs();
                        let dist2 = (wind_speed - self.wind_speeds[i + 1]).abs();
                        Ok(if dist1 <= dist2 {
                            values[i]
                        } else {
                            values[i + 1]
                        })
                    }
                    InterpolationMethod::Akima => self.akima_interpolate(i, wind_speed, values),
                };
            }
        }

        Err(CpCtError::InterpolationFailed(format!(
            "找不到风速 {} 的插值区间",
            wind_speed
        )))
    }

    fn cubic_spline_interpolate(
        &self,
        idx: usize,
        x: f64,
        values: &[f64],
    ) -> Result<f64, CpCtError> {
        // 简化实现：实际需要预计算样条系数
        // 这里使用线性插值作为占位符
        let t = (x - self.wind_speeds[idx]) / (self.wind_speeds[idx + 1] - self.wind_speeds[idx]);
        Ok(values[idx] * (1.0 - t) + values[idx + 1] * t)
    }

    fn akima_interpolate(&self, idx: usize, x: f64, values: &[f64]) -> Result<f64, CpCtError> {
        // 简化实现：实际需要Akima插值算法
        let t = (x - self.wind_speeds[idx]) / (self.wind_speeds[idx + 1] - self.wind_speeds[idx]);
        Ok(values[idx] * (1.0 - t) + values[idx + 1] * t)
    }

    // pub fn validate(&self) -> Result<(), CpCtError> {
    //     // 验证数据一致性
    //     if self.wind_speeds.len() < 2 {
    //         return Err(CpCtError::ValidationFailed("至少需要2个风速点".to_string()));
    //     }

    //     Ok(())
    // }

    // pub fn wind_speed_range(&self) -> (f64, f64) {
    //     (
    //         self.wind_speeds[0],
    //         self.wind_speeds[self.wind_speeds.len() - 1],
    //     )
    // }
}

// ============================================================================
// 多维表实现
// ============================================================================

/// 多维Cp/Ct表
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MultiDimTable {
    dimensions: Vec<DimensionRange>,
    cp_data: DynArray,
    ct_data: DynArray,
}

impl MultiDimTable {
    pub fn new(
        dimensions: Vec<DimensionRange>,
        cp_data: DynArray,
        ct_data: DynArray,
    ) -> Result<Self, CpCtError> {
        // 验证维度一致性
        Self::validate_dimensions(&dimensions, &cp_data, &ct_data)?;

        Ok(Self {
            dimensions: dimensions.to_vec(),
            cp_data,
            ct_data,
        })
    }
    fn validate_dimensions(
        dimensions: &[DimensionRange],
        cp_data: &DynArray,
        ct_data: &DynArray,
    ) -> Result<(), CpCtError> {
        Self::validate_inline_data(dimensions, cp_data, ct_data)
    }
    fn validate_inline_data(
        dimensions: &[DimensionRange],
        cp_data: &DynArray,
        ct_data: &DynArray,
    ) -> Result<(), CpCtError> {
        // 验证维度数量
        if dimensions.len() != cp_data.ndim() {
            return Err(CpCtError::DimensionMismatch(format!(
                "维度数量不匹配: 期望{}，实际{}",
                dimensions.len() + 1,
                cp_data.ndim()
            )));
        }

        // 验证数据形状
        if cp_data.shape() != ct_data.shape() {
            return Err(CpCtError::DimensionMismatch(
                "Cp和Ct数据形状不匹配".to_string(),
            ));
        }

        // 验证维度值数量与数据形状匹配
        for (i, dim) in dimensions.iter().enumerate() {
            if dim.values.len() != cp_data.shape()[i] {
                return Err(CpCtError::DimensionMismatch(format!(
                    "维度{}的值数量不匹配",
                    i
                )));
            }
        }

        Ok(())
    }
    // fn load_from_file(
    //     path: &Path,
    //     format: &FileFormat,
    // ) -> Result<(Vec<DimensionRange>, DynArray, DynArray), CpCtError> {
    //     match format {
    //         FileFormat::Csv => csv_loader::CsvDataLoader::load_multidimensional_data(path),
    //         FileFormat::Npy => NpyDataLoader::load_npy_file(path),
    //         _ => Err(CpCtError::FileFormat(format!(
    //             "不支持的文件格式: {:?}",
    //             format
    //         ))),
    //     }
    // }
    pub fn get_cp(&self, conditions: &TableConditions) -> Result<f64, CpCtError> {
        let query_point = self.build_query_point(conditions);
        let cp = self.interpolate_nd(&self.cp_data, &query_point)?;
        Ok(cp)
    }

    pub fn get_ct(&self, conditions: &TableConditions) -> Result<f64, CpCtError> {
        let query_point = self.build_query_point(conditions);
        let ct = self.interpolate_nd(&self.ct_data, &query_point)?;
        Ok(ct)
    }

    fn build_query_point(&self, conditions: &TableConditions) -> Vec<f64> {
        self.dimensions
            .iter()
            .map(|dim_range| match dim_range.dimension {
                Dimension::WindSpeed => conditions.wind_speed,
                Dimension::TurbulenceIntensity => conditions.turbulence_intensity.unwrap_or(0.06),
                Dimension::YawAngle => conditions.yaw_angle.unwrap_or(0.0),
                Dimension::AirDensity => conditions.air_density.unwrap_or(1.225),
                Dimension::BladePitch => conditions.blade_pitch.unwrap_or(0.0),
                Dimension::RotorSpeed => conditions.rotor_speed.unwrap_or(1.0),
                Dimension::WindDirection => conditions.wind_direction.unwrap_or(0.0),
                Dimension::ShearExponent => conditions.shear_exponent.unwrap_or(0.2),
                Dimension::InflowAngle => conditions.inflow_angle.unwrap_or(0.0),
                Dimension::WavePeriod => conditions.wave_period.unwrap_or(0.0),
                Dimension::WaveHeight => conditions.wave_height.unwrap_or(0.0),
                Dimension::Custom(ref name) => conditions
                    .custom_dimensions
                    .get(name)
                    .copied()
                    .unwrap_or(0.0),
            })
            .collect()
    }

    fn interpolate_nd(&self, data: &DynArray, point: &[f64]) -> Result<f64, CpCtError> {
        match self.dimensions.len() {
            0 => Err(CpCtError::DimensionMismatch("没有维度".to_string())),
            1 => self.interpolate_1d(data, point[0]),
            2 => self.interpolate_2d(data, point[0], point[1]),
            3 => self.interpolate_3d(data, point[0], point[1], point[2]),
            _ => panic!("不支持超过3维的插值"),
        }
    }

    fn interpolate_1d(&self, data: &DynArray, x: f64) -> Result<f64, CpCtError> {
        let values = &self.dimensions[0].values;

        if x <= values[0] {
            return Ok(data[[0]]);
        }

        if x >= values[values.len() - 1] {
            return Ok(data[[values.len() - 1]]);
        }

        for i in 0..values.len() - 1 {
            if x >= values[i] && x <= values[i + 1] {
                let t = (x - values[i]) / (values[i + 1] - values[i]);
                return Ok(data[[i]] * (1.0 - t) + data[[i + 1]] * t);
            }
        }

        Err(CpCtError::InterpolationFailed("一维插值失败".to_string()))
    }

    fn interpolate_2d(&self, data: &DynArray, x: f64, y: f64) -> Result<f64, CpCtError> {
        let x_values = &self.dimensions[0].values;
        let y_values = &self.dimensions[1].values;

        // 找到边界索引
        let x_idx = self.find_index(x_values, x)?;
        let y_idx = self.find_index(y_values, y)?;

        // 双线性插值
        let q11 = data[[x_idx, y_idx]];
        let q12 = data[[x_idx, y_idx + 1]];
        let q21 = data[[x_idx + 1, y_idx]];
        let q22 = data[[x_idx + 1, y_idx + 1]];

        let x1 = x_values[x_idx];
        let x2 = x_values[x_idx + 1];
        let y1 = y_values[y_idx];
        let y2 = y_values[y_idx + 1];

        let result = q11 * (x2 - x) * (y2 - y)
            + q21 * (x - x1) * (y2 - y)
            + q12 * (x2 - x) * (y - y1)
            + q22 * (x - x1) * (y - y1);

        Ok(result / ((x2 - x1) * (y2 - y1)))
    }

    fn interpolate_3d(&self, data: &DynArray, x: f64, y: f64, z: f64) -> Result<f64, CpCtError> {
        // 三线性插值实现
        let x_values = &self.dimensions[0].values;
        let y_values = &self.dimensions[1].values;
        let z_values = &self.dimensions[2].values;

        let x_idx = self.find_index(x_values, x)?;
        let y_idx = self.find_index(y_values, y)?;
        let z_idx = self.find_index(z_values, z)?;

        // 获取8个角点的值
        let c000 = data[[x_idx, y_idx, z_idx]];
        let c001 = data[[x_idx, y_idx, z_idx + 1]];
        let c010 = data[[x_idx, y_idx + 1, z_idx]];
        let c011 = data[[x_idx, y_idx + 1, z_idx + 1]];
        let c100 = data[[x_idx + 1, y_idx, z_idx]];
        let c101 = data[[x_idx + 1, y_idx, z_idx + 1]];
        let c110 = data[[x_idx + 1, y_idx + 1, z_idx]];
        let c111 = data[[x_idx + 1, y_idx + 1, z_idx + 1]];

        let x1 = x_values[x_idx];
        let x2 = x_values[x_idx + 1];
        let y1 = y_values[y_idx];
        let y2 = y_values[y_idx + 1];
        let z1 = z_values[z_idx];
        let z2 = z_values[z_idx + 1];

        let xd = (x - x1) / (x2 - x1);
        let yd = (y - y1) / (y2 - y1);
        let zd = (z - z1) / (z2 - z1);

        // 三线性插值公式
        let c00 = c000 * (1.0 - xd) + c100 * xd;
        let c01 = c001 * (1.0 - xd) + c101 * xd;
        let c10 = c010 * (1.0 - xd) + c110 * xd;
        let c11 = c011 * (1.0 - xd) + c111 * xd;

        let c0 = c00 * (1.0 - yd) + c10 * yd;
        let c1 = c01 * (1.0 - yd) + c11 * yd;

        let result = c0 * (1.0 - zd) + c1 * zd;

        Ok(result)
    }

    fn find_index(&self, values: &[f64], target: f64) -> Result<usize, CpCtError> {
        if target <= values[0] {
            return Ok(0);
        }

        if target >= values[values.len() - 1] {
            return Ok(values.len() - 2);
        }

        for i in 0..values.len() - 1 {
            if target >= values[i] && target <= values[i + 1] {
                return Ok(i);
            }
        }

        Err(CpCtError::InterpolationFailed(format!(
            "找不到值 {} 的索引",
            target
        )))
    }

    fn calculate_weight(&self, point: &[f64], indices: &[usize]) -> f64 {
        let mut weight = 1.0;

        for i in 0..self.dimensions.len() {
            let dim = &self.dimensions[i];
            let idx = indices[i];

            if idx >= dim.values.len() {
                return 0.0;
            }

            let value = dim.values[idx];
            let distance = (point[i] - value).abs();

            // 使用距离的倒数作为权重（简化实现）
            // 实际应该使用更复杂的权重函数
            weight *= 1.0 / (1.0 + distance);
        }
        weight
    }
}

/// 维度范围
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq)]
pub struct DimensionRange {
    pub dimension: Dimension,
    pub values: Vec<f64>,
    #[serde(default = "default_interpolation")]
    pub interpolation: InterpolationMethod,
}
fn default_interpolation() -> InterpolationMethod {
    InterpolationMethod::Linear
}

/// 维度类型
#[derive(Debug, Deserialize, Serialize, Clone, PartialEq, Eq, Hash)]
pub enum Dimension {
    #[serde(rename = "wind_speed")]
    WindSpeed,
    #[serde(rename = "turbulence_intensity")]
    TurbulenceIntensity,
    #[serde(rename = "yaw_angle")]
    YawAngle,
    #[serde(rename = "air_density")]
    AirDensity,
    #[serde(rename = "blade_pitch")]
    BladePitch,
    #[serde(rename = "rotor_speed")]
    RotorSpeed,
    #[serde(rename = "wind_direction")]
    WindDirection,
    #[serde(rename = "shear_exponent")]
    ShearExponent,
    #[serde(rename = "inflow_angle")]
    InflowAngle,
    #[serde(rename = "wave_period")]
    WavePeriod,
    #[serde(rename = "wave_height")]
    WaveHeight,
    #[serde(rename = "custom")]
    Custom(String),
}

/// 插值方法
#[derive(Debug, Deserialize, Serialize, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InterpolationMethod {
    #[serde(rename = "linear")]
    Linear,
    #[serde(rename = "cubic_spline")]
    CubicSpline,
    #[serde(rename = "nearest_neighbor")]
    NearestNeighbor,
    #[serde(rename = "akima")]
    Akima,
}

/// 文件格式
#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum FileFormat {
    #[serde(rename = "csv")]
    Csv,
    #[serde(rename = "npy")]
    Npy,
}

/// 表查询条件
#[derive(Debug, Clone, Builder)]
#[builder(setter(into))]
pub struct TableConditions {
    pub wind_speed: f64,
    pub turbulence_intensity: Option<f64>,
    pub yaw_angle: Option<f64>,
    pub air_density: Option<f64>,
    pub blade_pitch: Option<f64>,
    pub rotor_speed: Option<f64>,
    pub wind_direction: Option<f64>,
    pub shear_exponent: Option<f64>,
    pub inflow_angle: Option<f64>,
    pub wave_period: Option<f64>,
    pub wave_height: Option<f64>,
    pub custom_dimensions: HashMap<String, f64>,
}

impl TableConditions {
    pub fn builder() -> TableConditionsBuilder {
        TableConditionsBuilder::default()
    }
}

impl Default for TableConditions {
    fn default() -> Self {
        Self {
            wind_speed: 0.0,
            turbulence_intensity: Some(0.06),
            yaw_angle: Some(0.0),
            air_density: Some(1.225),
            blade_pitch: Some(0.0),
            rotor_speed: Some(1.0),
            wind_direction: Some(0.0),
            shear_exponent: Some(0.2),
            inflow_angle: Some(0.0),
            wave_period: Some(0.0),
            wave_height: Some(0.0),
            custom_dimensions: HashMap::new(),
        }
    }
}

struct Neighbor;
