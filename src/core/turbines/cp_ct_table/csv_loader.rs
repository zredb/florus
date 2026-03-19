// CSV数据加载器

use super::{CpCtError, Dimension, DimensionRange, InterpolationMethod, MultiDimTable};
use crate::types::Array2;
use crate::types::DynArray;
use csv::Reader;
use ordered_float::OrderedFloat;
use std::collections::{BTreeSet, HashSet};
use std::path::Path;

pub struct CsvDataLoader;

impl CsvDataLoader {
    pub fn load_multidimensional_data<P: AsRef<Path>>(
        file_path: P,
    ) -> Result<MultiDimTable, CpCtError> {
        let mut rdr = Reader::from_path(file_path)?;

        // 读取表头
        let headers = rdr.headers()?.clone();

        // 确定维度
        let mut dimensions = Self::infer_dimensions(&headers)?;

        // 读取数据并提取维度值
        let mut records = Vec::new();
        let mut dimension_values = vec![BTreeSet::new(); dimensions.len()];

        for result in rdr.records() {
            let record = result?;
            records.push(record.clone());

            // 提取维度值
            for (i, _) in dimensions.iter().enumerate() {
                if let Some(value_str) = record.get(i) {
                    if let Ok(value) = value_str.parse::<f64>() {
                        dimension_values[i].insert(OrderedFloat(value));
                    }
                }
            }
        }

        // 更新维度值
        for (i, dim_range) in dimensions.iter_mut().enumerate() {
            let values: Vec<f64> = dimension_values[i].iter().map(|v| v.0).collect();
            dim_range.values = values;
        }

        // 转换为多维数组
        let (cp_data, ct_data) = Self::convert_to_arrays(&dimensions, &records)?;

        Ok(MultiDimTable::new(dimensions, cp_data, ct_data)?)
    }
    fn infer_dimensions(headers: &csv::StringRecord) -> Result<Vec<DimensionRange>, CpCtError> {
        let mut dimensions = Vec::new();

        // 检查列数
        let len = headers.len();
        if len < 5 {
            return Err(CpCtError::InvalidCsvFormat("至少需要5列".to_string()));
        }
        let header_set: Vec<&str> = headers.iter().map(|h| h.trim()).collect();

        // 验证最后三列
        let last_three = &header_set[len - 3..];
        if last_three[0] != "ws"
            || last_three[1] != "power"
            || last_three[2] != "thrust_coefficient"
        {
            return Err(CpCtError::InvalidCsvFormat(
                "最后三列必须是ws, power, thrust_coefficient".to_string(),
            ));
        }

        // 前面的列是维度（根据你的CSV文件，是Tp和Hs, ws）
        for i in 0..len - 2 {
            let dim_name = &headers[i];
            let dimension = match dim_name.to_lowercase().as_str() {
                "tp" => Dimension::WavePeriod,
                "hs" => Dimension::WaveHeight,
                "ti" => Dimension::TurbulenceIntensity,
                "yaw" => Dimension::YawAngle,
                "ws" => Dimension::WindSpeed,
                _ => Dimension::Custom(dim_name.to_string()),
            };

            dimensions.push(DimensionRange {
                dimension: dimension.clone(),
                values: Vec::new(), // 稍后从数据填充
                interpolation: match dimension {
                    Dimension::WindSpeed => InterpolationMethod::Linear,
                    _ => InterpolationMethod::NearestNeighbor,
                },
            });
        }

        Ok(dimensions)
    }

    fn convert_to_arrays(
        dimensions: &[DimensionRange],
        records: &[csv::StringRecord],
    ) -> Result<(DynArray, DynArray), CpCtError> {
        if records.is_empty() {
            return Err(CpCtError::InvalidCsvFormat("没有数据记录".to_string()));
        }

        // 计算每个维度的长度
        let dim_lengths: Vec<usize> = dimensions.iter().map(|d| d.values.len()).collect();

        // 计算总数据点数量
        let total_points: usize = dim_lengths.iter().product();

        // 创建CP和CT数据数组
        let mut cp_data_vec = vec![0.0; total_points];
        let mut ct_data_vec = vec![0.0; total_points];

        // 创建索引映射：从多维索引到线性索引
        let mut index_map = std::collections::HashMap::new();

        for (record_idx, record) in records.iter().enumerate() {
            // 提取维度值
            let mut dim_indices = Vec::new();
            for (i, dim_range) in dimensions.iter().enumerate() {
                let value_str = record.get(i).ok_or_else(|| {
                    CpCtError::InvalidCsvFormat(format!("记录{}缺少第{}列", record_idx, i))
                })?;

                let value = value_str.parse::<f64>().map_err(|e| {
                    CpCtError::InvalidCsvFormat(format!("无法解析值'{}': {}", value_str, e))
                })?;

                // 找到该值在维度值列表中的索引
                let dim_index = dim_range
                    .values
                    .iter()
                    .position(|&v| (v - value).abs() < 1e-10)
                    .ok_or_else(|| {
                        CpCtError::InvalidCsvFormat(format!(
                            "值{}不在维度{}的有效值列表中",
                            value, i
                        ))
                    })?;

                dim_indices.push(dim_index);
            }

            // 提取功率和推力系数
            let power_idx = dimensions.len();
            let ct_idx = dimensions.len() + 1;

            let power = record
                .get(power_idx)
                .ok_or_else(|| {
                    CpCtError::InvalidCsvFormat(format!("记录{}缺少power列", record_idx))
                })?
                .parse::<f64>()
                .map_err(|e| CpCtError::InvalidCsvFormat(format!("无法解析power值: {}", e)))?;

            let thrust_coefficient = record
                .get(ct_idx)
                .ok_or_else(|| {
                    CpCtError::InvalidCsvFormat(format!(
                        "记录{}缺少thrust_coefficient列",
                        record_idx
                    ))
                })?
                .parse::<f64>()
                .map_err(|e| {
                    CpCtError::InvalidCsvFormat(format!("无法解析thrust_coefficient值: {}", e))
                })?;

            // 计算线性索引
            let linear_index = Self::calculate_linear_index(&dim_indices, &dim_lengths);

            cp_data_vec[linear_index] = power;
            ct_data_vec[linear_index] = thrust_coefficient;
            index_map.insert(dim_indices, linear_index);
        }

        // 转换为ndarray数组
        let cp_data = DynArray::from_shape_vec(dim_lengths.clone(), cp_data_vec)
            .map_err(|e| CpCtError::InvalidCsvFormat(format!("无法创建CP数组: {}", e)))?;

        let ct_data = DynArray::from_shape_vec(dim_lengths, ct_data_vec)
            .map_err(|e| CpCtError::InvalidCsvFormat(format!("无法创建CT数组: {}", e)))?;

        Ok((cp_data, ct_data))
    }

    fn calculate_linear_index(indices: &[usize], lengths: &[usize]) -> usize {
        let mut linear_index = 0;
        let mut stride = 1;

        for (i, &idx) in indices.iter().enumerate().rev() {
            linear_index += idx * stride;
            stride *= lengths[i];
        }

        linear_index
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_multidimensional_data() {
        let file = "src/turbine_library/iea_15MW_multi_dim_Tp_Hs.csv";

        // 加载数据
        let table = CsvDataLoader::load_multidimensional_data(file).unwrap();

        // 验证维度 - 根据CSV文件，应该是Tp和Hs两个维度
        assert_eq!(table.dimensions.len(), 3); // Tp, Hs,ws

        // 验证维度值
        let tp_values = &table.dimensions[0].values;
        let hs_values = &table.dimensions[1].values;
        let ws_values = &table.dimensions[2].values;

        // Tp应该有2个唯一值：2和4
        assert_eq!(tp_values.len(), 2);
        assert!(tp_values.contains(&2.0));
        assert!(tp_values.contains(&4.0));

        // Hs应该有2个唯一值：1和5
        assert_eq!(hs_values.len(), 2);
        assert!(hs_values.contains(&1.0));
        assert!(hs_values.contains(&5.0));

        // ws应该有54个唯一值：3和6
        assert_eq!(ws_values.len(), 54);
        assert!(ws_values.contains(&0.0));
        assert!(ws_values.contains(&2.9));
        assert!(ws_values.contains(&3.0));

        // 验证数组形状 - 应该是2x2x2的数组
        assert_eq!(table.cp_data.shape(), &[2, 2, 54]);
        assert_eq!(table.ct_data.shape(), &[2, 2, 54]);

        // 验证数据不为零
        assert!(table.cp_data.sum() > 0.0);
        assert!(table.ct_data.sum() > 0.0);
    }

    #[test]
    fn test_calculate_linear_index() {
        let lengths = vec![2, 2]; // 2x2数组
        let indices = vec![0, 0]; // 第一个元素
        assert_eq!(CsvDataLoader::calculate_linear_index(&indices, &lengths), 0);

        let indices = vec![0, 1]; // 第二列第一个元素
        assert_eq!(CsvDataLoader::calculate_linear_index(&indices, &lengths), 1);

        let indices = vec![1, 0]; // 第二行第一个元素
        assert_eq!(CsvDataLoader::calculate_linear_index(&indices, &lengths), 2);

        let indices = vec![1, 1]; // 最后一个元素
        assert_eq!(CsvDataLoader::calculate_linear_index(&indices, &lengths), 3);
    }
}
