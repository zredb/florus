/// Type definitions for FLORIS-RS
///
/// Corresponds to type_dec.py in Python implementation

use ndarray::{Array1 as NdArray1, Array2 as NdArray2, Array3 as NdArray3, Array4 as NdArray4};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Primary floating point type used throughout FLORIS
/// Corresponds to floris_float_type in Python (np.float64)
pub type Float = f64;

/// 1D array of floats
pub type Array1 = NdArray1<Float>;

/// 2D array of floats
pub type Array2 = NdArray2<Float>;

/// 3D array of floats
pub type Array3 = NdArray3<Float>;

/// 4D array of floats
pub type Array4 = NdArray4<Float>;

/// Numeric dictionary type for configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NumericDict {
    #[serde(flatten)]
    pub data: HashMap<String, ConfigValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ConfigValue {
    Float(Float),
    FloatArray(Vec<Float>),
    String(String),
    Bool(bool),
}

impl NumericDict {
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    pub fn get_scalar(&self, key: &str) -> Option<Float> {
        match self.data.get(key)? {
            ConfigValue::Float(v) => Some(*v),
            _ => None,
        }
    }

    pub fn get_array(&self, key: &str) -> Option<&[Float]> {
        match self.data.get(key)? {
            ConfigValue::FloatArray(v) => Some(v.as_slice()),
            _ => None,
        }
    }

    pub fn get_string(&self, key: &str) -> Option<&str> {
        match self.data.get(key)? {
            ConfigValue::String(v) => Some(v.as_str()),
            _ => None,
        }
    }

    pub fn get_bool(&self, key: &str) -> Option<bool> {
        match self.data.get(key)? {
            ConfigValue::Bool(v) => Some(*v),
            _ => None,
        }
    }
}

impl Default for NumericDict {
    fn default() -> Self {
        Self {
            data: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_numeric_dict() {
        let mut dict = NumericDict::new();
        dict.data.insert("key1".to_string(), ConfigValue::Float(1.0));
        assert_eq!(dict.get_scalar("key1"), Some(1.0));
    }
}
