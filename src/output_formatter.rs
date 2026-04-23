/// Output formatting utilities to match Python FLORIS output style
/// 
/// This module provides helper functions to format arrays and values
/// in a way that matches Python's numpy output.

use ndarray::{Array1, Array2, Array3};

/// Format a 1D array similar to Python's numpy print
pub fn format_array_1d(arr: &Array1<f64>) -> String {
    let mut result = String::from("[");
    for (i, val) in arr.iter().enumerate() {
        if i > 0 {
            result.push(' ');
        }
        result.push_str(&format!("{:.8}", val));
    }
    result.push(']');
    result
}

/// Format a 2D array similar to Python's numpy print
pub fn format_array_2d(arr: &Array2<f64>) -> String {
    let rows = arr.shape()[0];
    let cols = arr.shape()[1];
    let mut result = String::new();
    
    for i in 0..rows {
        if i == 0 {
            result.push('[');
        } else {
            result.push(' ');
        }
        result.push('[');
        for j in 0..cols {
            if j > 0 {
                result.push(' ');
            }
            result.push_str(&format!("{:.8}", arr[[i, j]]));
        }
        result.push(']');
        if i < rows - 1 {
            result.push('\n');
        } else {
            result.push(']');
        }
    }
    result
}

/// Format a 3D array shape description
pub fn format_shape_3d(arr: &Array3<f64>) -> String {
    format!("({}, {}, {})", arr.shape()[0], arr.shape()[1], arr.shape()[2])
}

/// Format a 2D array shape description
pub fn format_shape_2d<T>(arr: &ndarray::ArrayBase<T, ndarray::Dim<[usize; 2]>>) -> String 
where
    T: ndarray::Data,
{
    format!("({}, {})", arr.shape()[0], arr.shape()[1])
}

/// Format a 1D array shape description
pub fn format_shape_1d<T>(arr: &ndarray::ArrayBase<T, ndarray::Dim<[usize; 1]>>) -> String 
where
    T: ndarray::Data,
{
    format!("({},)", arr.shape()[0])
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr1;

    #[test]
    fn test_format_array_1d() {
        let arr = arr1(&[1.23456789, 2.34567890, 3.45678901]);
        let formatted = format_array_1d(&arr);
        assert!(formatted.starts_with('['));
        assert!(formatted.ends_with(']'));
    }

    #[test]
    fn test_format_shape() {
        let arr = Array2::<f64>::zeros((2, 3));
        let shape = format_shape_2d(&arr);
        assert_eq!(shape, "(2, 3)");
    }
}
