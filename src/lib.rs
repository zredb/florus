/// FLORIS-RS: Rust implementation of FLORIS wind farm wake modeling software
/// 
/// This is a translation of the Python FLORIS project (v4.6) to Rust for 
/// improved performance and safety.

pub mod types;
pub mod core;
pub mod floris_model;
pub mod utilities;
pub mod wind_data;
pub mod turbine;
pub mod aep;
pub mod optimization;

// Re-export commonly used types
pub use floris_model::FlorisModel;
pub use types::{Float, Array1, Array2, Array3, Array4};

/// Result type for FLORIS operations
pub type Result<T> = std::result::Result<T, anyhow::Error>;
