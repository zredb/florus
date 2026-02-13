pub mod aep;
pub mod core;
pub mod floris_config;
pub mod floris_model;
pub mod heterogeneous_map;
pub mod optimization;
/// FLORIS-RS: Rust implementation of FLORIS wind farm wake modeling software
///
/// This is a translation of Python FLORIS project (v4.6) to Rust for
/// improved performance and safety.
pub mod types;
pub mod utilities;

pub mod wind_data;

// Re-export commonly used types
pub use floris_config::FlorisConfig;
pub use floris_model::FlorisModel;
pub use types::{Array1, Array2, Array3, Array4, Float};

/// Result type for FLORIS operations
pub type Result<T> = std::result::Result<T, anyhow::Error>;

pub enum OneOrManyD1 {
    One(Float),
    Many(Array1),
}
pub enum OneOrManyD2 {
    One(Float),
    Many(Array2),
}
pub enum OneOrManyD3 {
    One(Float),
    Many(Array3),
}
pub enum OneOrManyD4 {
    One(Float),
    Many(Array4),
}
