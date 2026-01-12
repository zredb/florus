//! Wind data structures for FLORIS-RS
//!
//! Provides wind data objects to hold ambient wind conditions including:
//! - TimeSeries: Time series wind data
//! - WindRose: Aggregated wind statistics by direction/speed bins
//! - WindTIRose: Wind rose with TI as an additional dimension
//! - WindRoseWRG: Wind rose from WAsP WRG file
//! - WindRoseByTurbine: Wind rose with separate wind rose for each turbine
//!
//! Corresponds to wind_data.py in Python FLORIS v4.6

pub mod traits;
pub mod timeseries;
pub mod wind_rose;
pub mod wind_ti_rose;
pub mod wind_rose_wrg;

pub use traits::{TIParams, WindData};
pub use timeseries::TimeSeries;
pub use wind_rose::WindRose;
pub use wind_ti_rose::WindTIRose;
pub use wind_rose_wrg::{RegularGridInterpolant, WindRoseByTurbine, WindRoseWRG, WRGData};
