/// Wake turbulence models
///
/// Corresponds to wake_turbulence/ module in Python implementation

pub mod none;
pub mod crespo_hernandez;

pub use none::NoneTurbulence;
pub use crespo_hernandez::CrespoHernandez;
