# FLORIS-RS

Rust implementation of FLORIS - A controls-oriented wind farm wake modeling software.

This is a translation of the Python FLORIS project (v4.6) to Rust, providing improved performance and memory safety for wind farm simulations.

## Overview

FLORIS (FLOw Redirection and Induction in Steady State) is a wake modeling and wind farm controls software incorporating steady-state engineering wake models into a performance-focused framework. This Rust implementation aims to maintain compatibility with the original Python version while leveraging Rust's performance benefits.

## Project Status

This is a **work-in-progress** translation. Currently implemented:

### ✅ Completed
- Core data structures (`Farm`, `FlowField`, `Turbine`, `Grid`)
- Type system with ndarray support
- Utility functions (coordinate transformations, YAML loading)
- Wake models (Gaussian wake velocity, Jimenez deflection)
- Wind data structures (`TimeSeries`, `WindRose`)
- FlorisModel main interface
- Basic turbine operation models

### 🚧 In Progress
- Complete wake solver implementations
- Wake superposition logic
- Power and thrust calculations
- Grid initialization

### 📋 TODO
- Optimization modules
- Heterogeneous inflow handling
- Floating turbine support
- Parallel computation support
- Visualization tools
- Full test coverage
- Documentation
- Turbine library YAML files

## Structure

```
src/
├── lib.rs              # Library root with re-exports
├── main.rs             # Example binary
├── types.rs            # Type definitions and conversions
├── utilities.rs        # Utility functions
├── floris_model.rs     # Main FlorisModel interface
├── wind_data.rs        # Wind data structures
├── core/               # Core simulation components
│   ├── flow_field.rs   # Flow field representation
│   ├── farm.rs         # Wind farm layout and turbines
│   ├── turbine.rs      # Turbine model
│   ├── grid.rs         # Computational grids
│   ├── state.rs        # Solver state
│   └── base.rs         # Base traits
├── wake/               # Wake models
│   ├── wake_velocity.rs     # Wake deficit models
│   ├── wake_deflection.rs   # Wake steering models
│   ├── wake_turbulence.rs   # Turbulence models
│   └── wake_combination.rs  # Wake superposition
└── turbine/            # Turbine operations
    └── operation_models.rs  # Control models
```

## Building

```bash
# Build the project
cargo build --release

# Run tests
cargo test

# Run the example
cargo run --release

# Generate documentation
cargo doc --open
```

## Usage Example

```rust
use florus::{FlorisModel, Array1};

fn main() -> anyhow::Result<()> {
    // Load configuration from file
    let mut model = FlorisModel::from_file("config.yaml")?;
    
    // Set wind conditions
    let wind_speeds = Array1::from_vec(vec![8.0, 10.0, 12.0]);
    let wind_directions = Array1::from_vec(vec![270.0, 280.0, 290.0]);
    let turbulence_intensities = Array1::from_vec(vec![0.06, 0.08, 0.07]);
    
    model.set_wind_conditions(
        wind_speeds,
        wind_directions,
        turbulence_intensities,
    )?;
    
    // Run simulation
    model.run()?;
    
    // Get results
    let powers = model.get_turbine_powers();
    let farm_power = model.get_farm_power();
    
    println!("Farm power: {:?}", farm_power);
    
    Ok(())
}
```

## Dependencies

Key dependencies:
- `ndarray` - N-dimensional arrays (NumPy equivalent)
- `serde` / `serde_yaml` - Serialization and configuration
- `anyhow` / `thiserror` - Error handling
- `rayon` - Parallel computation
- `nalgebra` - Linear algebra

## Original Project

This is based on the Python FLORIS project:
- Repository: https://github.com/NREL/floris
- Documentation: https://nrel.github.io/floris
- Version: 4.6

## License

BSD-3-Clause (matching original FLORIS license)

## Contributing

This is an ongoing translation effort. Contributions are welcome! Key areas needing work:
1. Complete solver implementations
2. Port remaining wake models
3. Add comprehensive tests
4. Improve documentation
5. Performance optimization

## Differences from Python Version

- Uses Rust's type system for safety
- Leverages `ndarray` instead of NumPy
- Trait-based design for extensibility
- Compile-time guarantees for correctness
- Potential for better performance through zero-cost abstractions

## Contact

For questions about this Rust implementation, please open an issue on the repository.

For the original Python FLORIS, see the [NREL FLORIS repository](https://github.com/NREL/floris).
