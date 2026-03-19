# FLORIS-RS Agent Guidelines

## Build/Test Commands
- **Build**: `cargo build --release`
- **Single test**: `cargo test <test_name>` (e.g., `cargo test solver::tests::test_axial_induction_low_ct`)
- **Run all tests**: `cargo test`
- **Check formatting**: `cargo fmt -- --check`
- **Lint**: `cargo clippy -- -D warnings`

## Code Style

### Imports & Organization
- Use absolute crate paths: `crate::core::wake::WakeModelManager`
- Group imports: std → crates → local (`use std::...; use crate::...;`)
- Doc comments: `///` for public API, `//!` for module-level

### Naming Conventions
- Types: `PascalCase` (structs, enums, traits)
- Functions/variables: `snake_case`
- Constants: `SCREAMING_SNAKE_CASE`
- Array types: `Array1`, `Array2`, `Array3`, `Array4` (see `types.rs`)
- Float type alias: `Float = f64`

### Error Handling
- Return type: `crate::Result<T>` = `anyhow::Result<T>`
- Use `anyhow::bail!("...")` for early returns
- Avoid `unwrap()` in production code (use `?` operator)
- Never suppress errors with `as any`

### Array Indexing
- 4D arrays: `[findex, turbine, y, z]`
- Use `ndarray::s![]` macro for slicing
- Clone arrays before mutation to avoid borrow checker issues

### Testing
- Tests go in `#[cfg(test)]` module at end of file
- Use `approx::assert_relative_eq!` for floating-point comparisons
- Mock complex dependencies, don't mock primitives

### Wake Models
- Traits: `VelocityModel`, `DeflectionModel`, `TurbulenceModel`, `CombinationModel`
- Implement in `src/core/wake/{submodule}/`
- Use `Box<dyn Trait>` for dynamic dispatch
