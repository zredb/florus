//! Build script for FLORUS.
//!
//! It copies `default_inputs.yaml` into `OUT_DIR` so that it is available
//! next to the other build artifacts.
//!
//! No build-time code generation takes place: the TurbOPark lookup table is
//! stored as a pre-generated binary (`data/turbopark_overlap.bin`) and is
//! embedded into the crate with `include_bytes!` in
//! `src/core/wake/wake_velocity/turbopark_lookup_matlab.rs`. Building the
//! crate therefore needs neither Python nor SciPy.

use std::env;
use std::fs;
use std::path::PathBuf;

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let default_input = "default_inputs.yaml";
    let output_default_input_file = out_dir.join(default_input);
    if !output_default_input_file.exists() {
        fs::copy(default_input, &output_default_input_file)
            .expect("Failed to copy default_inputs.yaml");
    }

    println!("cargo::rerun-if-changed=default_inputs.yaml");
}
