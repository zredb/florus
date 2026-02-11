//! Build script to convert turbopark_lookup_table.mat to embedded Rust data.
//!
//! This script is called during cargo build to read the MATLAB lookup table
//! and generate a Rust source file with the data as const arrays.

use std::env;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let mat_file = "src/core/wake/wake_velocity/turbopark_lookup_table.mat";

    // Check if the .mat file exists
    if !Path::new(mat_file).exists() {
        println!("cargo:warning=turbopark_lookup_table.mat not found, using runtime calculation");
        return;
    }

    // Get the output directory
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let output_file = out_dir.join("turbopark_lookup.rs");

    // Get the modification times
    let mat_mtime = fs::metadata(mat_file)
        .and_then(|m| m.modified())
        .ok();

    let needs_rebuild = if let Ok(output_metadata) = fs::metadata(&output_file) {
        if let (Some(mat_m), Ok(output_m)) = (mat_mtime, output_metadata.modified()) {
            mat_m > output_m
        } else {
            true
        }
    } else {
        true
    };

    if needs_rebuild {
        println!("cargo:rerun-if-changed={}", mat_file);

        // Try to run Python with scipy
        let python_exes = vec![
            "python3",
            "python",
        ];

        let mut python_found = None;
        for &python in &python_exes {
            let output = Command::new(python)
                .args(&["--version"])
                .output();

            if let Ok(o) = output {
                if o.status.success() {
                    python_found = Some(python);
                    break;
                }
            }
        }

        if let Some(python) = python_found {
            println!("cargo:warning=Using {} to convert .mat file", python);

            // Run the conversion script with output directory as argument
            let status = Command::new(python)
                .arg("build_rs_contents.py")
                .arg("--out-dir")
                .arg(&out_dir)
                .current_dir(env::current_dir().unwrap())
                .status();

            match status {
                Ok(s) if s.success() => {
                    // The script writes to OUT_DIR/turbopark_lookup.rs
                    if output_file.exists() {
                        println!("cargo:warning=Successfully generated turbopark_lookup.rs");
                    } else {
                        println!("cargo:warning=Generated file not found, using runtime calculation");
                    }
                }
                Ok(_) => {
                    println!("cargo:warning=Python script failed, using runtime calculation");
                }
                Err(e) => {
                    println!("cargo:warning=Could not run Python: {}, using runtime calculation", e);
                }
            }
        } else {
            println!("cargo:warning=Python not found, using runtime calculation");
        }
    }
}
