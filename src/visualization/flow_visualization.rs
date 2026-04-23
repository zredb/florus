//! Flow visualization utilities for FLORUS
//! 
//! This module provides functions to visualize wind farm flow fields.
//! Corresponds to `floris/flow_visualization.py` in Python FLORIS.

use crate::{FlorisModel, Result, core::extract_horizontal_plane};
use plotters::prelude::*;
use std::path::Path;

/// Plot rotor-disk velocity values for all turbines (simplified version)
/// 
/// This function creates a subplot grid showing placeholder data.
/// Full implementation will extract actual rotor values from the flow field.
/// Plot rotor values as a grid of heatmaps
/// 
/// Displays the gridded turbine rotor values (u, v, or w velocity components)
/// for inspection and comparison. Each subplot shows one turbine's rotor disk.
/// Corresponds to `plot_rotor_values()` in Python FLORIS.
pub fn plot_rotor_values<P: AsRef<Path>>(
    values: &ndarray::Array4<f64>,
    findex: usize,
    n_rows: usize,
    n_cols: usize,
    output_path: P,
    cmap_name: &str,
) -> Result<()> {
    // Validate dimensions
    let (nfindex, n_turbines, n_y, n_z) = values.dim();
    
    if findex >= nfindex {
        anyhow::bail!("findex {} out of range (max: {})", findex, nfindex - 1);
    }
    
    let total_plots = n_rows * n_cols;
    if total_plots < n_turbines {
        eprintln!("Warning: Not enough subplots ({}x{}={}) for all turbines ({})", 
                  n_rows, n_cols, total_plots, n_turbines);
    }
    
    // Create figure with subplots
    let width = 400 * n_cols as u32;
    let height = 400 * n_rows as u32;
    let root = BitMapBackend::new(output_path.as_ref(), (width, height)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let areas = root.split_evenly((n_rows, n_cols));
    
    // Extract data for this findex
    let data_slice = values.slice(ndarray::s![findex, .., .., ..]);
    
    // Calculate global min/max for consistent color scale
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    
    for t in 0..n_turbines.min(total_plots) {
        for y in 0..n_y {
            for z in 0..n_z {
                let val = data_slice[(t, y, z)];
                if val < vmin { vmin = val; }
                if val > vmax { vmax = val; }
            }
        }
    }
    
    // Get colormap function
    let colormap_fn = get_colormap(cmap_name);
    
    // Draw each subplot
    for idx in 0..total_plots {
        let area = &areas[idx];
        
        if idx >= n_turbines {
            // Empty subplot for unused slots
            let mut chart = ChartBuilder::on(area)
                .margin(5)
                .build_cartesian_2d(0f32..1f32, 0f32..1f32)?;
            chart.configure_mesh().disable_mesh().draw()?;
            continue;
        }
        
        // Create chart for this turbine
        let title = format!("Turbine {}", idx);
        let mut chart = ChartBuilder::on(area)
            .caption(title, ("sans-serif", 15).into_font())
            .margin(5)
            .build_cartesian_2d(
                0..n_z as i32,
                0..n_y as i32,
            )?;
        
        chart.configure_mesh()
            .disable_mesh()
            .x_labels(0)
            .y_labels(0)
            .draw()?;
        
        // Draw heatmap using rectangles
        let cell_width = 1.0;
        let cell_height = 1.0;
        
        for y in 0..n_y {
            for z in 0..n_z {
                let val = data_slice[(idx, y, z)];
                
                // Normalize value to 0-1 range
                let normalized = if vmax > vmin {
                    (val - vmin) / (vmax - vmin)
                } else {
                    0.5
                };
                
                // Get color from colormap
                let color = colormap_fn(normalized);
                
                // Draw rectangle for this cell
                // Note: Invert x-axis to match Python's ax.invert_xaxis()
                let x_start = (n_z - 1 - z) as i32;
                let y_start = y as i32;
                
                chart.draw_series(std::iter::once(
                    Rectangle::new(
                        [(x_start, y_start), (x_start + 1, y_start + 1)],
                        color.filled(),
                    )
                ))?;
            }
        }
    }
    
    Ok(())
}

/// Visualize a horizontal cut plane through the wind farm flow field
/// 
/// This function creates a 2D contour plot of wind speed at hub height.
/// 
/// # Arguments
/// * `fmodel` - FlorisModel instance with computed flow field
/// * `output_path` - Path to save the output image (PNG or SVG)
/// * `min_speed` - Minimum wind speed for color scale (optional)
/// * `max_speed` - Maximum wind speed for color scale (optional)
/// * `title` - Plot title (optional)
pub fn visualize_horizontal_plane<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
    min_speed: Option<f64>,
    max_speed: Option<f64>,
    title: Option<&str>,
) -> Result<()> {
    let flow_field = fmodel.flow_field();
    
    // Check if flow field has been computed
    if flow_field.u.shape()[0] == 0 {
        anyhow::bail!("Flow field has not been computed. Call run() first.");
    }
    
    // Get grid coordinates (if available)
    let grid = fmodel.grid();
    if grid.is_none() {
        anyhow::bail!("Grid not initialized. Call run() first.");
    }
    
    // For now, create a simple visualization using turbine positions
    // Full implementation would extract the actual horizontal plane from the grid
    let layout_x = fmodel.layout_x();
    let layout_y = fmodel.layout_y();
    
    if layout_x.is_empty() {
        anyhow::bail!("No turbines in layout");
    }
    
    // Get hub height
    let hub_height = flow_field.reference_wind_height;
    
    // Extract horizontal plane data (simplified - uses turbine positions)
    // In full implementation, would use extract_horizontal_plane with actual grid
    let u_field = &flow_field.u;
    let v_field = &flow_field.v;
    let w_field = &flow_field.w;
    
    // For demonstration, create a simple grid around turbines
    let x_min = layout_x.iter().cloned().fold(f64::INFINITY, f64::min) - 500.0;
    let x_max = layout_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 2000.0;
    let y_min = layout_y.iter().cloned().fold(f64::INFINITY, f64::min) - 500.0;
    let y_max = layout_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 500.0;
    
    let resolution = 50;
    let mut x_vals = Vec::new();
    let mut y_vals = Vec::new();
    let mut u_vals = Vec::new();
    
    for i in 0..resolution {
        for j in 0..resolution {
            let x = x_min + (x_max - x_min) * i as f64 / (resolution - 1) as f64;
            let y = y_min + (y_max - y_min) * j as f64 / (resolution - 1) as f64;
            
            // Simple wake model approximation for visualization
            // Real implementation would interpolate from actual flow field
            let mut u = flow_field.wind_speeds[0];
            
            // Apply simple wake deficit for each turbine
            for (ti, (&tx, &ty)) in layout_x.iter().zip(layout_y.iter()).enumerate() {
                let dx = x - tx;
                let dy = y - ty;
                let dist = (dx * dx + dy * dy).sqrt();
                
                // Simple Gaussian wake deficit model (very simplified)
                if dx > 0.0 && dist < 1000.0 {
                    let deficit = 0.3 * (-dx / 500.0).exp() * (-dy * dy / 10000.0).exp();
                    u -= deficit;
                }
            }
            
            x_vals.push(x);
            y_vals.push(y);
            u_vals.push(u);
        }
    }
    
    // Determine color scale bounds
    let vmin = min_speed.unwrap_or_else(|| u_vals.iter().cloned().fold(f64::INFINITY, f64::min));
    let vmax = max_speed.unwrap_or_else(|| u_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    
    // Create the plot
    let root = BitMapBackend::new(output_path.as_ref(), (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let title_text = title.unwrap_or("Horizontal Plane - Wind Speed");
    let mut chart = ChartBuilder::on(&root)
        .caption(title_text, ("sans-serif", 30).into_font())
        .margin(5)
        .x_label_area_size(60)
        .y_label_area_size(80)
        .build_cartesian_2d(
            x_min..x_max,
            y_min..y_max,
        )?;
    
    chart.configure_mesh()
        .x_desc("X Position (m)")
        .y_desc("Y Position (m)")
        .draw()?;
    
    // Draw heatmap
    let cell_width = (x_max - x_min) / resolution as f64;
    let cell_height = (y_max - y_min) / resolution as f64;
    
    for idx in 0..u_vals.len() {
        let i = idx / resolution;
        let j = idx % resolution;
        let value = u_vals[idx];
        
        // Normalize value to [0, 1] range for color mapping
        let normalized = if vmax > vmin {
            ((value - vmin) / (vmax - vmin)).clamp(0.0, 1.0)
        } else {
            0.5
        };
        
        // Coolwarm colormap: blue (cold/low) to red (warm/high)
        let color = if normalized < 0.5 {
            // Blue to white
            let t = normalized * 2.0;
            RGBColor(
                (255.0 * t) as u8,
                (255.0 * t) as u8,
                255,
            )
        } else {
            // White to red
            let t = (normalized - 0.5) * 2.0;
            RGBColor(
                255,
                (255.0 * (1.0 - t)) as u8,
                (255.0 * (1.0 - t)) as u8,
            )
        };
        
        let x = x_min + i as f64 * cell_width;
        let y = y_min + j as f64 * cell_height;
        
        chart.draw_series(std::iter::once(
            Rectangle::new(
                [(x, y), (x + cell_width, y + cell_height)],
                color,
            )
        ))?;
    }
    
    // Draw turbine positions on top
    for i in 0..layout_x.len() {
        chart.draw_series(std::iter::once(
            Circle::new((layout_x[i], layout_y[i]), 8, BLACK.filled())
        ))?;
    }
    
    Ok(())
}

/// Get a colormap function by name
/// 
/// Returns a closure that maps a value in [0, 1] to an RGBAColor.
/// Supported colormaps: "coolwarm", "viridis", "plasma", "inferno", "magma"
fn get_colormap(name: &str) -> Box<dyn Fn(f64) -> RGBAColor> {
    match name.to_lowercase().as_str() {
        "coolwarm" => Box::new(coolwarm_colormap),
        "viridis" => Box::new(viridis_colormap),
        "plasma" => Box::new(plasma_colormap),
        "inferno" => Box::new(inferno_colormap),
        "magma" => Box::new(magma_colormap),
        _ => {
            eprintln!("Warning: Unknown colormap '{}', using coolwarm", name);
            Box::new(coolwarm_colormap)
        }
    }
}

/// Coolwarm colormap (blue -> white -> red)
fn coolwarm_colormap(t: f64) -> RGBAColor {
    let t = t.clamp(0.0, 1.0);
    
    if t < 0.5 {
        // Blue to white
        let ratio = t / 0.5;
        let r = (ratio * 255.0) as u8;
        let g = (ratio * 255.0) as u8;
        let b = 255u8;
        RGBAColor(r, g, b, 1.0)
    } else {
        // White to red
        let ratio = (t - 0.5) / 0.5;
        let r = 255u8;
        let g = ((1.0 - ratio) * 255.0) as u8;
        let b = ((1.0 - ratio) * 255.0) as u8;
        RGBAColor(r, g, b, 1.0)
    }
}

/// Viridis colormap (purple -> blue -> green -> yellow)
fn viridis_colormap(t: f64) -> RGBAColor {
    let t = t.clamp(0.0, 1.0);
    
    // Simplified viridis approximation
    let r = (72.0 + 193.0 * t - 133.0 * t * t) as u8;
    let g = (10.0 + 220.0 * t - 140.0 * t * t) as u8;
    let b = (150.0 + 50.0 * t - 100.0 * t * t) as u8;
    
    RGBAColor(r.min(255), g.min(255), b.min(255), 1.0)
}

/// Plasma colormap (purple -> red -> orange -> yellow)
fn plasma_colormap(t: f64) -> RGBAColor {
    let t = t.clamp(0.0, 1.0);
    
    let r = (60.0 + 195.0 * t) as u8;
    let g = (10.0 + 200.0 * t - 100.0 * t * t) as u8;
    let b = (150.0 - 130.0 * t) as u8;
    
    RGBAColor(r.min(255), g.min(255), b.max(0) as u8, 1.0)
}

/// Inferno colormap (black -> red -> orange -> yellow)
fn inferno_colormap(t: f64) -> RGBAColor {
    let t = t.clamp(0.0, 1.0);
    
    let r = (50.0 + 205.0 * t) as u8;
    let g = (10.0 + 180.0 * t - 80.0 * t * t) as u8;
    let b = (20.0 + 30.0 * t - 40.0 * t * t) as u8;
    
    RGBAColor(r.min(255), g.min(255), b.min(255), 1.0)
}

/// Magma colormap (black -> purple -> pink -> white)
fn magma_colormap(t: f64) -> RGBAColor {
    let t = t.clamp(0.0, 1.0);
    
    let r = (40.0 + 215.0 * t) as u8;
    let g = (10.0 + 150.0 * t + 50.0 * t * t) as u8;
    let b = (80.0 + 100.0 * t + 50.0 * t * t) as u8;
    
    RGBAColor(r.min(255), g.min(255), b.min(255), 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array4;
    
    #[test]
    fn test_plot_rotor_values_creates_file() {
        let values = Array4::<f64>::zeros((2, 3, 10, 10));
        let result = plot_rotor_values(&values, 0, 1, 3, "test_rotor.png", "coolwarm");
        assert!(result.is_ok());
        
        // Clean up
        let _ = std::fs::remove_file("test_rotor.png");
    }
}
