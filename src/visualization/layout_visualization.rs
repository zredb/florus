//! Layout visualization utilities for FLORUS (Simplified initial version)
//! 
//! This module provides basic functions to visualize wind farm layouts.
//! Corresponds to `floris/layout_visualization.py` in Python FLORIS.

use crate::{FlorisModel, Result};
use plotters::prelude::*;
use std::path::Path;

/// Plot turbine layout points
pub fn plot_turbine_points<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
    _turbine_indices: Option<&[usize]>,
    color: &str,
    marker_size: u32,
) -> Result<()> {
    let layout_x = fmodel.layout_x();
    let layout_y = fmodel.layout_y();
    
    if layout_x.is_empty() {
        anyhow::bail!("No turbines to plot");
    }
    
    // Determine plot bounds with some margin
    let x_min = layout_x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = layout_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = layout_y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = layout_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    let x_margin = (x_max - x_min) * 0.1;
    let y_margin = (y_max - y_min) * 0.1;
    
    // Create the plot
    let root = BitMapBackend::new(output_path.as_ref(), (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Wind Farm Layout", ("sans-serif", 25).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (x_min - x_margin)..(x_max + x_margin),
            (y_min - y_margin)..(y_max + y_margin),
        )?;
    
    chart.configure_mesh()
        .x_desc("X Position (m)")
        .y_desc("Y Position (m)")
        .draw()?;
    
    // Parse color
    let point_color = parse_color(color)?;
    
    // Plot turbine points
    for i in 0..layout_x.len() {
        chart.draw_series(std::iter::once(
            Circle::new((layout_x[i], layout_y[i]), marker_size, point_color)
        ))?;
    }
    
    Ok(())
}

/// Plot turbine labels on the layout
pub fn plot_turbine_labels<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
    turbine_names: Option<&[String]>,
    label_offset: Option<f64>,
    show_bbox: bool,
) -> Result<()> {
    let layout_x = fmodel.layout_x();
    let layout_y = fmodel.layout_y();
    let n_turbines = layout_x.len();
    
    if n_turbines == 0 {
        anyhow::bail!("No turbines to plot");
    }
    
    // Generate default turbine names if not provided
    let names: Vec<String> = if let Some(names) = turbine_names {
        if names.len() != n_turbines {
            anyhow::bail!(
                "Number of turbine names ({}) does not match number of turbines ({})",
                names.len(),
                n_turbines
            );
        }
        names.to_vec()
    } else {
        (0..n_turbines).map(|i| format!("{:03}", i)).collect()
    };
    
    // Calculate default offset based on rotor diameter
    let offset = label_offset.unwrap_or(10.0); // Default 10m offset
    
    // Determine plot bounds
    let x_min = layout_x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = layout_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = layout_y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = layout_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    let margin = (x_max - x_min).max(y_max - y_min) * 0.15;
    
    // Create the plot
    let root = BitMapBackend::new(output_path.as_ref(), (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Wind Farm Layout with Labels", ("sans-serif", 25).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (x_min - margin)..(x_max + margin),
            (y_min - margin)..(y_max + margin),
        )?;
    
    chart.configure_mesh()
        .x_desc("X Position (m)")
        .y_desc("Y Position (m)")
        .draw()?;
    
    // Plot turbine points and labels
    for i in 0..n_turbines {
        let x = layout_x[i];
        let y = layout_y[i];
        
        // Draw turbine point
        chart.draw_series(std::iter::once(
            Circle::new((x, y), 5, BLACK)
        ))?;
        
        // Add label
        let label_x = x + offset;
        let label_y = y + offset;
        
        // Use format! to create owned String that lives long enough
        let label_text = format!("{}", names[i]);
        
        chart.draw_series(std::iter::once(
            Text::new(
                label_text,
                (label_x, label_y),
                ("sans-serif", 15).into_font().color(&BLACK),
            )
        ))?;
    }
    
    Ok(())
}

/// Plot turbine rotors showing yaw angles
/// 
/// Draws lines representing turbine rotors oriented according to their yaw angles.
/// Corresponds to `plot_turbine_rotors()` in Python FLORIS.
pub fn plot_turbine_rotors<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
    color: &str,
    wd: Option<f64>,
) -> Result<()> {
    let layout_x = fmodel.layout_x();
    let layout_y = fmodel.layout_y();
    
    if layout_x.is_empty() {
        anyhow::bail!("No turbines to plot");
    }
    
    // Get yaw angles and wind direction
    let core = fmodel.core();
    let flow_field = fmodel.flow_field();
    
    let yaw_angles = &core.farm.yaw_angles;
    let ref_wd = wd.unwrap_or(flow_field.wind_directions[0]);
    
    // Get rotor diameters
    let rotor_diameters = &core.farm.rotor_diameters_sorted;
    
    // Calculate wind delta for rotation (270 - wd)
    let wind_delta_deg = 270.0 - ref_wd;
    
    // Determine plot bounds with margin for rotors
    let max_radius = rotor_diameters.iter().cloned().fold(f64::NEG_INFINITY, f64::max) / 2.0;
    let x_min = layout_x.iter().cloned().fold(f64::INFINITY, f64::min) - max_radius * 1.5;
    let x_max = layout_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + max_radius * 1.5;
    let y_min = layout_y.iter().cloned().fold(f64::INFINITY, f64::min) - max_radius * 1.5;
    let y_max = layout_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + max_radius * 1.5;
    
    // Create the plot
    let root = BitMapBackend::new(output_path.as_ref(), (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("Turbine Rotors (WD = {:.1}°)", ref_wd),
            ("sans-serif", 25).into_font()
        )
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(
            x_min..x_max,
            y_min..y_max,
        )?;
    
    chart.configure_mesh()
        .x_desc("X Position (m)")
        .y_desc("Y Position (m)")
        .draw()?;
    
    // Parse color
    let line_color = parse_color(color)?;
    
    // Draw rotor lines for each turbine (use first findex)
    let n_turbines = layout_x.len();
    for i in 0..n_turbines {
        let x = layout_x[i];
        let y = layout_y[i];
        let yaw = yaw_angles[(0, i)]; // First findex
        let d = rotor_diameters[(0, i)];
        let r = d / 2.0;
        
        // Rotate yaw angle to inertial frame for plotting
        // yaw_inertial = yaw - wind_delta
        let yaw_rad = (yaw - wind_delta_deg).to_radians();
        
        // Calculate rotor endpoints
        // The rotor is perpendicular to the yaw direction
        let x_0 = x + yaw_rad.sin() * r;
        let x_1 = x - yaw_rad.sin() * r;
        let y_0 = y - yaw_rad.cos() * r;
        let y_1 = y + yaw_rad.cos() * r;
        
        // Draw rotor line
        chart.draw_series(LineSeries::new(
            vec![(x_0, y_0), (x_1, y_1)],
            line_color,
        ))?;
    }
    
    Ok(())
}

/// Plot waking directions between turbines
/// 
/// Draws lines connecting turbines that can wake each other, with distance labels.
/// Corresponds to `plot_waking_directions()` in Python FLORIS.
pub fn plot_waking_directions<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
    limit_dist_d: Option<f64>,
    limit_num: Option<usize>,
) -> Result<()> {
    let layout_x = fmodel.layout_x();
    let layout_y = fmodel.layout_y();
    let n_turbines = layout_x.len();
    
    if n_turbines == 0 {
        anyhow::bail!("No turbines to plot");
    }
    
    // Get rotor diameter
    let core = fmodel.core();
    let d = core.farm.rotor_diameters_sorted[(0, 0)];
    
    // Calculate distances and angles between all turbine pairs
    let mut dists_m = vec![vec![0.0; n_turbines]; n_turbines];
    let mut angles_d = vec![vec![0.0; n_turbines]; n_turbines];
    
    for i in 0..n_turbines {
        for j in 0..n_turbines {
            let dx = layout_x[j] - layout_x[i];
            let dy = layout_y[j] - layout_y[i];
            dists_m[i][j] = (dx * dx + dy * dy).sqrt();
            
            // Calculate wake direction angle
            angles_d[i][j] = get_wake_direction(layout_x[i], layout_y[i], layout_x[j], layout_y[j]);
        }
    }
    
    // Apply distance limit if specified
    if let Some(limit_d) = limit_dist_d {
        let limit_m = limit_d * d;
        for i in 0..n_turbines {
            for j in 0..n_turbines {
                if dists_m[i][j] > limit_m {
                    dists_m[i][j] = f64::NAN;
                    angles_d[i][j] = f64::NAN;
                }
            }
        }
    }
    
    // Determine plot bounds
    let x_min = layout_x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = layout_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = layout_y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = layout_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    let margin = (x_max - x_min).max(y_max - y_min) * 0.1;
    
    // Create the plot
    let root = BitMapBackend::new(output_path.as_ref(), (1000, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Waking Directions", ("sans-serif", 25).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (x_min - margin)..(x_max + margin),
            (y_min - margin)..(y_max + margin),
        )?;
    
    chart.configure_mesh()
        .x_desc("X Position (m)")
        .y_desc("Y Position (m)")
        .draw()?;
    
    // Draw waking direction lines
    for i in 0..n_turbines {
        for j in 0..n_turbines {
            if !dists_m[i][j].is_nan() && dists_m[i][j] > 0.0 {
                // Check limit_num constraint
                if let Some(max_connections) = limit_num {
                    // Count valid connections for turbine i
                    let mut sorted_dists: Vec<_> = dists_m[i].clone();
                    sorted_dists.retain(|&d| !d.is_nan() && d > 0.0);
                    sorted_dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    
                    if dists_m[i][j] > sorted_dists.get(max_connections).copied().unwrap_or(f64::INFINITY) {
                        continue;
                    }
                }
                
                // Draw line
                chart.draw_series(LineSeries::new(
                    vec![(layout_x[i], layout_y[i]), (layout_x[j], layout_y[j])],
                    BLACK.mix(0.3),
                ))?;
                
                // Add distance label at midpoint
                let mid_x = (layout_x[i] + layout_x[j]) / 2.0;
                let mid_y = (layout_y[i] + layout_y[j]) / 2.0;
                let dist_d = dists_m[i][j] / d;
                
                chart.draw_series(std::iter::once(
                    Text::new(
                        format!("{:.1}D", dist_d),
                        (mid_x, mid_y),
                        ("sans-serif", 10).into_font().color(&BLACK),
                    )
                ))?;
            }
        }
    }
    
    // Plot turbine points on top
    for i in 0..n_turbines {
        chart.draw_series(std::iter::once(
            Circle::new((layout_x[i], layout_y[i]), 5, RED)
        ))?;
    }
    
    Ok(())
}

/// Calculate wake direction angle from turbine i to turbine j
/// 
/// Returns the wind direction (in degrees) at which turbine i would wake turbine j.
fn get_wake_direction(x_i: f64, y_i: f64, x_j: f64, y_j: f64) -> f64 {
    let dx = x_j - x_i;
    let dy = y_j - y_i;
    
    let angle_rad = dy.atan2(dx);
    let angle_deg = angle_rad.to_degrees();
    
    // Adjust for "from" direction and wrap to 0-360
    let wind_direction = (270.0 - angle_deg).rem_euclid(360.0);
    
    wind_direction
}

/// Shade a region defined by vertices
/// 
/// Fills a polygonal region defined by a set of vertices and optionally plots the vertices.
/// Corresponds to `shade_region()` in Python FLORIS.
pub fn shade_region<P: AsRef<Path>>(
    points: &[(f64, f64)],
    output_path: P,
    show_points: bool,
    region_color: &str,
    region_alpha: f64,
    point_color: &str,
) -> Result<()> {
    if points.len() < 3 {
        anyhow::bail!("Need at least 3 points to define a region");
    }
    
    // Determine plot bounds from points
    let x_min = points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let x_max = points.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    let y_min = points.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
    let y_max = points.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    
    let margin_x = (x_max - x_min) * 0.1;
    let margin_y = (y_max - y_min) * 0.1;
    
    // Create the plot
    let root = BitMapBackend::new(output_path.as_ref(), (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Shaded Region", ("sans-serif", 25).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (x_min - margin_x)..(x_max + margin_x),
            (y_min - margin_y)..(y_max + margin_y),
        )?;
    
    chart.configure_mesh()
        .x_desc("X Position (m)")
        .y_desc("Y Position (m)")
        .draw()?;
    
    // Parse colors
    let fill_color = parse_color(region_color)?;
    let pt_color = parse_color(point_color)?;
    
    // Create filled polygon with alpha
    let base_rgba = fill_color.to_rgba();
    let fill_color_with_alpha = RGBAColor(
        base_rgba.0,
        base_rgba.1,
        base_rgba.2,
        region_alpha,  // Already f64
    );
    
    // Draw filled polygon
    let polygon_points: Vec<_> = points.iter()
        .map(|&(x, y)| (x, y))
        .collect();
    
    chart.draw_series(std::iter::once(
        Polygon::new(polygon_points, fill_color_with_alpha)
    ))?;
    
    // Optionally draw vertex points
    if show_points {
        for &(x, y) in points {
            chart.draw_series(std::iter::once(
                Circle::new((x, y), 5, pt_color.filled())
            ))?;
        }
    }
    
    Ok(())
}

/// Plot farm terrain showing hub heights
/// 
/// Creates a visualization of turbine hub heights as a proxy for terrain.
/// Uses color-coded circles to represent different hub heights.
/// Corresponds to `plot_farm_terrain()` in Python FLORIS.
pub fn plot_farm_terrain<P: AsRef<Path>>(
    fmodel: &FlorisModel,
    output_path: P,
) -> Result<()> {
    let layout_x = fmodel.layout_x();
    let layout_y = fmodel.layout_y();
    
    if layout_x.is_empty() {
        anyhow::bail!("No turbines to plot");
    }
    
    // Get hub heights from farm
    let core = fmodel.core();
    let hub_heights = &core.farm.hub_heights;
    
    if hub_heights.len() != layout_x.len() {
        anyhow::bail!("Hub heights count doesn't match turbine count");
    }
    
    // Find min/max hub heights for color scaling
    let min_hh = hub_heights.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_hh = hub_heights.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let hh_range = max_hh - min_hh;
    
    // Determine plot bounds
    let x_min = layout_x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = layout_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = layout_y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = layout_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    let margin_x = (x_max - x_min) * 0.15;
    let margin_y = (y_max - y_min) * 0.15;
    
    // Create the plot
    let root = BitMapBackend::new(output_path.as_ref(), (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("Farm Terrain (Hub Heights: {:.0}-{:.0} m)", min_hh, max_hh),
            ("sans-serif", 20).into_font()
        )
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (x_min - margin_x)..(x_max + margin_x),
            (y_min - margin_y)..(y_max + margin_y),
        )?;
    
    chart.configure_mesh()
        .x_desc("X Position (m)")
        .y_desc("Y Position (m)")
        .draw()?;
    
    // Plot turbines with color-coded hub heights
    // Use RdBu_r colormap: low = red, mid = white, high = blue
    for i in 0..layout_x.len() {
        let x = layout_x[i];
        let y = layout_y[i];
        let hh = hub_heights[i];
        
        // Normalize hub height to 0-1 range
        let normalized = if hh_range > 0.0 {
            (hh - min_hh) / hh_range
        } else {
            0.5
        };
        
        // Map to RdBu_r colormap (reversed: low=red, high=blue)
        let color = if normalized < 0.5 {
            // Red to white
            let ratio = normalized / 0.5;
            let r = 255u8;
            let g = (ratio * 255.0) as u8;
            let b = (ratio * 255.0) as u8;
            RGBAColor(r, g, b, 1.0)
        } else {
            // White to blue
            let ratio = (normalized - 0.5) / 0.5;
            let r = ((1.0 - ratio) * 255.0) as u8;
            let g = ((1.0 - ratio) * 255.0) as u8;
            let b = 255u8;
            RGBAColor(r, g, b, 1.0)
        };
        
        // Draw circle with size proportional to hub height
        let radius = 8.0 + (normalized * 7.0); // 8-15 pixels
        
        chart.draw_series(std::iter::once(
            Circle::new((x, y), radius as u32, color.filled())
        ))?;
        
        // Add hub height label
        chart.draw_series(std::iter::once(
            Text::new(
                format!("{:.0}", hh),
                (x, y - 20.0),
                ("sans-serif", 12).into_font().color(&BLACK),
            )
        ))?;
    }
    
    Ok(())
}

/// Helper function to parse color strings
fn parse_color(color: &str) -> Result<plotters::style::RGBColor> {
    match color.to_lowercase().as_str() {
        "black" | "k" => Ok(BLACK),
        "white" | "w" => Ok(WHITE),
        "red" | "r" => Ok(RED),
        "green" | "g" => Ok(GREEN),
        "blue" | "b" => Ok(BLUE),
        "yellow" | "y" => Ok(YELLOW),
        "cyan" | "c" => Ok(CYAN),
        "magenta" | "m" => Ok(MAGENTA),
        _ => {
            eprintln!("Warning: Unknown color '{}', using black", color);
            Ok(BLACK)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse_color() {
        assert!(parse_color("black").is_ok());
        assert!(parse_color("red").is_ok());
        assert!(parse_color("unknown").is_ok()); // Should default to black
    }
}
