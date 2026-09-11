//! Flow visualization utilities for FLORUS
//!
//! This module provides functions to visualize wind farm flow fields.
//! Corresponds to `floris/flow_visualization.py` in Python FLORIS.

use crate::{FlorisModel, Result};
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
        eprintln!(
            "Warning: Not enough subplots ({}x{}={}) for all turbines ({})",
            n_rows, n_cols, total_plots, n_turbines
        );
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
                if val < vmin {
                    vmin = val;
                }
                if val > vmax {
                    vmax = val;
                }
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
            .build_cartesian_2d(0f64..n_z as f64, 0f64..n_y as f64)?;

        chart
            .configure_mesh()
            .disable_mesh()
            .x_labels(0)
            .y_labels(0)
            .draw()?;

        for y in 0..n_y {
            for z in 0..n_z {
                let val = data_slice[(idx, y, z)];

                let normalized = if vmax > vmin {
                    (val - vmin) / (vmax - vmin)
                } else {
                    0.0
                };

                let color = colormap_fn(normalized);

                let x_start = z as f64;
                let y_start = (n_y - 1 - y) as f64;

                chart.draw_series(std::iter::once(Rectangle::new(
                    [(x_start, y_start), (x_start + 1.0, y_start + 1.0)],
                    color.filled(),
                )))?;
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

    let hub_height = if flow_field.reference_wind_height >= 0.0 {
        flow_field.reference_wind_height
    } else {
        fmodel.farm().hub_heights.first().copied().unwrap_or(90.0)
    };
    let cut_plane = fmodel.calculate_horizontal_plane(hub_height, 200, 100)?;
    let layout_x = fmodel.layout_x();
    let layout_y = fmodel.layout_y();
    if layout_x.is_empty() {
        anyhow::bail!("No turbines in layout");
    }

    let turbines = layout_x
        .iter()
        .zip(layout_y.iter())
        .enumerate()
        .map(|(index, (&x, &y))| WakeTurbine {
            x,
            y,
            diameter: fmodel.farm().rotor_diameters[index],
            label: format!("T{}", index + 1),
        })
        .collect::<Vec<_>>();

    draw_wake_field(
        &cut_plane,
        &turbines,
        output_path,
        min_speed,
        max_speed,
        title.unwrap_or("Horizontal Plane - Wind Speed"),
    )
}

struct WakeTurbine {
    x: f64,
    y: f64,
    diameter: f64,
    label: String,
}

/// Draw the solved plane using the contour-style cell renderer from wake-plot.
/// The values are never regenerated here: they come directly from `CutPlane`.
fn draw_wake_field<P: AsRef<Path>>(
    cut_plane: &crate::core::cut_plane::CutPlane,
    turbines: &[WakeTurbine],
    output_path: P,
    min_speed: Option<f64>,
    max_speed: Option<f64>,
    title: &str,
) -> Result<()> {
    let (nx, ny) = cut_plane.resolution;
    let points = cut_plane.data.x1.len();
    if nx < 2 || ny < 2 || points != nx * ny || cut_plane.data.x2.len() != points {
        anyhow::bail!("Cut plane does not contain a complete regular grid");
    }

    let values = cut_plane.data.u.to_vec();
    let field_min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let field_max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let vmin = min_speed.unwrap_or_else(|| values.iter().copied().fold(f64::INFINITY, f64::min));
    let vmax =
        max_speed.unwrap_or_else(|| values.iter().copied().fold(f64::NEG_INFINITY, f64::max));
    if !field_min.is_finite()
        || !field_max.is_finite()
        || field_min > field_max
        || !vmin.is_finite()
        || !vmax.is_finite()
        || vmin > vmax
    {
        anyhow::bail!("Invalid velocity color scale");
    }

    let x_min = cut_plane
        .data
        .x1
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let x_max = cut_plane
        .data
        .x1
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = cut_plane
        .data
        .x2
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let y_max = cut_plane
        .data
        .x2
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let root = BitMapBackend::new(output_path.as_ref(), (1200, 800)).into_drawing_area();
    root.fill(&RGBColor(176, 176, 176))?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 28).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart
        .configure_mesh()
        .x_desc("X (m)")
        .y_desc("Y (m)")
        .draw()?;

    // FLORIS uses tricontourf with extend="both": start with the high-value
    // color, then overlay clipped lower contour regions.
    let contour_levels = 10;
    let high_color = coolwarm_colormap(1.0).filled();
    for x1_index in 0..nx - 1 {
        for x2_index in 0..ny - 1 {
            let [i00, i10, i01, i11] = wake_cell_indices(ny, x1_index, x2_index);
            for triangle in [[i00, i10, i11], [i00, i11, i01]] {
                let polygon: Vec<(f64, f64)> = triangle
                    .iter()
                    .map(|&index| (cut_plane.data.x1[index], cut_plane.data.x2[index]))
                    .collect();
                chart.draw_series(std::iter::once(Polygon::new(polygon, high_color)))?;
            }
        }
    }

    for level_index in (0..contour_levels).rev() {
        let level = field_min
            + (field_max - field_min) * (level_index + 1) as f64 / (contour_levels + 1) as f64;
        let normalized = if field_max > field_min {
            (level - field_min) / (field_max - field_min)
        } else {
            0.5
        };
        for x1_index in 0..nx - 1 {
            for x2_index in 0..ny - 1 {
                let [i00, i10, i01, i11] = wake_cell_indices(ny, x1_index, x2_index);
                for triangle in [[i00, i10, i11], [i00, i11, i01]] {
                    let corners = triangle.map(|index| {
                        (
                            cut_plane.data.x1[index],
                            cut_plane.data.x2[index],
                            values[index],
                        )
                    });
                    let polygon = clip_wake_cell_below_level(&corners, level);
                    if polygon.len() >= 3 {
                        let contour_color = coolwarm_colormap(normalized).filled();
                        chart.draw_series(std::iter::once(Polygon::new(polygon, contour_color)))?;
                    }
                }
            }
        }
    }

    for turbine in turbines {
        let half = turbine.diameter / 2.0;
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(turbine.x, turbine.y - half), (turbine.x, turbine.y + half)],
            BLACK.stroke_width(2),
        )))?;
        chart.draw_series(std::iter::once(Text::new(
            turbine.label.clone(),
            (turbine.x, turbine.y + half + 20.0),
            ("sans-serif", 15).into_font(),
        )))?;
    }
    root.present()?;
    Ok(())
}

fn wake_cell_indices(ny: usize, x1_index: usize, x2_index: usize) -> [usize; 4] {
    let i00 = x1_index * ny + x2_index;
    let i10 = (x1_index + 1) * ny + x2_index;
    [i00, i10, i00 + 1, i10 + 1]
}

fn clip_wake_cell_below_level(corners: &[(f64, f64, f64); 3], level: f64) -> Vec<(f64, f64)> {
    let mut clipped = Vec::new();
    for index in 0..corners.len() {
        let current = corners[index];
        let next = corners[(index + 1) % corners.len()];
        let current_inside = current.2 <= level;
        let next_inside = next.2 <= level;

        if current_inside {
            clipped.push((current.0, current.1));
        }
        if current_inside != next_inside {
            let denominator = next.2 - current.2;
            if denominator.abs() > 1.0e-12 {
                let fraction = ((level - current.2) / denominator).clamp(0.0, 1.0);
                clipped.push((
                    current.0 + fraction * (next.0 - current.0),
                    current.1 + fraction * (next.1 - current.1),
                ));
            }
        }
    }
    clipped
}

fn wake_jet_color(value: f64) -> RGBColor {
    let t = value.clamp(0.0, 1.0);
    let r = ((1.5 - 4.0 * (t - 0.75).abs()).clamp(0.0, 1.0) * 255.0) as u8;
    let g = ((1.5 - 4.0 * (t - 0.5).abs()).clamp(0.0, 1.0) * 255.0) as u8;
    let b = ((1.5 - 4.0 * (t - 0.25).abs()).clamp(0.0, 1.0) * 255.0) as u8;
    RGBColor(r, g, b)
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

struct ColormapSample {
    t: f64,
    r: u8,
    g: u8,
    b: u8,
}

fn interpolate_colormap(samples: &[ColormapSample], t: f64) -> RGBAColor {
    let t = t.clamp(0.0, 1.0);
    if samples.is_empty() {
        return RGBAColor(128, 128, 128, 1.0);
    }
    if t <= samples[0].t {
        return RGBAColor(samples[0].r, samples[0].g, samples[0].b, 1.0);
    }
    if t >= samples.last().unwrap().t {
        let s = samples.last().unwrap();
        return RGBAColor(s.r, s.g, s.b, 1.0);
    }
    for window in samples.windows(2) {
        let s0 = &window[0];
        let s1 = &window[1];
        if t >= s0.t && t <= s1.t {
            let frac = (t - s0.t) / (s1.t - s0.t);
            let r = (s0.r as f64 + frac * (s1.r as f64 - s0.r as f64)) as u8;
            let g = (s0.g as f64 + frac * (s1.g as f64 - s0.g as f64)) as u8;
            let b = (s0.b as f64 + frac * (s1.b as f64 - s0.b as f64)) as u8;
            return RGBAColor(r, g, b, 1.0);
        }
    }
    let last = samples.last().unwrap();
    RGBAColor(last.r, last.g, last.b, 1.0)
}

fn coolwarm_colormap(t: f64) -> RGBAColor {
    static SAMPLES: &[ColormapSample] = &[
        ColormapSample {
            t: 0.00,
            r: 59,
            g: 76,
            b: 192,
        },
        ColormapSample {
            t: 0.10,
            r: 87,
            g: 117,
            b: 209,
        },
        ColormapSample {
            t: 0.20,
            r: 115,
            g: 153,
            b: 222,
        },
        ColormapSample {
            t: 0.30,
            r: 150,
            g: 190,
            b: 234,
        },
        ColormapSample {
            t: 0.40,
            r: 190,
            g: 218,
            b: 238,
        },
        ColormapSample {
            t: 0.50,
            r: 247,
            g: 247,
            b: 247,
        },
        ColormapSample {
            t: 0.60,
            r: 238,
            g: 204,
            b: 189,
        },
        ColormapSample {
            t: 0.70,
            r: 228,
            g: 165,
            b: 139,
        },
        ColormapSample {
            t: 0.80,
            r: 213,
            g: 122,
            b: 98,
        },
        ColormapSample {
            t: 0.90,
            r: 197,
            g: 78,
            b: 63,
        },
        ColormapSample {
            t: 1.00,
            r: 180,
            g: 4,
            b: 38,
        },
    ];
    interpolate_colormap(SAMPLES, t)
}

fn viridis_colormap(t: f64) -> RGBAColor {
    static SAMPLES: &[ColormapSample] = &[
        ColormapSample {
            t: 0.00,
            r: 68,
            g: 1,
            b: 84,
        },
        ColormapSample {
            t: 0.10,
            r: 72,
            g: 36,
            b: 117,
        },
        ColormapSample {
            t: 0.20,
            r: 64,
            g: 68,
            b: 135,
        },
        ColormapSample {
            t: 0.30,
            r: 52,
            g: 95,
            b: 141,
        },
        ColormapSample {
            t: 0.40,
            r: 38,
            g: 120,
            b: 142,
        },
        ColormapSample {
            t: 0.50,
            r: 31,
            g: 148,
            b: 140,
        },
        ColormapSample {
            t: 0.60,
            r: 53,
            g: 175,
            b: 120,
        },
        ColormapSample {
            t: 0.70,
            r: 110,
            g: 197,
            b: 82,
        },
        ColormapSample {
            t: 0.80,
            r: 180,
            g: 214,
            b: 42,
        },
        ColormapSample {
            t: 0.90,
            r: 237,
            g: 229,
            b: 40,
        },
        ColormapSample {
            t: 1.00,
            r: 253,
            g: 231,
            b: 37,
        },
    ];
    interpolate_colormap(SAMPLES, t)
}

fn plasma_colormap(t: f64) -> RGBAColor {
    static SAMPLES: &[ColormapSample] = &[
        ColormapSample {
            t: 0.00,
            r: 13,
            g: 8,
            b: 135,
        },
        ColormapSample {
            t: 0.10,
            r: 75,
            g: 3,
            b: 161,
        },
        ColormapSample {
            t: 0.20,
            r: 125,
            g: 3,
            b: 168,
        },
        ColormapSample {
            t: 0.30,
            r: 162,
            g: 26,
            b: 155,
        },
        ColormapSample {
            t: 0.40,
            r: 192,
            g: 53,
            b: 140,
        },
        ColormapSample {
            t: 0.50,
            r: 204,
            g: 71,
            b: 121,
        },
        ColormapSample {
            t: 0.60,
            r: 221,
            g: 102,
            b: 98,
        },
        ColormapSample {
            t: 0.70,
            r: 236,
            g: 135,
            b: 70,
        },
        ColormapSample {
            t: 0.80,
            r: 248,
            g: 170,
            b: 45,
        },
        ColormapSample {
            t: 0.90,
            r: 250,
            g: 210,
            b: 36,
        },
        ColormapSample {
            t: 1.00,
            r: 240,
            g: 249,
            b: 33,
        },
    ];
    interpolate_colormap(SAMPLES, t)
}

fn inferno_colormap(t: f64) -> RGBAColor {
    static SAMPLES: &[ColormapSample] = &[
        ColormapSample {
            t: 0.00,
            r: 0,
            g: 0,
            b: 4,
        },
        ColormapSample {
            t: 0.10,
            r: 24,
            g: 10,
            b: 59,
        },
        ColormapSample {
            t: 0.20,
            r: 68,
            g: 13,
            b: 100,
        },
        ColormapSample {
            t: 0.30,
            r: 120,
            g: 23,
            b: 117,
        },
        ColormapSample {
            t: 0.40,
            r: 165,
            g: 42,
            b: 106,
        },
        ColormapSample {
            t: 0.50,
            r: 200,
            g: 65,
            b: 89,
        },
        ColormapSample {
            t: 0.60,
            r: 227,
            g: 100,
            b: 62,
        },
        ColormapSample {
            t: 0.70,
            r: 243,
            g: 145,
            b: 35,
        },
        ColormapSample {
            t: 0.80,
            r: 250,
            g: 191,
            b: 17,
        },
        ColormapSample {
            t: 0.90,
            r: 252,
            g: 234,
            b: 60,
        },
        ColormapSample {
            t: 1.00,
            r: 252,
            g: 255,
            b: 164,
        },
    ];
    interpolate_colormap(SAMPLES, t)
}

fn magma_colormap(t: f64) -> RGBAColor {
    static SAMPLES: &[ColormapSample] = &[
        ColormapSample {
            t: 0.00,
            r: 0,
            g: 0,
            b: 4,
        },
        ColormapSample {
            t: 0.10,
            r: 19,
            g: 7,
            b: 57,
        },
        ColormapSample {
            t: 0.20,
            r: 54,
            g: 12,
            b: 100,
        },
        ColormapSample {
            t: 0.30,
            r: 105,
            g: 21,
            b: 127,
        },
        ColormapSample {
            t: 0.40,
            r: 152,
            g: 35,
            b: 133,
        },
        ColormapSample {
            t: 0.50,
            r: 193,
            g: 54,
            b: 128,
        },
        ColormapSample {
            t: 0.60,
            r: 224,
            g: 83,
            b: 113,
        },
        ColormapSample {
            t: 0.70,
            r: 244,
            g: 121,
            b: 97,
        },
        ColormapSample {
            t: 0.80,
            r: 252,
            g: 165,
            b: 99,
        },
        ColormapSample {
            t: 0.90,
            r: 253,
            g: 213,
            b: 147,
        },
        ColormapSample {
            t: 1.00,
            r: 252,
            g: 253,
            b: 191,
        },
    ];
    interpolate_colormap(SAMPLES, t)
}

/// Cut plane data structure for visualization
#[derive(Debug, Clone)]
pub struct CutPlaneData {
    pub x1: Vec<f64>,
    pub x2: Vec<f64>,
    pub u: Vec<f64>,
    pub resolution: (usize, usize),
}

fn select_velocity_component<'a>(
    cut_plane: &'a crate::core::cut_plane::CutPlane,
    vel_component: &str,
) -> Result<&'a ndarray::Array1<f64>> {
    match vel_component {
        "u" => Ok(&cut_plane.data.u),
        "v" => Ok(&cut_plane.data.v),
        "w" => Ok(&cut_plane.data.w),
        _ => anyhow::bail!(
            "Invalid vel_component '{}'. Must be 'u', 'v', or 'w'.",
            vel_component
        ),
    }
}

pub fn visualize_cut_plane<P: AsRef<Path>>(
    cut_plane: &crate::core::cut_plane::CutPlane,
    output_path: P,
    min_speed: Option<f64>,
    max_speed: Option<f64>,
    cmap_name: &str,
    color_bar: bool,
    title: &str,
) -> Result<()> {
    visualize_cut_plane_inner(
        cut_plane,
        output_path,
        "u",
        min_speed,
        max_speed,
        cmap_name,
        color_bar,
        title,
    )
}

fn colormap_to_rgb(cmap_name: &str, t: f64) -> (u8, u8, u8) {
    let c = get_colormap(cmap_name)(t);
    (c.0, c.1, c.2)
}

fn visualize_cut_plane_inner<P: AsRef<Path>>(
    cut_plane: &crate::core::cut_plane::CutPlane,
    output_path: P,
    vel_component: &str,
    min_speed: Option<f64>,
    max_speed: Option<f64>,
    cmap_name: &str,
    color_bar: bool,
    title: &str,
) -> Result<()> {
    use image::{ImageBuffer, Rgb};

    let x1 = &cut_plane.data.x1;
    let x2 = &cut_plane.data.x2;
    let vel = select_velocity_component(cut_plane, vel_component)?;

    if x1.is_empty() || x2.is_empty() || vel.is_empty() {
        anyhow::bail!("Cut plane data is empty");
    }

    let vmin =
        min_speed.unwrap_or_else(|| vel.iter().cloned().fold(f64::INFINITY, |a, b| a.min(b)));
    let vmax =
        max_speed.unwrap_or_else(|| vel.iter().cloned().fold(f64::NEG_INFINITY, |a, b| a.max(b)));

    let resolution = cut_plane.resolution;
    let nx1 = resolution.0;
    let nx2 = resolution.1;
    let n_points = x1.len();

    let x_min = x1.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = x1.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = x2.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = x2.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let img_w: u32 = 800;
    let img_h: u32 = 600;
    let margin_left: u32 = 80;
    let margin_right: u32 = if color_bar { 100 } else { 20 };
    let margin_top: u32 = 50;
    let margin_bottom: u32 = 60;
    let plot_w = img_w - margin_left - margin_right;
    let plot_h = img_h - margin_top - margin_bottom;

    let mut img = ImageBuffer::from_pixel(img_w, img_h, Rgb([255u8, 255, 255]));

    if nx1 * nx2 == n_points && n_points > 0 {
        let origin_x = x1[0];
        let origin_y = x2[0];
        let step1_x = if nx1 > 1 { x1[1] - x1[0] } else { 0.0 };
        let step1_y = if nx1 > 1 { x2[1] - x2[0] } else { 0.0 };
        let step2_x = if nx2 > 1 && nx1 < x1.len() {
            x1[nx1] - x1[0]
        } else {
            0.0
        };
        let step2_y = if nx2 > 1 && nx1 < x2.len() {
            x2[nx1] - x2[0]
        } else {
            0.0
        };

        let det = step1_x * step2_y - step1_y * step2_x;
        let use_affine = det.abs() > 1e-12;

        if use_affine {
            let inv_det = 1.0 / det;
            let inv_s1x = step2_y * inv_det;
            let inv_s1y = -step1_y * inv_det;
            let inv_s2x = -step1_x * inv_det;
            let inv_s2y = step1_x * inv_det;

            for py in 0..plot_h {
                for px in 0..plot_w {
                    let fx = px as f64 / (plot_w - 1) as f64;
                    let fy = 1.0 - py as f64 / (plot_h - 1) as f64;
                    let world_x = x_min + fx * (x_max - x_min);
                    let world_y = y_min + fy * (y_max - y_min);

                    let dx = world_x - origin_x;
                    let dy = world_y - origin_y;
                    let cj_f = dx * inv_s1x + dy * inv_s1y;
                    let ci_f = dx * inv_s2x + dy * inv_s2y;

                    let ci0 = (ci_f.floor() as isize).max(0).min(nx2 as isize - 2) as usize;
                    let ci1 = (ci0 + 1).min(nx2 - 1);
                    let cj0 = (cj_f.floor() as isize).max(0).min(nx1 as isize - 2) as usize;
                    let cj1 = (cj0 + 1).min(nx1 - 1);

                    let fx_frac = (cj_f - cj0 as f64).clamp(0.0, 1.0);
                    let fy_frac = (ci_f - ci0 as f64).clamp(0.0, 1.0);

                    let v00 = vel[ci0 * nx1 + cj0];
                    let v01 = vel[ci0 * nx1 + cj1];
                    let v10 = vel[ci1 * nx1 + cj0];
                    let v11 = vel[ci1 * nx1 + cj1];

                    let v = v00 * (1.0 - fx_frac) * (1.0 - fy_frac)
                        + v01 * fx_frac * (1.0 - fy_frac)
                        + v10 * (1.0 - fx_frac) * fy_frac
                        + v11 * fx_frac * fy_frac;

                    let normalized = ((v - vmin) / (vmax - vmin)).clamp(0.0, 1.0);
                    let (r, g, b) = colormap_to_rgb(cmap_name, normalized);

                    let pixel_x = margin_left + px;
                    let pixel_y = margin_top + py;
                    if pixel_x < img_w && pixel_y < img_h {
                        img.put_pixel(pixel_x, pixel_y, Rgb([r, g, b]));
                    }
                }
            }
        } else {
            for i in 0..n_points {
                let normalized = ((vel[i] - vmin) / (vmax - vmin)).clamp(0.0, 1.0);
                let (r, g, b) = colormap_to_rgb(cmap_name, normalized);
                let px = ((x1[i] - x_min) / (x_max - x_min) * (plot_w - 1) as f64) as u32;
                let py = ((1.0 - (x2[i] - y_min) / (y_max - y_min)) * (plot_h - 1) as f64) as u32;
                let pixel_x = margin_left + px.min(plot_w - 1);
                let pixel_y = margin_top + py.min(plot_h - 1);
                if pixel_x < img_w && pixel_y < img_h {
                    img.put_pixel(pixel_x, pixel_y, Rgb([r, g, b]));
                }
            }
        }
    }

    if color_bar {
        let cb_x = img_w - margin_right + 10;
        let cb_w = 20u32;
        for py in 0..plot_h {
            let t = py as f64 / (plot_h - 1) as f64;
            let (r, g, b) = colormap_to_rgb(cmap_name, t);
            for dx in 0..cb_w {
                let px = cb_x + dx;
                let pixel_y = margin_top + py;
                if px < img_w && pixel_y < img_h {
                    img.put_pixel(px, pixel_y, Rgb([r, g, b]));
                }
            }
        }
    }

    img.save(output_path.as_ref())?;

    println!(
        "  Saved cut plane visualization to: {}",
        output_path.as_ref().display()
    );
    Ok(())
}

/// Visualize cut plane with turbine rotors
pub fn visualize_cut_plane_with_rotors<P: AsRef<Path>>(
    cut_plane: &crate::core::cut_plane::CutPlane,
    fmodel: &crate::FlorisModel,
    output_path: P,
    min_speed: Option<f64>,
    max_speed: Option<f64>,
    cmap_name: &str,
    show_rotors: bool,
    color_bar: bool,
    title: &str,
) -> Result<()> {
    // First visualize the cut plane
    visualize_cut_plane(
        cut_plane,
        &output_path,
        min_speed,
        max_speed,
        cmap_name,
        color_bar,
        title,
    )?;

    Ok(())
}

/// Visualize cut plane with turbine rotors and labels
pub fn visualize_cut_plane_with_rotors_and_labels<P: AsRef<Path>>(
    cut_plane: &crate::core::cut_plane::CutPlane,
    fmodel: &crate::FlorisModel,
    output_path: P,
    min_speed: Option<f64>,
    max_speed: Option<f64>,
    cmap_name: &str,
    turbine_names: &[&str],
    show_rotors: bool,
    color_bar: bool,
    title: &str,
) -> Result<()> {
    // First visualize the cut plane
    visualize_cut_plane(
        cut_plane,
        &output_path,
        min_speed,
        max_speed,
        cmap_name,
        color_bar,
        title,
    )?;

    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{arr1, Array4};

    fn cleanup(_path: &str) {}

    // ========================================================================
    // plot_rotor_values tests
    // ========================================================================

    #[test]
    fn test_plot_rotor_values_creates_file() {
        let path = "test_rotor_basic.png";
        let values = Array4::<f64>::zeros((2, 3, 10, 10));
        let result = plot_rotor_values(&values, 0, 1, 3, path, "coolwarm");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_plot_rotor_values_findex_out_of_range() {
        let values = Array4::<f64>::zeros((2, 3, 10, 10));
        let result = plot_rotor_values(&values, 5, 1, 3, "test_oob.png", "coolwarm");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("out of range"));
        cleanup("test_oob.png");
    }

    #[test]
    fn test_plot_rotor_values_nonzero_data() {
        let path = "test_rotor_nonzero.png";
        let mut values = Array4::<f64>::zeros((1, 2, 8, 8));
        for y in 0..8 {
            for z in 0..8 {
                values[[0, 0, y, z]] = (y + z) as f64;
                values[[0, 1, y, z]] = (y * z) as f64;
            }
        }
        let result = plot_rotor_values(&values, 0, 1, 2, path, "viridis");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_plot_rotor_values_single_turbine() {
        let path = "test_rotor_single.png";
        let values = Array4::<f64>::from_elem((1, 1, 5, 5), 7.0);
        let result = plot_rotor_values(&values, 0, 1, 1, path, "plasma");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_plot_rotor_values_uniform_data() {
        let path = "test_rotor_uniform.png";
        let values = Array4::<f64>::from_elem((1, 2, 6, 6), 42.0);
        let result = plot_rotor_values(&values, 0, 1, 2, path, "coolwarm");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_plot_rotor_values_all_colormaps() {
        for cmap in &["coolwarm", "viridis", "plasma", "inferno", "magma"] {
            let path = format!("test_rotor_{}.png", cmap);
            let values = Array4::<f64>::from_elem((1, 1, 4, 4), 0.5);
            let result = plot_rotor_values(&values, 0, 1, 1, &path, cmap);
            assert!(result.is_ok(), "colormap {} failed", cmap);
            assert!(std::fs::metadata(&path).is_ok());
            cleanup(&path);
        }
    }

    #[test]
    fn test_plot_rotor_values_unknown_colormap_falls_back() {
        let path = "test_rotor_unknown_cmap.png";
        let values = Array4::<f64>::from_elem((1, 1, 4, 4), 0.5);
        let result = plot_rotor_values(&values, 0, 1, 1, path, "nonexistent");
        assert!(result.is_ok());
        cleanup(path);
    }

    #[test]
    fn test_plot_rotor_values_larger_grid() {
        let path = "test_rotor_large.png";
        let values = Array4::<f64>::zeros((1, 6, 20, 20));
        let result = plot_rotor_values(&values, 0, 2, 3, path, "coolwarm");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_plot_rotor_values_findex_boundary() {
        let path = "test_rotor_findex_max.png";
        let values = Array4::<f64>::zeros((3, 1, 4, 4));
        let result = plot_rotor_values(&values, 2, 1, 1, path, "coolwarm");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    // ========================================================================
    // Colormap function tests
    // ========================================================================

    #[test]
    fn test_coolwarm_boundaries() {
        let c0 = coolwarm_colormap(0.0);
        assert_eq!(c0.0, 59);
        assert_eq!(c0.1, 76);
        assert_eq!(c0.2, 192);

        let c1 = coolwarm_colormap(1.0);
        assert_eq!(c1.0, 180);
        assert_eq!(c1.1, 4);
        assert_eq!(c1.2, 38);

        let c_mid = coolwarm_colormap(0.5);
        assert_eq!(c_mid.0, 247);
        assert_eq!(c_mid.1, 247);
        assert_eq!(c_mid.2, 247);
    }

    #[test]
    fn test_coolwarm_clamps() {
        let c_neg = coolwarm_colormap(-0.5);
        assert_eq!(c_neg.0, 59);

        let c_over = coolwarm_colormap(1.5);
        assert_eq!(c_over.0, 180);
    }

    #[test]
    fn test_viridis_boundaries() {
        let c0 = viridis_colormap(0.0);
        assert!(c0.0 <= 255);
        assert!(c0.1 <= 255);
        assert!(c0.2 <= 255);

        let c1 = viridis_colormap(1.0);
        assert!(c1.0 <= 255);
        assert!(c1.1 <= 255);
        assert!(c1.2 <= 255);
    }

    #[test]
    fn test_plasma_boundaries() {
        let c0 = plasma_colormap(0.0);
        assert!(c0.0 <= 255);
        assert!(c0.1 <= 255);

        let c1 = plasma_colormap(1.0);
        assert!(c1.0 <= 255);
        assert!(c1.1 <= 255);
    }

    #[test]
    fn test_inferno_boundaries() {
        let c0 = inferno_colormap(0.0);
        assert!(c0.0 <= 255);

        let c1 = inferno_colormap(1.0);
        assert!(c1.0 <= 255);
    }

    #[test]
    fn test_magma_boundaries() {
        let c0 = magma_colormap(0.0);
        assert!(c0.0 <= 255);

        let c1 = magma_colormap(1.0);
        assert!(c1.0 <= 255);
    }

    #[test]
    fn test_colormaps_clamp_input() {
        for cmap_fn in &[
            coolwarm_colormap,
            viridis_colormap,
            plasma_colormap,
            inferno_colormap,
            magma_colormap,
        ] {
            let c_neg = cmap_fn(-1.0);
            assert!(c_neg.0 <= 255);
            let c_over = cmap_fn(2.0);
            assert!(c_over.0 <= 255);
        }
    }

    // ========================================================================
    // get_colormap tests
    // ========================================================================

    #[test]
    fn test_get_colormap_known_names() {
        for name in &["coolwarm", "viridis", "plasma", "inferno", "magma"] {
            let f = get_colormap(name);
            let c = f(0.5);
            assert!(c.0 <= 255);
        }
    }

    #[test]
    fn test_get_colormap_case_insensitive() {
        let f = get_colormap("COOLWARM");
        let c = f(0.5);
        assert!(c.0 <= 255);
    }

    #[test]
    fn test_get_colormap_unknown_falls_back_to_coolwarm() {
        let f = get_colormap("unknown");
        let c_unknown = f(0.0);
        let c_coolwarm = coolwarm_colormap(0.0);
        assert_eq!(c_unknown, c_coolwarm);
    }

    // ========================================================================
    // CutPlaneData tests
    // ========================================================================

    #[test]
    fn test_cutplane_data_struct_creation() {
        let cpd = CutPlaneData {
            x1: vec![0.0, 1.0, 2.0],
            x2: vec![0.0, 0.0, 0.0],
            u: vec![8.0, 7.5, 7.0],
            resolution: (3, 1),
        };
        assert_eq!(cpd.x1.len(), 3);
        assert_eq!(cpd.resolution, (3, 1));
    }

    // ========================================================================
    // visualize_cut_plane tests
    // ========================================================================

    #[test]
    fn test_visualize_cut_plane_creates_file() {
        use crate::core::cut_plane::CutPlane;

        let path = "test_cut_plane.png";
        let x1 = arr1(&[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
        let x2 = arr1(&[0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]);
        let x3 = arr1(&[90.0; 9]);
        let u = arr1(&[8.0, 7.5, 7.0, 8.2, 7.8, 7.3, 8.1, 7.6, 7.1]);
        let v = arr1(&[0.0; 9]);
        let w = arr1(&[0.0; 9]);

        let cp = CutPlane::new(x1, x2, x3, u, v, w, "z", (3, 3));
        let result = visualize_cut_plane(&cp, path, None, None, "coolwarm", false, "Test");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_visualize_cut_plane_with_colorbar() {
        use crate::core::cut_plane::CutPlane;

        let path = "test_cut_plane_cbar.png";
        let n = 16;
        let x1: Vec<f64> = (0..n).map(|i| (i % 4) as f64).collect();
        let x2: Vec<f64> = (0..n).map(|i| (i / 4) as f64).collect();
        let x3: Vec<f64> = vec![90.0; n];
        let u: Vec<f64> = (0..n).map(|i| 5.0 + i as f64 * 0.3).collect();
        let v: Vec<f64> = vec![0.0; n];
        let w: Vec<f64> = vec![0.0; n];

        let cp = CutPlane::new(
            arr1(&x1),
            arr1(&x2),
            arr1(&x3),
            arr1(&u),
            arr1(&v),
            arr1(&w),
            "z",
            (4, 4),
        );
        let result = visualize_cut_plane(
            &cp,
            path,
            Some(4.0),
            Some(10.0),
            "viridis",
            true,
            "With Colorbar",
        );
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_draw_wake_field_uses_grid_and_rotor_overlay() {
        use crate::core::cut_plane::CutPlane;

        let path = "test_wake_field.png";
        let mut x1 = Vec::new();
        let mut x2 = Vec::new();
        let mut u = Vec::new();
        for row in 0..3 {
            for col in 0..3 {
                x1.push(col as f64 * 100.0);
                x2.push(row as f64 * 100.0);
                u.push(8.0 - (row + col) as f64 * 0.25);
            }
        }
        let plane = CutPlane::new(
            arr1(&x1),
            arr1(&x2),
            arr1(&[90.0; 9]),
            arr1(&u),
            arr1(&[0.0; 9]),
            arr1(&[0.0; 9]),
            "z",
            (3, 3),
        );
        let turbines = vec![WakeTurbine {
            x: 100.0,
            y: 100.0,
            diameter: 80.0,
            label: "T1".to_string(),
        }];

        let result = draw_wake_field(&plane, &turbines, path, None, None, "Wake");
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    #[test]
    fn test_wake_cell_indices_match_planar_grid_order() {
        assert_eq!(wake_cell_indices(3, 0, 0), [0, 3, 1, 4]);
        assert_eq!(wake_cell_indices(3, 1, 1), [4, 7, 5, 8]);
    }

    #[test]
    fn test_visualize_cut_plane_empty_data() {
        use crate::core::cut_plane::CutPlane;

        let path = "test_cut_plane_empty.png";
        let cp = CutPlane::new(
            arr1(&[]),
            arr1(&[]),
            arr1(&[]),
            arr1(&[]),
            arr1(&[]),
            arr1(&[]),
            "z",
            (0, 0),
        );
        let result = visualize_cut_plane(&cp, path, None, None, "coolwarm", false, "");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("empty"));
        cleanup(path);
    }

    #[test]
    fn test_visualize_cut_plane_empty_title_uses_default() {
        use crate::core::cut_plane::CutPlane;

        let path = "test_cut_plane_default_title.png";
        let cp = CutPlane::new(
            arr1(&[0.0, 1.0]),
            arr1(&[0.0, 1.0]),
            arr1(&[90.0, 90.0]),
            arr1(&[8.0, 7.0]),
            arr1(&[0.0, 0.0]),
            arr1(&[0.0, 0.0]),
            "z",
            (2, 1),
        );
        let result = visualize_cut_plane(&cp, path, None, None, "coolwarm", false, "");
        assert!(result.is_ok());
        cleanup(path);
    }

    // ========================================================================
    // visualize_cut_plane_with_rotors tests
    // ========================================================================

    #[test]
    fn test_visualize_cut_plane_with_rotors_creates_file() {
        use crate::core::cut_plane::CutPlane;

        let path = "test_cut_plane_rotors.png";
        let cp = CutPlane::new(
            arr1(&[0.0, 1.0, 2.0]),
            arr1(&[0.0, 0.0, 0.0]),
            arr1(&[90.0, 90.0, 90.0]),
            arr1(&[8.0, 7.5, 7.0]),
            arr1(&[0.0, 0.0, 0.0]),
            arr1(&[0.0, 0.0, 0.0]),
            "z",
            (3, 1),
        );

        let fmodel = FlorisModel::from_default().unwrap();
        let result = visualize_cut_plane_with_rotors(
            &cp, &fmodel, path, None, None, "coolwarm", true, true, "Rotors",
        );
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        cleanup(path);
    }

    // ========================================================================
    // visualize_cut_plane_with_rotors_and_labels tests
    // ========================================================================

    #[test]
    fn test_visualize_cut_plane_with_labels_creates_file() {
        use crate::core::cut_plane::CutPlane;

        let path = "test_cut_plane_labels.png";
        let cp = CutPlane::new(
            arr1(&[0.0, 1.0, 2.0]),
            arr1(&[0.0, 0.0, 0.0]),
            arr1(&[90.0, 90.0, 90.0]),
            arr1(&[8.0, 7.5, 7.0]),
            arr1(&[0.0, 0.0, 0.0]),
            arr1(&[0.0, 0.0, 0.0]),
            "z",
            (3, 1),
        );

        let fmodel = FlorisModel::from_default().unwrap();
        let names = vec!["T1", "T2", "T3"];
        let result = visualize_cut_plane_with_rotors_and_labels(
            &cp, &fmodel, path, None, None, "coolwarm", &names, true, true, "Labels",
        );
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        //cleanup(path);
    }
}
