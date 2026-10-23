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
    let cut_plane = fmodel.calculate_horizontal_plane(hub_height, 200, 100, None, None)?;
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
            t: 0.0,
            r: 59,
            g: 76,
            b: 192,
        },
        ColormapSample {
            t: 0.05,
            r: 73,
            g: 97,
            b: 210,
        },
        ColormapSample {
            t: 0.1,
            r: 89,
            g: 119,
            b: 227,
        },
        ColormapSample {
            t: 0.15,
            r: 106,
            g: 139,
            b: 239,
        },
        ColormapSample {
            t: 0.2,
            r: 123,
            g: 159,
            b: 249,
        },
        ColormapSample {
            t: 0.25,
            r: 141,
            g: 176,
            b: 254,
        },
        ColormapSample {
            t: 0.3,
            r: 158,
            g: 190,
            b: 255,
        },
        ColormapSample {
            t: 0.35,
            r: 175,
            g: 202,
            b: 252,
        },
        ColormapSample {
            t: 0.4,
            r: 192,
            g: 212,
            b: 245,
        },
        ColormapSample {
            t: 0.45,
            r: 207,
            g: 218,
            b: 234,
        },
        ColormapSample {
            t: 0.5,
            r: 221,
            g: 220,
            b: 220,
        },
        ColormapSample {
            t: 0.55,
            r: 233,
            g: 213,
            b: 203,
        },
        ColormapSample {
            t: 0.6,
            r: 242,
            g: 203,
            b: 183,
        },
        ColormapSample {
            t: 0.65,
            r: 246,
            g: 189,
            b: 162,
        },
        ColormapSample {
            t: 0.7,
            r: 247,
            g: 172,
            b: 142,
        },
        ColormapSample {
            t: 0.75,
            r: 244,
            g: 152,
            b: 122,
        },
        ColormapSample {
            t: 0.8,
            r: 238,
            g: 132,
            b: 104,
        },
        ColormapSample {
            t: 0.85,
            r: 227,
            g: 108,
            b: 85,
        },
        ColormapSample {
            t: 0.9,
            r: 214,
            g: 82,
            b: 68,
        },
        ColormapSample {
            t: 0.95,
            r: 197,
            g: 51,
            b: 52,
        },
        ColormapSample {
            t: 1.0,
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
            t: 0.0,
            r: 68,
            g: 1,
            b: 84,
        },
        ColormapSample {
            t: 0.05,
            r: 71,
            g: 19,
            b: 101,
        },
        ColormapSample {
            t: 0.1,
            r: 72,
            g: 36,
            b: 117,
        },
        ColormapSample {
            t: 0.15,
            r: 70,
            g: 52,
            b: 128,
        },
        ColormapSample {
            t: 0.2,
            r: 65,
            g: 68,
            b: 135,
        },
        ColormapSample {
            t: 0.25,
            r: 59,
            g: 82,
            b: 139,
        },
        ColormapSample {
            t: 0.3,
            r: 53,
            g: 95,
            b: 141,
        },
        ColormapSample {
            t: 0.35,
            r: 47,
            g: 108,
            b: 142,
        },
        ColormapSample {
            t: 0.4,
            r: 42,
            g: 120,
            b: 142,
        },
        ColormapSample {
            t: 0.45,
            r: 37,
            g: 132,
            b: 142,
        },
        ColormapSample {
            t: 0.5,
            r: 33,
            g: 145,
            b: 140,
        },
        ColormapSample {
            t: 0.55,
            r: 30,
            g: 156,
            b: 137,
        },
        ColormapSample {
            t: 0.6,
            r: 34,
            g: 168,
            b: 132,
        },
        ColormapSample {
            t: 0.65,
            r: 47,
            g: 180,
            b: 124,
        },
        ColormapSample {
            t: 0.7,
            r: 68,
            g: 191,
            b: 112,
        },
        ColormapSample {
            t: 0.75,
            r: 94,
            g: 201,
            b: 98,
        },
        ColormapSample {
            t: 0.8,
            r: 122,
            g: 209,
            b: 81,
        },
        ColormapSample {
            t: 0.85,
            r: 155,
            g: 217,
            b: 60,
        },
        ColormapSample {
            t: 0.9,
            r: 189,
            g: 223,
            b: 38,
        },
        ColormapSample {
            t: 0.95,
            r: 223,
            g: 227,
            b: 24,
        },
        ColormapSample {
            t: 1.0,
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
            t: 0.0,
            r: 13,
            g: 8,
            b: 135,
        },
        ColormapSample {
            t: 0.05,
            r: 42,
            g: 5,
            b: 147,
        },
        ColormapSample {
            t: 0.1,
            r: 65,
            g: 4,
            b: 157,
        },
        ColormapSample {
            t: 0.15,
            r: 86,
            g: 1,
            b: 164,
        },
        ColormapSample {
            t: 0.2,
            r: 106,
            g: 0,
            b: 168,
        },
        ColormapSample {
            t: 0.25,
            r: 126,
            g: 3,
            b: 168,
        },
        ColormapSample {
            t: 0.3,
            r: 143,
            g: 13,
            b: 164,
        },
        ColormapSample {
            t: 0.35,
            r: 161,
            g: 27,
            b: 155,
        },
        ColormapSample {
            t: 0.4,
            r: 177,
            g: 42,
            b: 144,
        },
        ColormapSample {
            t: 0.45,
            r: 191,
            g: 57,
            b: 132,
        },
        ColormapSample {
            t: 0.5,
            r: 204,
            g: 71,
            b: 120,
        },
        ColormapSample {
            t: 0.55,
            r: 214,
            g: 85,
            b: 109,
        },
        ColormapSample {
            t: 0.6,
            r: 225,
            g: 100,
            b: 98,
        },
        ColormapSample {
            t: 0.65,
            r: 234,
            g: 116,
            b: 87,
        },
        ColormapSample {
            t: 0.7,
            r: 242,
            g: 132,
            b: 75,
        },
        ColormapSample {
            t: 0.75,
            r: 248,
            g: 149,
            b: 64,
        },
        ColormapSample {
            t: 0.8,
            r: 252,
            g: 166,
            b: 54,
        },
        ColormapSample {
            t: 0.85,
            r: 254,
            g: 186,
            b: 44,
        },
        ColormapSample {
            t: 0.9,
            r: 252,
            g: 206,
            b: 37,
        },
        ColormapSample {
            t: 0.95,
            r: 247,
            g: 228,
            b: 37,
        },
        ColormapSample {
            t: 1.0,
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
            t: 0.0,
            r: 0,
            g: 0,
            b: 4,
        },
        ColormapSample {
            t: 0.05,
            r: 7,
            g: 5,
            b: 27,
        },
        ColormapSample {
            t: 0.1,
            r: 22,
            g: 11,
            b: 57,
        },
        ColormapSample {
            t: 0.15,
            r: 43,
            g: 11,
            b: 87,
        },
        ColormapSample {
            t: 0.2,
            r: 66,
            g: 10,
            b: 104,
        },
        ColormapSample {
            t: 0.25,
            r: 87,
            g: 16,
            b: 110,
        },
        ColormapSample {
            t: 0.3,
            r: 106,
            g: 23,
            b: 110,
        },
        ColormapSample {
            t: 0.35,
            r: 127,
            g: 30,
            b: 108,
        },
        ColormapSample {
            t: 0.4,
            r: 147,
            g: 38,
            b: 103,
        },
        ColormapSample {
            t: 0.45,
            r: 168,
            g: 46,
            b: 95,
        },
        ColormapSample {
            t: 0.5,
            r: 188,
            g: 55,
            b: 84,
        },
        ColormapSample {
            t: 0.55,
            r: 204,
            g: 66,
            b: 72,
        },
        ColormapSample {
            t: 0.6,
            r: 221,
            g: 81,
            b: 58,
        },
        ColormapSample {
            t: 0.65,
            r: 234,
            g: 99,
            b: 42,
        },
        ColormapSample {
            t: 0.7,
            r: 243,
            g: 120,
            b: 25,
        },
        ColormapSample {
            t: 0.75,
            r: 249,
            g: 142,
            b: 9,
        },
        ColormapSample {
            t: 0.8,
            r: 252,
            g: 165,
            b: 10,
        },
        ColormapSample {
            t: 0.85,
            r: 251,
            g: 190,
            b: 35,
        },
        ColormapSample {
            t: 0.9,
            r: 246,
            g: 215,
            b: 70,
        },
        ColormapSample {
            t: 0.95,
            r: 241,
            g: 239,
            b: 117,
        },
        ColormapSample {
            t: 1.0,
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
            t: 0.0,
            r: 0,
            g: 0,
            b: 4,
        },
        ColormapSample {
            t: 0.05,
            r: 6,
            g: 5,
            b: 26,
        },
        ColormapSample {
            t: 0.1,
            r: 20,
            g: 14,
            b: 54,
        },
        ColormapSample {
            t: 0.15,
            r: 37,
            g: 18,
            b: 85,
        },
        ColormapSample {
            t: 0.2,
            r: 59,
            g: 15,
            b: 112,
        },
        ColormapSample {
            t: 0.25,
            r: 81,
            g: 18,
            b: 124,
        },
        ColormapSample {
            t: 0.3,
            r: 100,
            g: 26,
            b: 128,
        },
        ColormapSample {
            t: 0.35,
            r: 120,
            g: 34,
            b: 129,
        },
        ColormapSample {
            t: 0.4,
            r: 140,
            g: 41,
            b: 129,
        },
        ColormapSample {
            t: 0.45,
            r: 161,
            g: 48,
            b: 126,
        },
        ColormapSample {
            t: 0.5,
            r: 183,
            g: 55,
            b: 121,
        },
        ColormapSample {
            t: 0.55,
            r: 202,
            g: 62,
            b: 114,
        },
        ColormapSample {
            t: 0.6,
            r: 222,
            g: 73,
            b: 104,
        },
        ColormapSample {
            t: 0.65,
            r: 237,
            g: 90,
            b: 95,
        },
        ColormapSample {
            t: 0.7,
            r: 247,
            g: 112,
            b: 92,
        },
        ColormapSample {
            t: 0.75,
            r: 252,
            g: 137,
            b: 97,
        },
        ColormapSample {
            t: 0.8,
            r: 254,
            g: 159,
            b: 109,
        },
        ColormapSample {
            t: 0.85,
            r: 254,
            g: 183,
            b: 126,
        },
        ColormapSample {
            t: 0.9,
            r: 254,
            g: 207,
            b: 146,
        },
        ColormapSample {
            t: 0.95,
            r: 253,
            g: 231,
            b: 169,
        },
        ColormapSample {
            t: 1.0,
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

/// The turbines that are drawn on top of a cut plane.
///
/// The rotors are drawn as in Python's `plot_turbine_rotors`: a line through
/// the center of the turbine, perpendicular to the yaw direction in the
/// inertial frame.
#[derive(Debug, Clone)]
pub struct CutPlaneTurbines {
    /// The x coordinates of the turbines.
    pub x: Vec<f64>,
    /// The y coordinates of the turbines.
    pub y: Vec<f64>,
    /// The yaw angles of the turbines in degrees (the first condition).
    pub yaw_angles: Vec<f64>,
    /// The rotor diameters of the turbines.
    pub rotor_diameters: Vec<f64>,
    /// The labels of the turbines (the turbine indices are used if None).
    pub labels: Option<Vec<String>>,
    /// The wind direction in degrees (the first condition).
    pub wind_direction: f64,
}

impl CutPlaneTurbines {
    /// The turbines of the model (the first condition).
    ///
    /// `labels` are the labels of the turbines; the turbine indices are used if
    /// they are not given.
    pub fn from_model(fmodel: &FlorisModel, labels: Option<&[&str]>) -> Result<Self> {
        let farm = &fmodel.core().farm;
        let x = farm.layout_x.to_vec();
        let y = farm.layout_y.to_vec();
        let n_turbines = x.len();

        let yaw_angles = (0..n_turbines)
            .map(|i| {
                if farm.yaw_angles.nrows() > 0 && farm.yaw_angles.ncols() > i {
                    farm.yaw_angles[[0, i]]
                } else {
                    0.0
                }
            })
            .collect::<Vec<f64>>();
        let rotor_diameters = (0..n_turbines)
            .map(|i| farm.rotor_diameters.get(i).copied().unwrap_or(0.0))
            .collect::<Vec<f64>>();

        let labels = match labels {
            Some(names) => {
                if names.len() != n_turbines {
                    anyhow::bail!(
                        "Number of turbine names ({}) does not match number of turbines ({})",
                        names.len(),
                        n_turbines
                    );
                }
                Some(names.iter().map(|name| name.to_string()).collect())
            }
            None => None,
        };

        let wind_direction = fmodel
            .flow_field()
            .wind_directions
            .first()
            .copied()
            .unwrap_or(270.0);

        Ok(Self {
            x,
            y,
            yaw_angles,
            rotor_diameters,
            labels,
            wind_direction,
        })
    }
}

/// Visualize a cut plane.
///
/// The flow field is drawn as a filled contour: the cells of the grid are
/// filled with the color of the mean velocity of the cell, as Python's
/// `tricontourf` fills the triangles of the triangulation of the points.
/// Corresponds to `visualize_cut_plane()` in Python FLORIS.
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
        None,
        false,
        false,
    )
}

fn colormap_to_rgb(cmap_name: &str, t: f64) -> (u8, u8, u8) {
    let c = get_colormap(cmap_name)(t);
    (c.0, c.1, c.2)
}

/// The scale of the rendering of the cut planes.
///
/// The figure is rendered at a higher resolution and scaled down afterwards,
/// which smooths the edges of the filled contours (matplotlib anti-aliases
/// them).
const RENDER_SCALE: u32 = 2;

/// The size of the figures of the cut planes.
const FIGURE_SIZE: (u32, u32) = (800, 600);

/// The aspect ratio (width over height) of the plot area of the cut plane
/// figures: the figure is 800 x 600 pixels, with the margins (10 pixels), the
/// caption (about 35 pixels) and the axis labels (50 and 60 pixels).
const PLOT_AREA_ASPECT: f64 = 720.0 / 495.0;

/// Extend the ranges of the axes so that the scale is the same in both
/// directions (Python's `axis("equal")`).
fn equal_aspect(
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    aspect: f64,
) -> (f64, f64, f64, f64) {
    let x_span = (x_max - x_min).max(f64::EPSILON);
    let y_span = (y_max - y_min).max(f64::EPSILON);
    if x_span / y_span > aspect {
        // The x direction is too wide: the y range is extended
        let new_y_span = x_span / aspect;
        let center = 0.5 * (y_min + y_max);
        (
            x_min,
            x_max,
            center - 0.5 * new_y_span,
            center + 0.5 * new_y_span,
        )
    } else {
        // The y direction is too wide: the x range is extended
        let new_x_span = y_span * aspect;
        let center = 0.5 * (x_min + x_max);
        (
            center - 0.5 * new_x_span,
            center + 0.5 * new_x_span,
            y_min,
            y_max,
        )
    }
}

/// The color of a band of the filled contour.
///
/// The color of a band is the color of the midpoint of the band, mapped
/// through the normalization of the velocities (Python's `tricontourf`). The
/// region below the first level and the region above the last level are the
/// "under" and the "over" regions of the contour.
fn color_of_band(levels: &[f64], vmin: f64, vmax: f64, band: usize, cmap_name: &str) -> RGBColor {
    let t = if band == 0 {
        0.0
    } else if band >= levels.len() {
        1.0
    } else {
        let midpoint = 0.5 * (levels[band - 1] + levels[band]);
        if vmax > vmin {
            ((midpoint - vmin) / (vmax - vmin)).clamp(0.0, 1.0)
        } else {
            0.0
        }
    };
    let (r, g, b) = colormap_to_rgb(cmap_name, t);
    RGBColor(r, g, b)
}

/// Clip a polygon of the flow field against a level of the contour.
///
/// The vertices of the polygon are the positions and the values of the
/// velocity; the value is interpolated linearly along the edges of the
/// polygon, as Python's `tricontourf` interpolates the contours within the
/// triangles of the points.
fn clip_polygon(polygon: &[(f64, f64, f64)], level: f64, keep_above: bool) -> Vec<(f64, f64, f64)> {
    let mut clipped = Vec::with_capacity(polygon.len() + 1);
    let n = polygon.len();
    for i in 0..n {
        let current = polygon[i];
        let next = polygon[(i + 1) % n];
        let current_inside = if keep_above {
            current.2 >= level
        } else {
            current.2 < level
        };
        let next_inside = if keep_above {
            next.2 >= level
        } else {
            next.2 < level
        };
        if current_inside {
            clipped.push(current);
        }
        if current_inside != next_inside {
            // The intersection of the edge with the level
            let denominator = next.2 - current.2;
            let t = if denominator.abs() > f64::EPSILON {
                (level - current.2) / denominator
            } else {
                0.0
            };
            clipped.push((
                current.0 + t * (next.0 - current.0),
                current.1 + t * (next.1 - current.1),
                level,
            ));
        }
    }
    clipped
}

/// The contour levels of the flow field.
///
/// Python's `tricontourf` uses matplotlib's `MaxNLocator` when no levels are
/// given: the step of the levels is the first of the "extended steps" (1, 1.5,
/// 2, 2.5, 3, 4, 5, 6, 8, 10 times a power of ten) that is not smaller than
/// the raw step (the span of the values divided by the number of the levels).
fn contour_levels(vmin: f64, vmax: f64, n_levels: usize) -> Vec<f64> {
    const EXTENDED_STEPS: [f64; 10] = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0];

    let span = (vmax - vmin).abs().max(f64::EPSILON);
    let raw_step = span / (n_levels + 1) as f64;
    let magnitude = 10f64.powf(raw_step.log10().floor());
    let step = EXTENDED_STEPS
        .iter()
        .map(|candidate| candidate * magnitude)
        .find(|candidate| *candidate >= raw_step)
        .unwrap_or(10.0 * magnitude);

    // The levels are the multiples of the step that span the values
    let start = (vmin / step).floor() * step;
    let end = (vmax / step).ceil() * step;
    let n = ((end - start) / step).round().max(0.0) as usize;
    let ticks: Vec<f64> = (0..=n).map(|k| start + k as f64 * step).collect();

    // Matplotlib only uses the levels that are inside the range of the
    // values: the values below the first level and the values above the last
    // level are the "under" and the "over" regions of the contour.
    let interior: Vec<f64> = ticks
        .iter()
        .copied()
        .filter(|&level| level > vmin && level < vmax)
        .collect();
    if interior.is_empty() {
        vec![vmin, vmax]
    } else {
        interior
    }
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
    turbines: Option<&CutPlaneTurbines>,
    show_rotors: bool,
    show_labels: bool,
) -> Result<()> {
    let x1 = &cut_plane.data.x1;
    let x2 = &cut_plane.data.x2;
    let vel = select_velocity_component(cut_plane, vel_component)?;

    if x1.is_empty() || x2.is_empty() || vel.is_empty() {
        anyhow::bail!("Cut plane data is empty");
    }

    let vmin = min_speed.unwrap_or_else(|| vel.iter().copied().fold(f64::INFINITY, f64::min));
    let vmax = max_speed.unwrap_or_else(|| vel.iter().copied().fold(f64::NEG_INFINITY, f64::max));

    let resolution = cut_plane.resolution;
    let nx1 = resolution.0;
    let nx2 = resolution.1;
    let n_points = x1.len();

    let x_min = x1.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = x1.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_min = x2.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = x2.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    // The color bar is drawn to the right of the flow field
    let x_span = (x_max - x_min).max(1.0);
    let color_bar_width = if color_bar { 0.14 * x_span } else { 0.0 };

    // The axes have an equal aspect (Python's `axis("equal")`): the ranges are
    // extended so that the scale is the same in both directions
    let (x_min, x_max, y_min, y_max) = equal_aspect(x_min, x_max, y_min, y_max, PLOT_AREA_ASPECT);
    let x_span = (x_max - x_min).max(1.0);

    let title_text = if title.is_empty() {
        "Flow Visualization".to_string()
    } else {
        title.to_string()
    };

    let scale = RENDER_SCALE;
    let root = BitMapBackend::new(
        output_path.as_ref(),
        (FIGURE_SIZE.0 * scale, FIGURE_SIZE.1 * scale),
    )
    .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(title_text, ("sans-serif", 25 * scale).into_font())
        .margin(10 * scale)
        .x_label_area_size(50 * scale)
        .y_label_area_size(60 * scale)
        .build_cartesian_2d(x_min..(x_max + color_bar_width), y_min..y_max)?;

    // The descriptions of the axes depend on the normal vector of the plane
    // (Python's cut plane holds the coordinates of the plane in x1 and x2)
    let (x1_desc, x2_desc) = match cut_plane.normal_vector.as_str() {
        "y" => ("X Position (m)", "Z Position (m)"),
        "x" => ("Y Position (m)", "Z Position (m)"),
        _ => ("X Position (m)", "Y Position (m)"),
    };

    chart
        .configure_mesh()
        .x_desc(x1_desc)
        .y_desc(x2_desc)
        .label_style(("sans-serif", 12 * scale).into_font())
        .axis_desc_style(("sans-serif", 15 * scale).into_font())
        // Matplotlib does not draw a grid by default
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    // The color of a velocity: the field is drawn as a filled contour, so the
    // velocity is quantized to the contour levels and the color of the band is
    // the color of the midpoint of the band (Python's `tricontourf`).
    let levels = contour_levels(vmin, vmax, 7);
    let color = |value: f64| -> RGBColor {
        let band = levels.partition_point(|&level| level <= value);
        let t = if band == 0 {
            0.0
        } else if band >= levels.len() {
            1.0
        } else {
            let midpoint = 0.5 * (levels[band - 1] + levels[band]);
            if vmax > vmin {
                ((midpoint - vmin) / (vmax - vmin)).clamp(0.0, 1.0)
            } else {
                0.0
            }
        };
        let (r, g, b) = colormap_to_rgb(cmap_name, t);
        RGBColor(r, g, b)
    };

    if nx1 > 1 && nx2 > 1 && nx1 * nx2 == n_points {
        // The cells of the (rotated) grid. The cells are clipped against the
        // levels of the contour, so that the boundaries of the bands are
        // interpolated (Python's `tricontourf` interpolates the contours
        // within the triangles of the points).
        //
        // The cells of a band are drawn together, so that the seams of the
        // rasterization of the polygons (the pixels on the edges of the
        // polygons are not filled) have the same color on both sides.
        let under_color = color_of_band(&levels, vmin, vmax, 0, cmap_name);
        let over_color = color_of_band(&levels, vmin, vmax, levels.len(), cmap_name);
        let band_colors: Vec<RGBColor> = (1..levels.len())
            .map(|band| color_of_band(&levels, vmin, vmax, band, cmap_name))
            .collect();

        // The cells of the grid
        let cells: Vec<[(f64, f64, f64); 4]> = (0..(nx2 - 1))
            .flat_map(|j| {
                (0..(nx1 - 1)).map(move |i| {
                    let p00 = j * nx1 + i;
                    let p10 = p00 + 1;
                    let p01 = p00 + nx1;
                    let p11 = p01 + 1;
                    [
                        (x1[p00], x2[p00], vel[p00]),
                        (x1[p10], x2[p10], vel[p10]),
                        (x1[p11], x2[p11], vel[p11]),
                        (x1[p01], x2[p01], vel[p01]),
                    ]
                })
            })
            .collect();

        // The region below the first level
        for cell in &cells {
            let under = clip_polygon(cell, levels[0], false);
            if under.len() >= 3 {
                chart.draw_series(std::iter::once(Polygon::new(
                    under.iter().map(|&(x, y, _)| (x, y)).collect::<Vec<_>>(),
                    under_color.filled(),
                )))?;
            }
        }

        // The region above the last level
        for cell in &cells {
            let over = clip_polygon(cell, levels[levels.len() - 1], true);
            if over.len() >= 3 {
                chart.draw_series(std::iter::once(Polygon::new(
                    over.iter().map(|&(x, y, _)| (x, y)).collect::<Vec<_>>(),
                    over_color.filled(),
                )))?;
            }
        }

        // The bands between the levels
        for (band, level) in levels.iter().enumerate().take(levels.len() - 1) {
            for cell in &cells {
                let above = clip_polygon(cell, *level, true);
                let band_polygon = clip_polygon(&above, levels[band + 1], false);
                if band_polygon.len() >= 3 {
                    chart.draw_series(std::iter::once(Polygon::new(
                        band_polygon
                            .iter()
                            .map(|&(x, y, _)| (x, y))
                            .collect::<Vec<_>>(),
                        band_colors[band].filled(),
                    )))?;
                }
            }
        }
    } else {
        // The points are not on a regular grid: they are drawn as markers
        for k in 0..n_points {
            chart.draw_series(std::iter::once(Circle::new(
                (x1[k], x2[k]),
                3,
                color(vel[k]).filled(),
            )))?;
        }
    }

    // The color bar
    if color_bar {
        let bar_x0 = x_max + 0.09 * x_span;
        let bar_x1 = x_max + 0.12 * x_span;
        let y_span = y_max - y_min;
        for k in 0..100 {
            let t = k as f64 / 100.0;
            let y0 = y_min + t * y_span;
            let y1 = y_min + (k + 1) as f64 / 100.0 * y_span;
            chart.draw_series(std::iter::once(Rectangle::new(
                [(bar_x0, y0), (bar_x1, y1)],
                color(vmin + t * (vmax - vmin)).filled(),
            )))?;
        }
    }

    // The turbine rotors and their labels
    if let Some(turbines) = turbines {
        // The yaw angles are rotated to the inertial frame for the plotting
        // (Python's `plot_turbine_rotors`)
        let wind_delta = 270.0 - turbines.wind_direction;
        for i in 0..turbines.x.len() {
            let radius = turbines.rotor_diameters[i] / 2.0;
            if show_rotors {
                let yaw = (turbines.yaw_angles[i] - wind_delta).to_radians();
                let x_0 = turbines.x[i] + yaw.sin() * radius;
                let x_1 = turbines.x[i] - yaw.sin() * radius;
                let y_0 = turbines.y[i] - yaw.cos() * radius;
                let y_1 = turbines.y[i] + yaw.cos() * radius;
                chart.draw_series(LineSeries::new(
                    vec![(x_0, y_0), (x_1, y_1)],
                    BLACK.stroke_width(2),
                ))?;
            }
            if show_labels {
                let label = match &turbines.labels {
                    Some(labels) => labels[i].clone(),
                    None => format!("{:03}", i),
                };
                // The labels are offset by a fraction of the rotor radius
                // (Python's `plot_turbine_labels`)
                let offset = radius / 8.0;
                chart.draw_series(std::iter::once(Text::new(
                    label,
                    (turbines.x[i] + offset, turbines.y[i] + offset),
                    ("sans-serif", 15 * scale).into_font().color(&BLACK),
                )))?;
            }
        }
    }

    root.present()?;

    // The figure is scaled down to smooth the edges of the filled contours
    if scale > 1 {
        let image = image::open(output_path.as_ref())?.to_rgba8();
        let resized = image::imageops::resize(
            &image,
            FIGURE_SIZE.0,
            FIGURE_SIZE.1,
            image::imageops::FilterType::Triangle,
        );
        resized.save(output_path.as_ref())?;
    }

    println!(
        "  Saved cut plane visualization to: {}",
        output_path.as_ref().display()
    );
    Ok(())
}

/// Visualize a cut plane with the turbine rotors.
///
/// Corresponds to `visualize_cut_plane()` followed by `plot_turbine_rotors()`
/// in Python FLORIS.
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
    let turbines = CutPlaneTurbines::from_model(fmodel, None)?;

    visualize_cut_plane_inner(
        cut_plane,
        output_path,
        "u",
        min_speed,
        max_speed,
        cmap_name,
        color_bar,
        title,
        Some(&turbines),
        show_rotors,
        false,
    )
}

/// Visualize a cut plane with the turbine rotors and their labels.
///
/// Corresponds to `visualize_cut_plane()` followed by `plot_turbine_rotors()`
/// and `plot_turbine_labels()` in Python FLORIS.
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
    let turbines = CutPlaneTurbines::from_model(fmodel, Some(turbine_names))?;

    visualize_cut_plane_inner(
        cut_plane,
        output_path,
        "u",
        min_speed,
        max_speed,
        cmap_name,
        color_bar,
        title,
        Some(&turbines),
        show_rotors,
        true,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
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
        assert_eq!(c_mid.0, 221);
        assert_eq!(c_mid.1, 220);
        assert_eq!(c_mid.2, 220);
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
        let names: Vec<String> = (0..fmodel.layout_x().len())
            .map(|i| format!("T{}", i + 1))
            .collect();
        let name_refs: Vec<&str> = names.iter().map(|name| name.as_str()).collect();
        let result = visualize_cut_plane_with_rotors_and_labels(
            &cp, &fmodel, path, None, None, "coolwarm", &name_refs, true, true, "Labels",
        );
        assert!(result.is_ok());
        assert!(std::fs::metadata(path).is_ok());
        //cleanup(path);
    }

    /// The clipping of a cell of the flow field against the levels of the
    /// contour.
    #[test]
    fn test_clip_polygon() {
        // A cell of the grid of the horizontal plane of the FLORIS example
        // 002 (the values are the velocities at its corners)
        let cell = [
            (2697.054789, 846.069230, 7.138),
            (2704.979000, 867.840000, 6.8852),
            (2733.074000, 857.614000, 6.9074),
            (2725.150490, 835.843232, 7.1502),
        ];

        // The region above the level 6 is the whole cell
        let above = clip_polygon(&cell, 6.0, true);
        assert_eq!(above.len(), 4);

        // The band between 6 and 7 is the region around the two low corners
        let band = clip_polygon(&above, 7.0, false);
        assert_eq!(band.len(), 4);
        assert_relative_eq!(band[0].2, 7.0, epsilon = 1e-9);
        assert_relative_eq!(band[1].2, 6.8852, epsilon = 1e-9);
        assert_relative_eq!(band[2].2, 6.9074, epsilon = 1e-9);
        assert_relative_eq!(band[3].2, 7.0, epsilon = 1e-9);

        // The region below the level 6 is empty
        let under = clip_polygon(&cell, 6.0, false);
        assert!(under.is_empty());

        // The region above the level 7 is the region around the two high
        // corners
        let over = clip_polygon(&cell, 7.0, true);
        assert_eq!(over.len(), 4);
    }

    /// The contour levels of the flow field follow matplotlib's
    /// `MaxNLocator`, which Python's `tricontourf` uses.
    #[test]
    fn test_contour_levels_match_matplotlib() {
        // The levels of the horizontal plane of the FLORIS example 002 (the
        // velocities span 0.4418 .. 8.0 m/s): matplotlib only uses the levels
        // that are inside the range of the values
        let levels = contour_levels(0.44175156926122305, 8.0, 7);
        let expected = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(levels.len(), expected.len());
        for (level, expected) in levels.iter().zip(expected) {
            assert_relative_eq!(*level, expected, epsilon = 1e-9);
        }
    }

    /// The ranges of the axes are extended so that the scale is the same in
    /// both directions (Python's `axis("equal")`).
    #[test]
    fn test_equal_aspect_extends_the_ranges() {
        let (x_min, x_max, y_min, y_max) = equal_aspect(0.0, 100.0, 0.0, 100.0, 2.0);
        assert_relative_eq!(x_min, -50.0, epsilon = 1e-9);
        assert_relative_eq!(x_max, 150.0, epsilon = 1e-9);
        assert_relative_eq!(y_min, 0.0, epsilon = 1e-9);
        assert_relative_eq!(y_max, 100.0, epsilon = 1e-9);

        let (x_min, x_max, y_min, y_max) = equal_aspect(0.0, 200.0, 0.0, 100.0, 2.0);
        assert_relative_eq!(x_min, 0.0, epsilon = 1e-9);
        assert_relative_eq!(x_max, 200.0, epsilon = 1e-9);
        assert_relative_eq!(y_min, 0.0, epsilon = 1e-9);
        assert_relative_eq!(y_max, 100.0, epsilon = 1e-9);
    }
}
