# Compare the Rust examples with their Python FLORIS counterparts.
#
# Usage (from the repository root):
#   pwsh -File tools/compare_examples.ps1
#   pwsh -File tools/compare_examples.ps1 -Append
#   pwsh -File tools/compare_examples.ps1 -Filter layout_optimization
#
# Each pair is run twice: the Rust example binary (built with
# `cargo build --release --examples`) from the repository root and the Python
# example from the folder it lives in (with MPLBACKEND=Agg so matplotlib does
# not block on `plt.show()`).
#
# The stdout of both sides is printed side by side; stderr (warnings) is shown
# separately and is not part of the comparison.
#
# The Python examples need the `floris` package of the repository: it is not
# installed in the virtual environment, so PYTHONPATH points at the `floris`
# folder of the repository.

param(
    [string]$OutFile = "output/examples_comparison.txt",
    [switch]$Append,
    # Only run the pairs whose Rust or Python example matches this pattern
    [string]$Filter = ""
)

$ErrorActionPreference = "Continue"

# The examples of the root folder of the Python examples
$root = "floris/examples"

$pairs = @(
    @{ Rust = "001_opening_floris_computing_power"; Py = "001_opening_floris_computing_power.py"; Cwd = $root }
    @{ Rust = "002_visualization"; Py = "002_visualizations.py"; Cwd = $root }
    @{ Rust = "003_wind_data_objects"; Py = "003_wind_data_objects.py"; Cwd = $root }
    @{ Rust = "004_set"; Py = "004_set.py"; Cwd = $root }
    @{ Rust = "005_getting_power"; Py = "005_getting_power.py"; Cwd = $root }
    @{ Rust = "006_get_farm_aep"; Py = "006_get_farm_aep.py"; Cwd = $root }
    @{ Rust = "007_sweeping_variables"; Py = "007_sweeping_variables.py"; Cwd = $root }
    @{ Rust = "008_uncertain_models"; Py = "008_uncertain_models.py"; Cwd = $root }
    @{ Rust = "009_parallel_models"; Py = "009_parallel_models.py"; Cwd = $root }
    @{ Rust = "010_compare_farm_power_with_neighbor"; Py = "010_compare_farm_power_with_neighbor.py"; Cwd = $root }
    @{ Rust = "011_extract_wind_speed_at_turbines"; Py = "001_extract_wind_speed_at_turbines.py"; Cwd = "floris/examples/examples_get_flow" }
    @{ Rust = "control_optimization_001_opt_yaw_single_ws"; Py = "001_opt_yaw_single_ws.py"; Cwd = "floris/examples/examples_control_optimization" }
    @{ Rust = "control_optimization_002_opt_yaw_single_ws_uncertain"; Py = "002_opt_yaw_single_ws_uncertain.py"; Cwd = "floris/examples/examples_control_optimization" }
    @{ Rust = "control_optimization_003_opt_yaw_multiple_ws"; Py = "003_opt_yaw_multiple_ws.py"; Cwd = "floris/examples/examples_control_optimization" }
    @{ Rust = "control_optimization_004_optimize_yaw_aep"; Py = "004_optimize_yaw_aep.py"; Cwd = "floris/examples/examples_control_optimization" }
    @{ Rust = "control_optimization_005_optimize_yaw_aep_parallel"; Py = "005_optimize_yaw_aep_parallel.py"; Cwd = "floris/examples/examples_control_optimization" }
    @{ Rust = "control_optimization_006_compare_yaw_optimizers"; Py = "006_compare_yaw_optimizers.py"; Cwd = "floris/examples/examples_control_optimization" }
    @{ Rust = "control_optimization_007_optimize_yaw_with_neighbor_farms"; Py = "007_optimize_yaw_with_neighbor_farms.py"; Cwd = "floris/examples/examples_control_optimization" }
    @{ Rust = "control_types_001_derating_control"; Py = "001_derating_control.py"; Cwd = "floris/examples/examples_control_types" }
    @{ Rust = "control_types_002_disable_turbines"; Py = "002_disable_turbines.py"; Cwd = "floris/examples/examples_control_types" }
    @{ Rust = "control_types_003_setting_yaw_and_disabling"; Py = "003_setting_yaw_and_disabling.py"; Cwd = "floris/examples/examples_control_types" }
    @{ Rust = "control_types_004_helix_active_wake_mixing"; Py = "004_helix_active_wake_mixing.py"; Cwd = "floris/examples/examples_control_types" }
    @{ Rust = "control_types_005_peak_shaving"; Py = "005_peak_shaving.py"; Cwd = "floris/examples/examples_control_types" }
    @{ Rust = "emgauss_001_empirical_gauss_velocity_deficit_parameters"; Py = "001_empirical_gauss_velocity_deficit_parameters.py"; Cwd = "floris/examples/examples_emgauss" }
    @{ Rust = "emgauss_002_empirical_gauss_deflection_parameters"; Py = "002_empirical_gauss_deflection_parameters.py"; Cwd = "floris/examples/examples_emgauss" }
    @{ Rust = "floating_001_floating_turbine_models"; Py = "001_floating_turbine_models.py"; Cwd = "floris/examples/examples_floating" }
    @{ Rust = "floating_002_floating_vs_fixedbottom_farm"; Py = "002_floating_vs_fixedbottom_farm.py"; Cwd = "floris/examples/examples_floating" }
    @{ Rust = "floating_003_tilt_driven_vertical_wake_deflection"; Py = "003_tilt_driven_vertical_wake_deflection.py"; Cwd = "floris/examples/examples_floating" }
    @{ Rust = "get_flow_001_get_flow"; Py = "001_extract_wind_speed_at_turbines.py"; Cwd = "floris/examples/examples_get_flow" }
    @{ Rust = "get_flow_002_get_flow_at_turbines"; Py = "002_extract_wind_speed_at_points.py"; Cwd = "floris/examples/examples_get_flow" }
    @{ Rust = "get_flow_003_get_flow_on_grid"; Py = "003_extract_turbulence_intensity_at_points.py"; Cwd = "floris/examples/examples_get_flow" }
    @{ Rust = "get_flow_004_get_flow_with_custom_planes"; Py = "004_plot_velocity_deficit_profiles.py"; Cwd = "floris/examples/examples_get_flow" }
    @{ Rust = "heterogeneous_001_heterogeneous_inflow"; Py = "001_heterogeneous_inflow_single.py"; Cwd = "floris/examples/examples_heterogeneous" }
    @{ Rust = "heterogeneous_002_het_map_from_file"; Py = "002_heterogeneous_using_wind_data.py"; Cwd = "floris/examples/examples_heterogeneous" }
    @{ Rust = "heterogeneous_003_het_speedup_by_wd_ws"; Py = "003_heterogeneous_speedup_by_wd_and_ws.py"; Cwd = "floris/examples/examples_heterogeneous" }
    @{ Rust = "heterogeneous_004_het_with_wind_rose"; Py = "004_heterogeneous_2d_and_3d.py"; Cwd = "floris/examples/examples_heterogeneous" }
    @{ Rust = "layout_optimization_001_layout_optimization"; Py = "001_optimize_layout.py"; Cwd = "floris/examples/examples_layout_optimization" }
    @{ Rust = "layout_optimization_002_layout_optimization_with_wind_rose"; Py = "002_optimize_layout_with_heterogeneity.py"; Cwd = "floris/examples/examples_layout_optimization" }
    @{ Rust = "layout_optimization_003_genetic_random_search"; Py = "003_genetic_random_search.py"; Cwd = "floris/examples/examples_layout_optimization" }
    @{ Rust = "layout_optimization_004_layout_optimization_gridded"; Py = "004_generate_gridded_layout.py"; Cwd = "floris/examples/examples_layout_optimization" }
    @{ Rust = "layout_optimization_005_layout_optimization_complex_boundary"; Py = "005_layout_optimization_complex_boundary.py"; Cwd = "floris/examples/examples_layout_optimization" }
    @{ Rust = "load_optimization_001_load_optimization"; Py = "001_lti_and_voc.py"; Cwd = "floris/examples/examples_load_optimization" }
    @{ Rust = "load_optimization_002_load_optimization_with_wind_rose"; Py = "002_row_opt_example.py"; Cwd = "floris/examples/examples_load_optimization" }
    @{ Rust = "multidim_001_multi_dimensional_cp_ct"; Py = "001_multi_dimensional_cp_ct.py"; Cwd = "floris/examples/examples_multidim" }
    @{ Rust = "multidim_002_multi_dimensional_cp_ct_2Hs"; Py = "002_multi_dimensional_cp_ct_2Hs.py"; Cwd = "floris/examples/examples_multidim" }
    @{ Rust = "multidim_003_multi_dimensional_cp_ct_TI"; Py = "003_multi_dimensional_cp_ct_TI.py"; Cwd = "floris/examples/examples_multidim" }
    @{ Rust = "operation_models_001_operation_model"; Py = "001_compare_yaw_loss.py"; Cwd = "floris/examples/examples_operation_models" }
    @{ Rust = "turbine_001_turbine_library"; Py = "001_reference_turbines.py"; Cwd = "floris/examples/examples_turbine" }
    @{ Rust = "turbine_002_turbine_interaction"; Py = "002_multiple_turbine_types.py"; Cwd = "floris/examples/examples_turbine" }
    @{ Rust = "turbine_003_custom_turbine"; Py = "003_specify_turbine_power_curve.py"; Cwd = "floris/examples/examples_turbine" }
    @{ Rust = "turbopark_001_turbopark_model"; Py = "001_compare_turbopark_implementations.py"; Cwd = "floris/examples/examples_turbopark" }
    @{ Rust = "uncertain_001_uncertain_floris_model"; Py = "001_uncertain_model_params.py"; Cwd = "floris/examples/examples_uncertain" }
    @{ Rust = "uncertain_002_approx_floris_model"; Py = "002_approx_floris_model.py"; Cwd = "floris/examples/examples_uncertain" }
    @{ Rust = "uncertain_003_parallel_uncertain"; Py = "003_uncertain_model_with_parallelization.py"; Cwd = "floris/examples/examples_uncertain" }
    @{ Rust = "visualizations_001_layout_visualizations"; Py = "001_layout_visualizations.py"; Cwd = "floris/examples/examples_visualizations" }
    @{ Rust = "visualizations_002_visualize_y_cut_plane"; Py = "002_visualize_y_cut_plane.py"; Cwd = "floris/examples/examples_visualizations" }
    @{ Rust = "visualizations_003_visualize_cross_plane"; Py = "003_visualize_cross_plane.py"; Cwd = "floris/examples/examples_visualizations" }
    @{ Rust = "visualizations_004_visualize_rotor_values"; Py = "004_visualize_rotor_values.py"; Cwd = "floris/examples/examples_visualizations" }
    @{ Rust = "visualizations_005_visualize_flow_by_sweeping_turbines"; Py = "005_visualize_flow_by_sweeping_turbines.py"; Cwd = "floris/examples/examples_visualizations" }
    @{ Rust = "wind_data_001_wind_data_comparisons"; Py = "001_wind_data_comparisons.py"; Cwd = "floris/examples/examples_wind_data" }
    @{ Rust = "wind_data_002_generate_ti"; Py = "002_generate_ti.py"; Cwd = "floris/examples/examples_wind_data" }
    @{ Rust = "wind_data_003_generate_value"; Py = "003_generate_value.py"; Cwd = "floris/examples/examples_wind_data" }
    @{ Rust = "wind_resource_grid_001_wind_resource_grid"; Py = "001_wind_rose_wrg.py"; Cwd = "floris/examples/examples_wind_resource_grid" }
    @{ Rust = "wind_resource_grid_002_wrg_from_file"; Py = "002_set_floris_with_wrg.py"; Cwd = "floris/examples/examples_wind_resource_grid" }
    @{ Rust = "wind_resource_grid_003_wrg_layout_optimization"; Py = "003_wrg_compar_layout_optimization.py"; Cwd = "floris/examples/examples_wind_resource_grid" }
)

if ($Filter -ne "") {
    $pairs = $pairs | Where-Object { $_.Rust -match $Filter -or $_.Py -match $Filter }
}

$python = "floris\.venv\Scripts\python.exe"
$repoRoot = Split-Path -Parent $PSScriptRoot

function Get-Normalized([string[]]$lines) {
    $lines | ForEach-Object { $_.TrimEnd() } | Where-Object { $_.Trim().Length -gt 0 }
}

$out = New-Object System.Collections.Generic.List[string]
$out.Add("=" * 100)
$out.Add(" FLORIS (Rust) vs FLORIS (Python) example output comparison")
$out.Add(" Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')")
$out.Add(" Rust examples:    cargo build --release --examples, run from repo root")
$out.Add(" Python examples:  $python (MPLBACKEND=Agg), run from their example folder")
$out.Add("")
$out.Add(' "STDOUT IDENTICAL" = stdout equal after trimming per-line trailing')
$out.Add(" whitespace and dropping empty lines. stderr (warnings/errors) is shown")
$out.Add(" separately and not part of the comparison.")
$out.Add("=" * 100)

$index = 0
foreach ($pair in $pairs) {
    $index++
    $exe = Join-Path $repoRoot "target\release\examples\$($pair.Rust).exe"
    $out.Add("")
    $out.Add("=" * 100)
    $out.Add("[$index/$($pairs.Count)] EXAMPLE: $($pair.Rust)  vs  $($pair.Py)")

    # The examples are run once, with the output and the error messages in
    # separate files
    $rustOutFile = [System.IO.Path]::GetTempFileName()
    $rustErrFile = [System.IO.Path]::GetTempFileName()
    if (Test-Path $exe) {
        $process = Start-Process -FilePath $exe -WorkingDirectory $repoRoot -NoNewWindow -Wait -PassThru `
            -RedirectStandardOutput $rustOutFile -RedirectStandardError $rustErrFile
        $rustOut = @(Get-Content $rustOutFile -ErrorAction SilentlyContinue)
        $rustErr = @(Get-Content $rustErrFile -ErrorAction SilentlyContinue)
        $out.Add("  Rust   exit=$($process.ExitCode)")
    } else {
        $out.Add("  (rust example not built: $exe)")
    }
    Remove-Item $rustOutFile, $rustErrFile -ErrorAction SilentlyContinue

    $pyOutFile = [System.IO.Path]::GetTempFileName()
    $pyErrFile = [System.IO.Path]::GetTempFileName()
    $env:MPLBACKEND = "Agg"
    $env:PYTHONPATH = (Join-Path $repoRoot "floris")
    $process = Start-Process -FilePath (Join-Path $repoRoot $python) -ArgumentList $pair.Py `
        -WorkingDirectory (Join-Path $repoRoot $pair.Cwd) -NoNewWindow -Wait -PassThru `
        -RedirectStandardOutput $pyOutFile -RedirectStandardError $pyErrFile
    $pyOut = @(Get-Content $pyOutFile -ErrorAction SilentlyContinue)
    $pyErr = @(Get-Content $pyErrFile -ErrorAction SilentlyContinue)
    Remove-Item $pyOutFile, $pyErrFile -ErrorAction SilentlyContinue
    $out.Add("  Python exit=$($process.ExitCode)  ($($pair.Py))")

    # The comparison of the two outputs
    $rustNormalized = Get-Normalized $rustOut
    $pyNormalized = Get-Normalized $pyOut
    $identical = ($rustNormalized -join "`n") -eq ($pyNormalized -join "`n")
    $out.Add("  STDOUT IDENTICAL: $(if ($identical) { 'yes' } else { 'no' })")

    $out.Add("--- RUST STDOUT ---")
    $out.Add(($rustOut -join "`n"))
    $out.Add("--- PYTHON STDOUT ---")
    $out.Add(($pyOut -join "`n"))
    if ($rustErr.Count -gt 0) {
        $out.Add("--- RUST STDERR (last 10 lines) ---")
        $out.Add((($rustErr | Select-Object -Last 10) -join "`n"))
    }
    if ($pyErr.Count -gt 0) {
        $out.Add("--- PYTHON STDERR (last 10 lines) ---")
        $out.Add((($pyErr | Select-Object -Last 10) -join "`n"))
    }
}

if ($Append) {
    Add-Content -Path $OutFile -Value ($out -join "`r`n")
    Write-Output "Appended $($pairs.Count) comparisons to $OutFile"
} else {
    $out -join "`r`n" | Set-Content -Path $OutFile
    Write-Output "Wrote $($pairs.Count) comparisons to $OutFile"
}
