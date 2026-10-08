# Re-export scripts/gastruloids_lib.jl from the notebook, so the library and the
# notebook cannot drift apart.
#
#   cd code-claude; powershell -File scripts/export_lib.ps1
#
# The cell ids below are the function/model/parameter cells. If you add a new
# function cell to the notebook, add its id here.

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
# Read UTF-8 explicitly: PowerShell 5.1 falls back to the ANSI codepage on a
# BOM-less file and mangles every non-ASCII character on the way through.
$nb = [System.IO.File]::ReadAllText(
        (Join-Path $root "squeezing-gastruloids-claude.ipynb"),
        [System.Text.Encoding]::UTF8) | ConvertFrom-Json

$ids = @(
  '51154174',  # color_map
  'c17fa1b1',  # initialize_growth
  'f2f50d7d',  # grow_size!
  'f818a62d',  # initialize_diff!
  '9154a7b7',  # fraction_B, differentiate!
  'f06d6c5e',  # stabilize!
  '462afe86',  # plot_pancake
  'ba36859e',  # record_pancake
  '4cd7bbef',  # save_frames / load_frames / plot_frames / plot_pancake_frames
  '6b8e0c35',  # plot_gallery (3D)
  'a7f3c2d1',  # final figures: save_panel, plot_panels, plot_phase_space, record_pancake_frames, replot_figures
  '9557e7da',  # cluster_sizes / island_stats / classify_islands
  '756cfff2',  # island_timeline / polar_stats / confinement_report
  '7c1d9e42',  # thickness_timeline
  '2d4f6a81',  # B/C pattern readout: b_clusters, pattern_stats, classify_pattern
  '047d722c',  # adhesion energetics
  'cad66be3',  # fusion assay
  '8a8800de',  # PSETS, wall_settings, run_confined, phase_scan
  '5e7c1a90',  # run_unconfined
  'ff3415f8',  # plot_phase / critical_N0
  '99c84e68',  # THE MODEL
  '6ff203b2',  # parameters_tfm
  '93111bf3'   # parameters_do
)

$header = @'
# =============================================================================
# gastruloids_lib.jl -- model, parameters and analysis for the squeezing project
#
# EXPORTED FROM squeezing-gastruloids-claude.ipynb. Do not edit by hand; edit the
# notebook and re-export with scripts/export_lib.ps1 so the two cannot drift.
#
# Usage:
#   using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
#   using CellBasedModels; import CellBasedModels: update!
#   using JLD2, CairoMakie
#   include("gastruloids_lib.jl")
#   model = model_new
#
# `differentiate!` calls IJulia.clear_output only when show_progress = true, so
# batch scripts must pass show_progress = false (run_confined already does).
# =============================================================================

'@

function Export-Cells($cellIds, $hdr, $dst) {
  $out = foreach ($id in $cellIds) {
    $c = $nb.cells | Where-Object { $_.id -eq $id }
    if (-not $c) { Write-Warning "cell $id not found in notebook"; continue }
    "# ---- notebook cell $id " + ("-" * (60 - $id.Length))
    ($c.source -join '')
    ""
  }
  $noBom = New-Object System.Text.UTF8Encoding($false)
  [System.IO.File]::WriteAllText($dst, ($hdr + ($out -join "`n")), $noBom)
  Write-Output "wrote $dst"
}

Export-Cells $ids $header (Join-Path $PSScriptRoot "gastruloids_lib.jl")

# Analysis-only subset: no ABM(), so it loads without CellBasedModels and
# without paying the model-compilation cost. Enough for plotting and scoring
# saved trajectories.
$analysisIds = @(
  '51154174',  # color_map
  '4cd7bbef',  # save/load/plot_frames, plot_pancake_frames
  '6b8e0c35',  # plot_gallery (3D)
  'a7f3c2d1',  # final figures: save_panel, plot_panels, plot_phase_space, record_pancake_frames, replot_figures
  '9557e7da',  # cluster_sizes / island_stats / classify_islands
  '756cfff2',  # island_timeline / polar_stats / confinement_report
  '7c1d9e42',  # thickness_timeline
  '2d4f6a81',  # B/C pattern readout: b_clusters, pattern_stats, classify_pattern
  '047d722c',  # adhesion energetics
  'ff3415f8'   # plot_phase / critical_N0
)
$analysisHeader = @'
# =============================================================================
# gastruloids_analysis.jl -- scoring and plotting only, NO model.
#
# EXPORTED FROM squeezing-gastruloids-claude.ipynb by scripts/export_lib.ps1.
# Do not edit by hand.
#
# Loads without CellBasedModels, so it is fast and works on any machine that can
# read the saved trajectories:
#   using Statistics, Printf, JLD2, CairoMakie
#   include("gastruloids_analysis.jl")
#
# (`save_frames` is the one function here that needs a live Community; the rest
# operate on the array format from `load_frames`.)
# =============================================================================

'@
Export-Cells $analysisIds $analysisHeader (Join-Path $PSScriptRoot "gastruloids_analysis.jl")
