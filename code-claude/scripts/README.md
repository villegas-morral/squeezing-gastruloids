# Scripts: what each one is for

Run everything from the repo root with `julia --project=. code-claude/scripts/<script>`. The model and
analysis code lives in the notebook `../squeezing-gastruloids-claude.ipynb`; the two libraries below are
exported from it and are never edited by hand.

Tags: **core** = needed to reproduce the final results; **sweep** = parameter screening and tuning, only
needed to redo that history; **july** = left over from July, not used by anything current.

| file | tag | what it does |
|---|---|---|
| `gastruloids_lib.jl` | core | model, parameters, runs, analysis, plotting (exported from the notebook) |
| `gastruloids_analysis.jl` | core | the analysis and plotting part only, loads without CellBasedModels (exported) |
| `export_lib.ps1` | core | notebook → the two libraries above |
| `set_cell.ps1`, `repair_notebook.ps1` | core | replace one notebook cell from a file (with a backup); restore from the backup |
| `smoke_test.jl` | core | quick check that every function runs |
| `stage5_run.jl` | core | batch runner: a job list (`kind,cand,core,layers,N0,seed`) → trajectories + one CSV row per run |
| `stage5_pool.sh` | core | keeps the machine busy with `stage5_run.jl` queue workers |
| `stage5_common.jl` | core | shared by the analysis scripts: read the rows, re-score runs with the current readouts (cached), pancake diameters |
| `stage5_phase.jl` | core | the phase-space job lists, and its summary table |
| `stage6_figures.jl` | core | picks the runs and writes `figures/data/`, then redraws every final figure from it, and the movies |
| `stage7_presentation.jl` | core | the overview figure (`figures/presentation.pdf`) |
| `stage4_candidates.jl` | core | the named parameter sets (`stage5_run.jl` needs it) |
| `msd_monolayer.jl` | core | November: monolayer runs with one frame per 4 min, MSD, diffusion coefficient, Damköhler number |
| `grid_compare.jl` | study | Nov 3: the f_att BB × τ_div grid, six phase spaces (claude-answers/09 §2) |
| `polarisation_time.jl` | study | Nov 3: P(t) of unconfined runs by size (why large aggregates do not polarise; claude-answers/10) |
| `pol_screen.jl` | study | Nov 3: the polarisation screens (pol_nov3, pol2_nov3): time to P 0.3, P and pole-vs-dumbbell at fB 0.12, fusion (claude-answers/10 §2.4) |
| `pol_large.jl` | study | Nov 4: the large unconfined runs of the candidate sets and λ × 0.5: time to P 0.3, P, "end", and the readout frames as slices (claude-answers/10 §2.7) |
| `pol_shapes.jl` | study | Nov 3: mid-plane cross-sections of the screen's unconfined runs at fB 0.12 (pole vs dumbbell; claude-answers/10 §2.4) |
| `wall_compare.jl` | study | Nov 4: plate friction 3 and 5 on the cohesive sets (wall_nov4) against friction 1.5, per point, + readout frames (claude-answers/10 §2.8) |
| `grid_compare_b.jl` | study | Nov 3 night: do the sets that polarise large aggregates keep the confined patterns? six coarse phase spaces, seed 1 (claude-answers/10 §2.5) |
| `phase_compare.jl` | study | Nov 4: full phase spaces (3 seeds) of the candidate sets (phase_F, G, H) next to now's, on the same points (claude-answers/10 §2.6) |
| `size_table.jl` | study | Nov 3: aggregate diameter against N0, sphere and pancakes (claude-answers/10 §3) |

### `archive/` (moved there on Nov 3)

Not used by anything current; kept to redo that history. Every path in them was rewritten to point one
level up (`@__DIR__` → `joinpath(@__DIR__, "..")`), so they still run from `archive/`.

| file | tag | what it does |
|---|---|---|
| `tdiv_report.jl`, `tdiv_phase.jl` | sweep | November: 13 h division for B and C against the final runs; its phase-space sketch (answer 07) |
| `phaseE_compare.jl` | sweep | November: the f_att BB 5.8 phase space against the final one (answer 07 §8) |
| `stage2_ic.jl`, `stage2_walls.jl`, `stage2_wall_figs.jl` | sweep | October stage 2: initial condition and plate models |
| `stage4_screen.jl`, `stage4_gallery.jl`, `stage4_minimap.jl` | sweep | October stage 4: parameter screen and its galleries |
| `stage5_tune_report.jl`, `stage5_answer_figs.jl` | sweep | October stage 5a: tuning summary and the figures of claude-answers/03 |
| `make_figures.jl`, `run_phase_map.jl`, `run_checks.jl`, `probe_stability.jl` | july | the July pipeline (A/B/C model, island verdicts); older copies in `../../../archive/code-july-2026/scripts/` |

## The final cleanup (when the results are definitive)

A possible layout:

```
code-claude/
  squeezing-gastruloids-claude.ipynb
  lib/        gastruloids_lib.jl, gastruloids_analysis.jl, export_lib.ps1, set_cell.ps1
  run.jl      the batch runner (stage5_run + pool)
  phase_space.jl   job lists and tables (stage5_phase + stage5_common)
  figures.jl  every final figure and the overview (stage6 + stage7)
  diffusion.jl     (msd_monolayer)
  archive/    the sweep scripts, with a README saying which report section each belongs to
```

The July scripts can then go: they are in `../july-2026/scripts/` and in git history. The July-era functions in
the library (`island_stats`, `island_timeline`, `island_verdict`, `phase_scan`, `plot_phase`) could go with them,
once nothing calls them.
