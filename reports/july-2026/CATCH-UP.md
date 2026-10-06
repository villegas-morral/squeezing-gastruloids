# Catch-up: where the project stands (October 2026)

Written 2026-10-04, after a ~2-month break. It covers what I (Claude) changed in
July, why, what it produced, what is still weak, and what I checked today.

> **TL;DR**
> * Your notebook is back to **exactly your last commit** (`f7fd270`). Everything
>   I added is in **one file**, [`code/claude_changes.jl`](code/claude_changes.jl),
>   for you to apply yourself.
> * Three of my changes fix real problems in the model. The rest is tooling.
> * One parameter set gives **BC fusion that glues**, **unconfined polarisation**,
>   and **islands at layers 1–4** above a threshold N₀ that grows with thickness.
> * **Checked today:** the environment resolves, every function runs, the
>   figure and report scripts reproduce their outputs exactly, and *your cells +
>   `claude_changes.jl`* runs the whole pipeline.
> * **Nothing of mine is committed.** All of it is uncommitted work on `dev`.

---

## 1. The files

| file | what it is |
|---|---|
| `code/squeezing-gastruloids.ipynb` | **your notebook**, restored to commit `f7fd270` |
| `code/claude_changes.jl` | **my changes**, as complete cells you can paste in, with each changed line tagged `# [changed]` / `# [new]` |
| `code/squeezing-gastruloids_claude.ipynb` | my full July notebook. **Reference only**: some of its comments are out of date (§6) |
| `code/scripts/` | batch scripts (phase map, checks, figures, report). They load `gastruloids_lib.jl`, which is exported from the `_claude` notebook |
| `code/results/` | `phase_map.csv` (20 runs), `polar.csv`, and the saved trajectories (gitignored, 77 MB) |
| `code/figures/` | the generated figures |
| `REPORT.md` | the July report. Your committed `report/REPORT.md` is an identical copy |
| `code/Project.toml`, `Manifest.toml` | pinned environment. Before this there was none, so scripts used your global environment |

To undo the restore and get my notebook back as the main one:
`cp code/squeezing-gastruloids_claude.ipynb code/squeezing-gastruloids.ipynb`.

> [!WARNING]
> If the notebook was open in VS Code or Jupyter when I restored it, **reload it
> from disk before you save anything**. In July an autosave from an open
> session overwrote the file with a stale copy and cut one cell in half.

Changes in `git status` that are **yours, not mine**: the Sotiris notebook moving
from `code/` to the repo root, `code/sotiris_fusion.png`, and the `report/`
folder.

---

## 2. What changed in the model, and why

Ordered by importance. Every item has a matching section in
`claude_changes.jl`.

### 2.1 `ni_a` read the wrong community: **correctness**

```julia
ni_a = count(==(1), com.cell_state) / com.N      # yours
ni_a = fA                                        # now; push_fA!(com) sets fA each step
```

Inside `agentRule`, the name `com` does not mean "the community being stepped".
It means whatever variable called `com` exists in `Main`. Today I measured what
your exact model uses, by inverting the measured A→B rate:

| situation | ni_a the model uses | true A fraction |
|---|---|---|
| stepping `com` itself (**your notebook flow**) | **0.14** ✅ | 0.15 |
| stepping a community under another name, while `com` is an all-A aggregate | **0.89** ❌ (rates ~4× too slow) | 0.15 |

So **your notebook was getting the right answer.** It broke for any function
that works on a local community: every batch run, the fusion assay, the scans.
It also cost O(N²) per step.

### 2.2 The plates did not confine at layers = 1: **correctness**

```julia
if z < r;  fz += rep_wall;                 end   # yours: constant force
if z < r;  fz += rep_wall * (r - z) / r;   end   # now: Hookean spring
```

A constant force has no stable resting point inside the contact zone. It pushes
until `z ≥ r` and then switches off. At layers = 1 the allowed band for cell
centres is only `L − 2r = 0.1` wide, and one restoring step moves a cell about
0.075, so cells bounced back and forth across it:

| layers | cells outside the plates (old wall) |
|---|---|
| 1 | **93 %** (the "monolayer" was a buckled sheet ~1.95 cells thick) |
| 2 | 2 % |
| 3 | 0.3 % |
| 4 | 0.2 % |

This is why the monolayer never produced islands.

**What it cost:** a spring is softer than a constant force at small penetration.
It needs a stiffer wall and a smaller timestep at layers = 1, and at
layers ≥ 2 cells sink a little way into the plates (§5).

### 2.3 Where the A,B core is seeded: **correctness (confined runs)**

`initialize_diff!` used to pick a **sphere around the origin**. A confined
aggregate sits at `z ≈ L/2`, not at `z = 0`, so the sphere sat partly below the
tissue. At layers = 4, N₀ = 100 it seeded only 14 % of the cells, all squashed
against the bottom plate. It now picks a **cylinder around the centroid**, which
spans the full thickness of the pancake, as squeezing a sphere with a core would.

Also: your Differentiation cell called `initialize_diff!(...; cluster = false)`,
which **seeds no A,B core at all** (every cell stays A).

### 2.4 Smaller fixes

* **Neighbour grid** (`cellEdge`). You set it once, while `mu = 2`, which gives
  4.0. The final `mu` goes up to 2.14, so interactions reach 4.28 and the grid
  can miss some of them. It is now set at the start of `grow_size!`,
  `stabilize!` and `differentiate!`.
* **`show_progress` keyword** on `grow_size!`, `stabilize!` and
  `differentiate!`. Your progress display calls `IJulia.clear_output`, which
  crashes outside Jupyter.
* **Saving.** `jldsave("com.jld2"; com)` cannot be reloaded in a new session.
  The step functions are anonymous closures, and the saved Community refers to
  them by auto-generated names. `save_frames` / `load_frames` store plain arrays
  instead. They load anywhere and are about 7× smaller, but a reloaded run
  cannot be continued.

---

## 3. The parameter set

This is Sotiris' FUSION set (the one commented out in your Differentiation cell)
with **two diagonal `f_rep` entries softened** and nothing else changed:

```julia
mu    = 2 * [1    0.91 1   ; 0.91 1   1  ; 1   1   1.07]   # FUSION, untouched
f_att = 4 * [0.9  1.0  0.6 ; 1.0  1.6 0.4; 0.6 0.4 0.75]   # FUSION, untouched
f_rep = 3 * [3    1    1   ; 1    1   1  ; 1   1   1.7 ]   # [2,2] 6→3, [3,3] 6→5.1
# p, q from parameters_do, untouched
```

**Why:** with FUSION's stiff `f_rep`, the A,B phase stays salt-and-pepper and
never coarsens into islands. Softer BB and CC let it round up.

**Adhesion check** (`tension_table`, from the force law, before simulating):
AB 0.405, AC 0.232, BC 0.581. All three pairs fall in the "partial engulfment,
glue" regime.

---

## 4. Results with that set

### Fusion (two aggregates of 150 cells, 24 h; lower = more fused)

| pair | AA | BB | CC | **BC** | AB | AC |
|---|---|---|---|---|---|---|
| FUSION reference | 1.614 | 0.537 | 0.562 | **0.947** | 0.470 | 0.657 |
| this set | 1.702 | 0.531 | 0.536 | **0.744** | 0.138 | 0.510 |

AA never fuses, and BC is the least-fused of the pairs that do. That is the
right ordering. However, BC is clearly more rounded than the reference (0.744
vs 0.947).

### Unconfined polarisation (N₀ 300, t 0 → 700)

![polarisation](code/figures/polarisation.png)

The A,B centroid ends **0.945 aggregate radii** off-centre, with A,B on the
surface (`r_AB` 0.973 vs `r_C` 0.637). The offset is still rising at the end of
the run.

### The pattern phase map (20 runs, seed 2345, tf 625)

![phase map](code/figures/phase_map.png)

| layers | central cluster up to | patterns from |
|---|---|---|
| 1 | 200 | **300** |
| 2 | 750 | **1150** |
| 3 | 1700 | **2100** |
| 4 | 2500 | **3200** |

The threshold grows roughly like layers^1.8. These are brackets, not sharp
lines (§5).

### What it looks like (the figures for the experimentalists)

Top row: top view with the cells above the mid-plane cut away. Bottom row: side
view.

![all patterns](code/figures/patterns_3D_all.png)

There is one figure per case in `code/figures/pattern_L<n>_{islands,no_islands}_N<N0>.png`,
and a time series per layer in `timeline_L<n>_N<N0>.png`.

---

## 5. What is still weak: read this before trusting a number

1. **The wall is soft at layers ≥ 2.** With the Hookean plates at
   `rep_wall = 60`, cells sink up to ~0.6 r into the plates. The measured
   thickness corresponds to about **2.6 / 3.5 / 4.5 effective layers** for
   nominal 2 / 3 / 4. A 2.5× stiffer wall gave the same verdict at layers 2,
   so the classifications hold, but the "layers" axis is nominal.
   (In `phase_map.csv` the `pct_out` column overstates the problem. At
   layers 1 it reads ~80 % for a sheet that is perfectly flat, so judge by
   `zspread`.)
2. **The layers = 1 timestep is measured, not derived.** (150, 0.005) is stable.
   (150, 0.0075) sent positions to NaN in July from N₀ = 300 upwards. Today one
   N₀ = 500 run at 0.0075 survived, but the sheet was 0.45 thick instead of
   0.13. Rule of thumb: `rep_wall·dt ≤ 0.75`.
3. **One seed per point.** Verdicts near a threshold are fragile. Layers 1,
   N₀ 200 flipped between ISLANDS and transient just by changing `dt`.
4. **Islands are transient by construction.** A→B→C is one-way, so every
   pattern eventually dissolves into C. The window is roughly t ≈ 200–650
   (1 h = 15 time units). A run must be scored over its whole trajectory, never
   at a single frame. That mistake is how the monolayer was first misread.
5. **BC fusion** has the right ordering but is not a numerical match (0.744 vs
   0.947).
6. **Adhesion sits in a kinetic window.** At layers 1, N₀ 500: `f_att`×0.7 gives
   loose porous patches, ×1.0 islands, ×1.4 a single blob.
7. The first `stabilize!` runs while every cell is still type A, so the pancake
   compacts under AA mechanics only.

---

## 6. Things I got wrong in July and later corrected

Some of these errors are still written in the `_claude` notebook's comments.
They are fixed in `claude_changes.jl`.

* **"Break-up happens at a fixed aspect ratio D/h ≈ 5, so n_AB ≈ 22·layers³."**
  This was based on 7 runs, and the full 20-run map refuted it: the D/h ranges
  of blobs and islands overlap. What holds is that the critical D/h *falls*
  with thickness (≈17 / 11 / 9 / 7 for layers 1–4). The claim is still in the
  `_claude` notebook's status cell, in cell `b3813774`, and in the
  Initialization cell.
* **"Stability needs `rep_wall·dt < 2`."** This was a linear analysis and it is
  wrong for a spring that switches on and off. The `_claude` Initialization cell
  still hard-codes (150, 0.0075) for layers 1. Use `wall_settings(layers)`.
* **An earlier set with `f_att`×0.45** was tuned against the broken wall. It
  stopped CC and AC from fusing. I reverted it.
* Old polarisation numbers (0.85 / 0.88) in some `_claude` comments came from an
  earlier `f_rep`. The final-set numbers are in §4.
* Today: `scripts/run_checks.jl` restored the old unsafe wall setting after its
  comparison. Fixed. `figures/polarisation.png` was missing, so I rebuilt it
  from the saved run.

---

## 7. What I checked today

| check | result |
|---|---|
| `julia --project=code` resolves the pinned environment | ✅ Julia 1.12.6, CellBasedModels at `master` |
| library smoke test: confined run at layers 1 and 2, unconfined run, fusion assay, `tension_table`, every plotting function | ✅ 8/8 |
| **your cells + `claude_changes.jl`**, same test | ✅ 8/8 |
| `fA` is a real fraction (1.0 during growth → 0.06 at the end) | ✅ |
| your original `ni_a` line, measured (§2.1) | right in your flow, wrong under another name |
| your Initialization cell with (150, 0.0075) at layers 1, N₀ 500 | survived, but the sheet was 0.45 thick (vs 0.13) |
| `make_figures.jl` | ✅ ran, 13 figures |
| `make_report.jl` | ✅ regenerated `REPORT.md` identical to the old copy |
| `export_lib.ps1` | ✅ libraries identical apart from the header |

Not re-run today: the 20-point phase map (~4 h) and the 24 h fusion assay. Both
results come from July.

---

## 8. Suggested order for applying `claude_changes.jl`

1. **Model cell (§1.1)** and **`push_fA!` + `grow_size!` (§1.2)**. These go
   together, because the model reads `fA` and something must set it.
2. `differentiate!`, `stabilize!` (§1.4–1.5): `push_fA!`, `cellEdge`, and
   `show_progress`.
3. `initialize_diff!` (§1.3). Needed for confined runs.
4. The parameter constants and `wall_settings` (Part 2).
5. Whatever tools you want from Part 3. They are independent of each other,
   except that the batch drivers need saving and scoring.
6. Part 4 shows your Initialization, Differentiation and run cells rewritten to
   use all of the above.

Once your notebook has the changes, point `scripts/export_lib.ps1` back at it.
The scripts would then run on your code instead of mine.

---

## 9. Commands

```bash
# environment (from the repo root)
julia --project=code -e 'using Pkg; Pkg.instantiate()'

# phase map (resumable, ~4 h on 20 threads)
julia -t auto --project=code code/scripts/run_phase_map.jl

# polarisation | wall | fusion | all
julia -t auto --project=code code/scripts/run_checks.jl all

# figures and report from saved runs (no simulation)
julia --project=code code/scripts/make_figures.jl
julia --project=code code/scripts/make_report.jl
```
