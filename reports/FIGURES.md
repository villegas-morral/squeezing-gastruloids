# How every final figure was made

*The final figures live in `code-claude/figures/`, and the data behind them in `code-claude/figures/data/`.
Every figure can be redrawn from that data alone, without simulating anything:*

```julia
# in the notebook (last cell), or from the repo root:
#   julia --project=. -e 'include("code-claude/scripts/gastruloids_analysis.jl"); replot_figures("code-claude/figures")'
figs = replot_figures("code-claude/figures")
```

*The physics is in [MODEL.md](MODEL.md); the decisions and numbers behind the parameter choice are in
[REPORT.md](REPORT.md).*

---

## The data format

For each figure, `figures/data/<figure>/` holds:

| file | contents |
|---|---|
| `title.txt` | the figure title |
| `panels.tsv` | one row per panel: `panel, row, col, views, label, file` (tab-separated; `\|` in a label is a line break) |
| `<panel>.csv` | **one row per cell**: `x, y, z, r, cell_state` (2 = B, 3 = C) of the frame shown |
| `<panel>.jld2` | the same arrays plus `meta`: the run file, frame, time, B fraction, number of layers, N0, seed, parameter set |

The phase space has `figures/data/phase_space/phase_space.csv` (one row per run) and its three insets
(`inset_*.csv/.jld2`). In its `panels.tsv` the `row` and `col` columns hold the run's layers and N0.

The movies are made from whole trajectories, which are kept in `figures/data/movies/*.jld2` (the
`save_frames` format: per frame `x, y, z, r, cell_state, t, N`).

Views:
- **top**: seen from above, with the cells above the mid-plane left out (a cutaway), because the C roof would
  hide the B pattern
- **side**: seen from the side, all cells (shows that the slab keeps its thickness)
- **front**: the side view turned 90°

In polarised panels the aggregate is rotated so that B sits on the right (+x). Only the display is affected.

## The procedure every simulation follows

Shared by every run (`run_confined` / `run_unconfined` in the notebook; batch: `scripts/stage5_run.jl`):

1. **Growth.** One C cell grows by division to N0 cells, with the growth parameters (`parameters_tfm`). When
   confined, it grows between rigid plates set `layers` cell layers apart (gap $L = 2.1\,r \times$ layers).
2. **Relaxation**, 100 time units, with the B/C mechanics of the chosen parameter set (`parameters_do` plus
   `PSETS[set]`), still all C.
3. **Seeding** (`initialize_diff!`). 60 % of the cells become B:
   - **sharp core**: the innermost 60 % (a cylinder through the slab when confined, a sphere when free)
   - **diluted core** (`f_in = 0.75`): 75 % of the B cells inside that core, 25 % outside, densest near it

   The time is reset to 0 here.
4. A short relaxation of 25 time units, then **differentiation**:
   - B → C at the constant rate `q`, with division on
   - run until the B fraction drops to 0.10, or `tf` = 900 time units
   - one frame saved per hour (15 time units)
5. **Readout**: the first saved frame with B fraction ≤ 0.2 (≈ 96 h in the experiments), classified by
   `classify_pattern` (MODEL.md §7).

**Settings:**
- Rigid plates and `dt` = 0.015.
- **Parameter set `PSETS["final"]`** (Oct 8: `mu` BB 1.9) with plate friction `lambda_wall` = 1.5
  (`"july_frB4.5_frC6_muB1.9@lw1.5"` in the scripts), unless a figure says otherwise.
- The figures of the Oct 7 set (`PSETS["final_oct7"]`, `mu` BB 2) are kept in `code-claude/figures-oct7/`. They
  were made the same way, with `"july_frB4.5_frC6@lw1.5" sharp phase_B`.
- **The random seed fixes growth, seeding and differentiation.** So two parameter sets run with the same
  seed start from the same grown aggregate and the same B core.

---

## The figures

All of them are written by `scripts/stage6_figures.jl` (arguments: the set, the core, the phase-space folder):

```bash
julia --project=. code-claude/scripts/stage6_figures.jl "july_frB4.5_frC6_muB1.9@lw1.5" sharp phase_D all
julia --project=. code-claude/scripts/stage7_presentation.jl      # the presentation PDF, from figures/data
```

- `data` picks the runs and writes `figures/data/`.
- `plot` redraws the PNGs from `figures/data/` alone.
- `movies` makes the mp4s from the trajectories kept in `figures/data/movies/`.

Nothing is simulated at this step. The runs come from `code-claude/results/stage5/` (gitignored):
- `phase_D/`: the phase space of the final set, its unconfined runs (the "free" row), and the screen runs of
  the set copied in from `thick/`
- `tune/`: the fusion assay (`fusion_july_frB4.5_frC6_muB1.9.jld2`, copied from `thick/`) and the
  diluted-core runs of the Oct 7 set
- `phase_B/`: the phase space of the Oct 7 set (for `figures-oct7/`)

### fig1_initial_condition

**What:** the seeded aggregate (B = 60 % of the cells), sharp vs diluted core, squeezed and free.

| panel | run | frame |
|---|---|---|
| `ic_sharp_confined` | 2 layers, N0 2400, sharp core, seed 1 | the seeding frame (t = 0, fB = 0.6) |
| `ic_sharp_free` | unconfined, N0 600, sharp core, seed 1 | the seeding frame |
| `ic_fin0.75_confined` | 2 layers, N0 2400, diluted core (`f_in = 0.75`), seed 1 | the seeding frame |

- **No diluted unconfined panel:** the diluted core was never run unconfined with these mechanics.
- **Views:** confined panels show the top (mid-plane cutaway) and the side; free panels show the top cut and
  the front.
- **Source runs:**
  - the sharp ones are runs of the final set
  - the diluted ones are friction-free runs of the Oct 7 set (`july_frB4.5_frC6`) from tuning, since the
    diluted core was not run with the final set

  The seeding depends neither on the friction nor on the B–B range (the 100-unit relaxation before it is all
  C), so the same seed gives the same aggregate and B cells up to small differences in that relaxation.
  Each panel's `.jld2` names its run.

### fig2_behaviours

**What:** your panel B, the same aggregate size (N0 600) under no, medium and strong confinement, through
time. One run per row:

| row | run | outcome |
|---|---|---|
| no confinement | unconfined, N0 600, seed 1 (`phase_D/com_free_N600_july_frB4.5_frC6_muB1.9_w0_s1.jld2`) | polarises |
| medium confinement | 4 layers, N0 600, seed 1 | central cluster that shrinks |
| strong confinement | 1 layer, N0 600, seed 1 | islands |

- **How the confined rows are picked:** they are the phase-space examples of their class, the first point
  of a preference list whose majority over seeds is that class (its lowest seed with the class):
  - central cluster: (4, 600), (3, 600), (2, 150), …
  - pattern: (1, 600), (1, 1200), (2, 2400), …
  - polarised (used for the phase-space inset and the movie): (8, 300), (6, 150), …
  - pattern in a thick slab (inset and movie only, none if no point qualifies): (3, 2400), (4, 2400),
    (3, 4800), …
- **Columns:** the frames where fB first drops to 0.6 (seeding), 0.45, 0.3, 0.2 (**readout**) and 0.12.
  The label gives the time in hours (15 time units = 1 h) and fB.
- **Views:** top (mid-plane cutaway) and side; the unconfined row shows the top (mid-plane cut) and front.
- **Orientation:** the unconfined row is rotated in 3D so that the final B side points to +x. This is for
  display only.
- **Source runs:** each panel's `.jld2` holds its source run in `meta`.

### fig3_patterns

**What:** the B pattern at the readout (fB = 0.2), seed 1, at points whose majority is PATTERN.
- **Which points:** every thickness that has patterns, with its largest N0, then its smallest, then the rest,
  up to 8 points. Thinnest slabs first.
- **View:** top cutaway only.
- **Labels:** the class, the number of islands and the looseness S.

### fig4_unconfined

**What:** the polarisation check, and the "free" row of the phase space. It shows unconfined aggregates of
N0 150, 300, 600, 1200 and 2400 (seed 1, final mechanics, no plates so no friction) at the readout
fB = 0.2:
- **Labels:** the class and the B offset.
- **Views:** top (cut at the mid-plane) and front (the whole aggregate).
- **Orientation:** rotated so that B points to +x.
- **Result:** B polarises up to N0 600; at 1200 and 2400 it is still central when B runs out (both sets).

### fig5_fusion

**What:** the fusion assay (`fusion_assay` in the notebook) with the final mechanics.
- **Setup:** two aggregates of 150 cells, each grown and relaxed alone as a pure type, are placed just
  touching.
- **Run:** 24 h with no division and no differentiation.
- **Panels:** BB, CC and BC at 24 h.
- **Metric:** centre separation / (Rg1 + Rg2), about 1.6 at contact and 0.5 for one sphere.
- **Views:** top and front.

### phase_space

**What:** the phase space of the final set (`results/stage5/phase_D`).
- **Grid:** layers 1, 2, 3, 4, 6, 8, 12 × N0 75, 150, 300, 600, 1200, 2400, plus N0 4800 at 1–4 layers.
- **Seeds:** 2 per point, 1 at N0 4800, and a third seed at the points next to a threshold
  (`phase_D_refine4.csv`).
- **Refinement:** points between neighbours of different majority class, 2 seeds each:
  - N0 89, 106, 178, 212, 357, 424, 714, 849, 1009, 1427, 1697, and layers 5 and 7
  - two rounds (`scripts/stage5_phase.jl table` writes them; job lists `phase_D_refine1.csv` and
    `phase_D_refine3.csv`)
  - plus 4 layers N0 3394 and unconfined N0 849 and 4800 (`phase_D_refine1/2.csv`)
- **Seeds 1 at the screen points** are the screen runs themselves (copied in from `thick/`).

The Oct 7 set (`figures-oct7/`) used 3 seeds, 5 at the smallest points, and three refinement rounds.
- **Pies:** each point is a pie of its seeds' outcomes at the readout fB = 0.2:
  - central cluster
  - pattern (ring or islands)
  - polarised
  - other
- **Background:** the majority outcome of the nearest point.
- **Lines:** where patterns start (dashed) and up to where the aggregate polarises (dotted).
- **Insets:** the 3D panels on the right are the readout frames of the class examples (top + side), all
  seed 1:
  - central (4 layers, N0 600)
  - pattern (1 layer, N0 600)
  - pattern in a thick slab (3 layers, N0 2400; final set only)
  - polarised (8 layers, N0 300)
- **Data:** `phase_space.csv`, one row per run:
  - `layers`, `N0`, `seed`
  - the classes at fB 0.4 / 0.3 / 0.2 / 0.15 / 0.12
  - the readout statistics (`n_isl`, `f_isl`, `f_iso`, `f_max`, `S`, `spread`, `fill`, `offset`,
    `offset_max`)
  - the final thickness `layers_eff` and the cell count `N_end`
  - the run file

**The "free" row** at the top of the chart is the same mechanics without plates:
- unconfined runs at the grid N0 values from 150 up, 2 seeds each (3 at N0 849 and 1200, 1 at 4800; `free.csv`)
- what a gap much wider than the aggregate tends to

**Four versions** (Oct 8). `replot_figures` draws every phase-space folder four times, with the same runs,
pies, background and thresholds:

| file | x | y |
|---|---|---|
| `phase_space.png` | log₁₀ N0 (as in your sketch) | plate gap, in cell layers (log scale) |
| `phase_space_N0.png` | N0 | plate gap |
| `phase_space_diameter.png` | log₁₀ N0 | pancake diameter at squeezing (log scale) |
| `phase_space_diameter_N0.png` | N0 | pancake diameter at squeezing |

- **The pancake diameter** is measured on each run at the seeding frame, after growth and the 100-unit
  relaxation, which is the moment of squeezing in the model. It is 2√(2⟨ρ²⟩) + 2r, where ρ is the in-plane
  distance of the cell centres from the centre (exact for a uniform disc), divided by the cell diameter 2r.
  It is stored per run as column `D` of `phase_space.csv` (`with_diameter!` in `scripts/stage5_common.jl`,
  cached in `diameter_cache.csv` of the results folder).
- **Pie positions.** A point sits at the mean diameter of its seeds. Points of different thickness can then
  sit close together, so a pie is smaller where its neighbours are close.
- **The unconfined row** stays on top in every version, above a white break.

### phase_space_late

**What:** the same runs and chart, read later, at fB = 0.12, where the class column is `class_012`
(`readout.txt`). The insets are the same examples at fB = 0.12. It shows which outcomes are still on their
way at the readout:
- rings that break into islands
- big aggregates that polarise late

### presentation.pdf

**What:** your sketch (`notes/3d-simulations presentation idea.jpg`) as an A4 PDF to finish in Inkscape. It is
made by `scripts/stage7_presentation.jl` from `figures/data` alone (no runs are picked there):

```bash
julia --project=. code-claude/scripts/stage7_presentation.jl
```

- **A, the model ingredients:**
  - adhesion, with the range μ
  - proliferation, with τ_div
  - differentiation T+ → T− at rate q
  - friction

  All drawn as vector shapes.
- **B, the same aggregate under no, medium and strong confinement:**
  - Runs: the three rows of `fig2_behaviours`.
  - Columns: frames at fB 0.6 (seeding), 0.3 and 0.2 (readout). The time axis gives their times in hours,
    taken from the medium row; the B fraction falls at the same rate in every run.
  - Views: the unconfined row is seen from above, the confined rows from a three-quarter view. All are cut
    at the mid-plane so the B cells show.
- **C, the phase space:**
  - Data: `phase_space/phase_space.csv` and `free.csv`.
  - Axes: log₁₀ N0 against the plate gap, with the unconfined row on top.
  - Colours: the sketch's regions (grey polarisation, cyan no polarisation, pink patterns), with the region
    names at their middle.
  - Pies: smaller (`pie_scale = 0.55`).
  - Insets: the readout frames of the three rows of B.
- **Colours:** T+ (B) green, T− (C) orange, as in the sketch.
- **Editing:** each 3D panel is one embedded image (`rasterize`), and the rest is vector with real text, so
  everything can be moved and edited in Inkscape. `presentation.png` is a preview.

### movies

`figures/movies/{central,pattern,pattern_thick,polarised,unconfined}.mp4`:
- **Runs:** central (4 layers, N0 600), pattern (1 layer, N0 600), pattern in a thick slab (3 layers,
  N0 2400), polarised (8 layers, N0 300) and unconfined (N0 600), all seed 1.
- **Frames:** one per hour, from seeding to fB = 0.10.
- **Views:** top (mid-plane cutaway) and side.
- **Made by** `record_pancake_frames`, the array version of your `record_pancake`.
- **Data:** the trajectories, from seeding on, are in `figures/data/movies/*.jld2` (`load_frames` reads them).
