# Mathematical and computational modelling of embryonic organoids

![Language: Julia](https://img.shields.io/badge/language-Julia-ed207b)
![Jupyter Notebook](https://img.shields.io/badge/Notebook-Jupyter-f37821)

Agent-based model of gastruloid squeezing, built on
[CellBasedModels.jl](https://github.com/dsb-lab/CellBasedModels.jl).

Two cell types: B (T/Bra+), which starts as a core holding 60 % of the cells, and C
(T/Bra−). B differentiates into C. We want one parameter set that gives:

* **patterns**: squeezed between plates, B breaks into sparse islands in a C sea
  once the aggregate is wide enough. Thicker aggregates need more cells.
* **a central cluster**: below that threshold, B stays in the centre and
  differentiates away.
* **polarisation**: with a large enough gap the aggregate behaves as if
  unconfined and polarises.
* **BC fusion**: two aggregates, one B and one C, glue together instead of
  rounding into a sphere.

The living report with the current state is
**[reports/REPORT.md](reports/REPORT.md)**. The physics and every parameter are in
[reports/MODEL.md](reports/MODEL.md), and how each final figure was made is in
[reports/FIGURES.md](reports/FIGURES.md). The key reasons behind the modelling choices are summarised in
[reports/APPENDIX.md](reports/APPENDIX.md).

**The parameter set** (`PSETS["final"]` in the notebook, Oct 8): the July/FUSION mechanics with B stiffness
`f_rep` BB = 4.5, C stiffness `f_rep` CC = 6 and a shorter B–B adhesion range `mu` BB = 1.9 (patterns in
slabs of 3+ layers), friction with the plates `lambda_wall` = 1.5, rigid plates, and a sharp B core. The BC
entries are the FUSION values, unchanged. The Oct 7 set (`mu` BB = 2, patterns only in 1–2 layers) is kept as
`PSETS["final_oct7"]`, with its figures in `code-claude/figures-oct7/`.

## Layout

```
code-original/        the original code, as committed in f7fd270
  squeezing-gastruloids.ipynb
  july-laptop/        your July notebook from the laptop (committed Oct 8)
  plots/              original figures and movies
  preamble/           earlier function notebooks
  sotiris/            Sotiris' fusion notebook and figure (reference only)
code-claude/          Claude's code, built on the original notebook
  squeezing-gastruloids-claude.ipynb
  scripts/            batch runs (library exported from the notebook); README.md says what each one is for
  figures/            the final 3D figures, movies and the presentation PDF
    data/             the data behind every final figure (cells: x, y, z, r, cell_state)
  figures-oct7/       the same figures for the Oct 7 parameter set (kept)
  results/            raw simulation runs (large, gitignored)
  july-2026/          superseded July code, figures and runs
reports/
  REPORT.md           living report: constraints, decisions, results
  MODEL.md            the model: equations, physics, parameters
  FIGURES.md          how each final figure was made (initial condition, procedure, data)
  APPENDIX.md         the key reasons behind the modelling choices (short)
  2026-07 july/       the July reports (superseded), with the LaTeX report in report/
  2026-10 october/    the October LaTeX report
  2026-11 november/   the November LaTeX report (in progress)
claude-answers/       Claude's answers, as markdown with images
experiments/         local only, gitignored
  pattern-examples/   experimental pattern examples
  presentations/      experimental presentations
Project.toml, Manifest.toml   pinned Julia environment
```

## Setup

```julia
using Pkg
Pkg.add(url = "https://github.com/dsb-lab/CellBasedModels.jl")
```

then, from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Project.toml` and `Manifest.toml` are a snapshot of the environment the results
were produced in. If you would rather install by hand:

```
IJulia
https://github.com/dsb-lab/CellBasedModels.jl
DifferentialEquations, Distributions, Random, Statistics, GLMakie, CairoMakie,
Printf, Dates, Glob, NBInclude, MathTeXEngine, JLD2, ProgressMeter
CSV, WriteVTK, DataFrames, Distances, Clustering, GeometryBasics, Colors,
StaticArrays
```

## A note on saving runs

`jldsave("com.jld2"; com)` does **not** round-trip. `ABM()` evaluates the
generated step functions as anonymous closures in `Main`, and the saved
`Community` refers to them by their generated names, so reloading in a
different session fails with

```
MethodError: Cannot `convert` an object of type
             JLD2.ReconstructedSingleton{Symbol("#19#20")} to ... Function
```

Use `save_frames` / `load_frames` (plain arrays of positions and cell types)
instead. They load in any session, even with the model undefined.
