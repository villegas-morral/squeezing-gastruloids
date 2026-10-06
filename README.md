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
**[reports/REPORT.md](reports/REPORT.md)**.

## Layout

```
code-original/        the original code, as committed in f7fd270
  squeezing-gastruloids.ipynb
  plots/              original figures and movies
  preamble/           earlier function notebooks
  sotiris/            Sotiris' fusion notebook and figure (reference only)
code-claude/          Claude's code, built on the original notebook
  squeezing-gastruloids-claude.ipynb
  scripts/            batch runs (library exported from the notebook)
  july-2026/          superseded July code, figures and runs
reports/
  REPORT.md           living report: constraints, decisions, results
  july-2026/          the July reports (superseded)
claude-answers/       Claude's answers, as markdown with images
experiments/
  pattern-examples/   experimental pattern examples
  presentations/      experimental presentations (local only, gitignored)
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
