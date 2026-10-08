# Figures of the Oct 7 set (kept)

These are the final figures of the set chosen on Oct 7: `PSETS["final_oct7"]`, "july_frB4.5_frC6@lw1.5" in
the stage 5 scripts. In that set, slabs of 3 or more layers always end as one central blob.

On Oct 8 they were redrawn with:
- the four phase-space versions (layers or pancake diameter, against N0 or log₁₀ N0)
- the presentation PDF

Then `../figures/` was regenerated for the new set (`mu` BB 1.9), which also patterns in thick slabs. See
`reports/REPORT.md` (stage 7) and `claude-answers/04-july-thick-slabs-phase-plots-pdf.md`.

Everything here can be redrawn from `data/` alone:

```julia
include("code-claude/scripts/gastruloids_analysis.jl"); replot_figures("code-claude/figures-oct7")
```

The runs behind them are in `results/stage5/phase_B/` and `results/stage5/tune/` (gitignored).
