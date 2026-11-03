# Equal division times (13 h for B and C, November 2026 test): a sketch of the phase space.
#
#   julia --project=. -t 8 code-claude/scripts/tdiv_phase.jl
#
# Reads results/stage5/tdiv (tau_div = 195 = 13 h for every cell, seed 1; layers 1, 2, 3,
# 4, 8 and unconfined x N0 150-2400) and the same points of results/stage5/phase_D (tau_div
# B 20 h, C 50 h; seed 1), and draws both charts side by side with the figures' readout
# (class_exp50_02) -> claude-answers/img/07-tdiv-phase-space.png.

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5   = joinpath(@__DIR__, "..", "results", "stage5")
const OUT  = joinpath(@__DIR__, "..", "..", "claude-answers", "img", "07-tdiv-phase-space.png")
const CAND = "july_frB4.5_frC6_muB1.9"
const CCOL = "class_exp50_02"

new = with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, "tdiv"))))
pts = Set((r["kind"] == "free" ? 0.0 : Float64(r["layers"]), Int(r["N0"])) for r in new)
old = with_readout!(filter(r -> r["ok"] == true && r["seed"] == 1 &&
						   ((r["kind"] == "confined" && r["cand"] == CAND * "@lw1.5") || (r["kind"] == "free" && r["cand"] == CAND)),
						   vcat(read_rows(joinpath(S5, "phase_D")), read_rows(joinpath(S5, "tune")))))
filter!(r -> (r["kind"] == "free" ? 0.0 : Float64(r["layers"]), Int(r["N0"])) in pts, old)

conf(rs) = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), class = String(r[CCOL]), late = String(r["class_exp50_012"]))
			for r in rs if r["kind"] == "confined"]
free(rs) = [(N0 = Int(r["N0"]), class = String(r[CCOL]), late = String(r["class_exp50_012"])) for r in rs if r["kind"] == "free"]

fig = Figure(size = (1700, 820), figure_padding = 20, fonts = TEX_FONTS)
for (j, (rs, ttl)) in enumerate(((old, L"\tau_\mathrm{div}^B = 20\,\mathrm{h},\ \tau_\mathrm{div}^C = 50\,\mathrm{h}\ \mathrm{(now)}"),
								 (new, L"\tau_\mathrm{div}^B = \tau_\mathrm{div}^C = 13\,\mathrm{h}")))
	phase_space!(fig[1, j], conf(rs); free = free(rs), title = ttl, pies = true)
end
els = [PolyElement(color = PHASE_COLORS[c]) for c in ("POLARISED", "CENTRAL", "PATTERN", "OTHER")]
Legend(fig[2, 1:2], els, [PHASE_NAMES[c] for c in ("POLARISED", "CENTRAL", "PATTERN", "OTHER")];
	   orientation = :horizontal, framevisible = false, labelsize = 24)
save(OUT, fig)
for (rs, nm) in ((old, "now"), (new, "13 h"))
	println(nm, ": ", join(sort([@sprintf("%s/%d %s", r["kind"] == "free" ? "free" : @sprintf("%g", r["layers"]), r["N0"], phase_class(String(r[CCOL])))
								 for r in rs]), ", "))
end
println("wrote ", OUT)
