# The night grid of Nov 3 (claude-answers/10 §2.5): do the sets that polarise large aggregates
# keep the confined patterns? Six phase spaces on the same coarse map (layers 1, 2, 3, 4, 8 +
# unconfined; N0 150-2400, + 4800 at 3, 4 and 8 layers and 9600 at 3 and 4 layers where run;
# seed 1), read with the figures' readout (class_read, fB 0.12).
#
#   julia --project=. -t 8 code-claude/scripts/grid_compare_b.jl
#
# Sets: now (phase_D), f_att BB 6.1 + 13 h (grid_nov3, the grid winner of 09 §2), friction x 0.5
# (lam_nov3 + pol_nov3), f_att BB 7.5 + CC 3.6 with and without 13 h (grid_nov3b, 9600 in
# grid_nov3c), and the half-way set f_att BB 7 + CC 3.3 (grid_nov3c, free in grid_nov3b). Writes
# claude-answers/img/10-grid-b.png and results/stage5/grid_nov3b/compare.txt.

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie                 # 2D charts only
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5  = joinpath(@__DIR__, "..", "results", "stage5")
const IMG = joinpath(@__DIR__, "..", "..", "claude-answers", "img")
const B   = "july_frB4.5_frC6_muB1.9"
const LAYERS = (1.0, 2.0, 3.0, 4.0, 8.0)
const NS = (150, 300, 600, 1200, 2400, 4800, 9600)
inmap(L, N) = (L in LAYERS && N in NS[1:5]) || (L in (3.0, 4.0, 8.0) && N == 4800) || (L in (3.0, 4.0) && N == 9600)

const DIRS = ("phase_D", "tune", "grid_nov3", "phase_F", "grid_nov3b", "grid_nov3c", "lam_nov3", "pol_nov3")
all_ = Dict(d => (isdir(joinpath(S5, d)) ? with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, d)))) : Dict{String,Any}[]) for d in DIRS)
# (short name, title, confined candidate, unconfined candidate, folders)
const S75 = B * "_faB7.5_faC3.6"
const S70 = B * "_faB7_faC3.3"
SETS = [("now", L"\mathrm{now}", B * "@lw1.5", B, ["phase_D", "tune"]),
		("6.1, 13 h", L"f_\mathrm{att}^{BB} = 6.1,\ \tau_\mathrm{div} = 13\,\mathrm{h}", B * "_faB6.1@lw1.5@tdiv195", B * "_faB6.1@tdiv195", ["grid_nov3", "phase_F"]),
		("lambda x 0.5", L"\lambda \times 0.5\ \mathrm{(diagnostic)}", B * "@lw1.5@lam0.5", B * "@lam0.5", ["lam_nov3", "pol_nov3"]),
		("7.5/3.6", L"f_\mathrm{att}^{BB} = 7.5,\ f_\mathrm{att}^{CC} = 3.6", S75 * "@lw1.5", S75, ["grid_nov3b", "grid_nov3c"]),
		("7.5/3.6, 13 h", L"f_\mathrm{att}^{BB} = 7.5,\ f_\mathrm{att}^{CC} = 3.6,\ \tau_\mathrm{div} = 13\,\mathrm{h}", S75 * "@lw1.5@tdiv195", S75 * "@tdiv195", ["grid_nov3b"]),
		("7/3.3", L"f_\mathrm{att}^{BB} = 7,\ f_\mathrm{att}^{CC} = 3.3", S70 * "@lw1.5", S70, ["grid_nov3b", "grid_nov3c"])]

function pick_rows(i)
	_, _, cc, fc, dirs = SETS[i]
	rs = vcat([all_[d] for d in dirs]...)
	conf = filter(r -> r["kind"] == "confined" && String(r["cand"]) == cc && r["seed"] == 1 && inmap(Float64(r["layers"]), Int(r["N0"])), rs)
	free = filter(r -> r["kind"] == "free" && String(r["cand"]) == fc && r["seed"] == 1 && Int(r["N0"]) in NS, rs)
	u(rs, k) = [first(filter(r -> k(r) == x, rs)) for x in unique(k.(rs))]     # one row per point
	(u(conf, r -> (r["layers"], r["N0"])), u(free, r -> r["N0"]))
end

fig = Figure(size = (2400, 1800), figure_padding = 20, fonts = TEX_FONTS)
abbr = Dict("POLARISED" => "o", "CENTRAL" => "c", "PATTERN" => "X", "OTHER" => "n")
lines_ = String["o Polarisation, c B-centroid, X Clusters (seed 1, readout class_read at fB 0.12)", ""]
for (i, s) in enumerate(SETS)
	conf, free = pick_rows(i)
	cr = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), class = String(r["class_read"]), late = String(r["class_read_late"])) for r in conf]
	fr = [(N0 = Int(r["N0"]), class = String(r["class_read"]), late = String(r["class_read_late"])) for r in free]
	isempty(cr) && continue
	phase_space!(fig[1 + (i - 1) ÷ 3, 1 + (i - 1) % 3], cr; free = fr, y = :layers_lin, pies = true, title = s[2], fontsize = 22)
	push!(lines_, "== $(s[1]) ($(length(conf)) confined, $(length(free)) free)")
	for L in sort(unique(Float64(r["layers"]) for r in conf))
		push!(lines_, @sprintf("  %-2g  ", L) * join([@sprintf("%d:%s", r["N0"], abbr[phase_class(String(r["class_read"]))])
												 for r in sort(filter(r -> r["layers"] == L, conf); by = r -> r["N0"])], "  "))
	end
	push!(lines_, "  free " * join([@sprintf("%d:%s", r["N0"], abbr[phase_class(String(r["class_read"]))]) for r in sort(free; by = r -> r["N0"])], "  "))
end
els = [PolyElement(color = PHASE_COLORS[c]) for c in ("POLARISED", "CENTRAL", "PATTERN")]
Legend(fig[3, 1:3], els, [PHASE_NAMES[c] for c in ("POLARISED", "CENTRAL", "PATTERN")]; orientation = :horizontal,
	   framevisible = false, labelsize = 26, tellwidth = false, tellheight = true)
save(joinpath(IMG, "10-grid-b.png"), fig)
# the extra unconfined seeds and the half-way set
for r in sort(filter(r -> r["kind"] == "free" && r["seed"] > 1, all_["grid_nov3b"]); by = r -> (String(r["cand"]), r["N0"], r["seed"]))
	push!(lines_, @sprintf("free %s N0 %d seed %d: %s", r["cand"], r["N0"], r["seed"], phase_class(String(r["class_read"]))))
end
for r in filter(r -> r["kind"] == "fusion", all_["grid_nov3b"])
	push!(lines_, @sprintf("fusion %s: BB %.2f CC %.2f BC %.2f", r["cand"], r["bb24"], r["cc24"], r["bc24"]))
end
mkpath(joinpath(S5, "grid_nov3b"))
write(joinpath(S5, "grid_nov3b", "compare.txt"), join(lines_, "\n") * "\n")
println(join(lines_, "\n"))
