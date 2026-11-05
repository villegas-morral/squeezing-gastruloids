# The f_att BB x tau_div grid of Nov 3 (claude-answers/09 §2): six phase spaces on the same
# coarse map (layers 1, 2, 3, 4, 8 + unconfined; N0 150-2400, + 4800 at 3 and 4 layers;
# seed 1), read with the figures' readout (class_read, fB 0.12).
#
#   julia --project=. -t 8 code-claude/scripts/grid_compare.jl
#
# Runs: f_att BB 6.4 / now = phase_D, 5.8 / now = phase_E, 6.4 / 13 h = tdiv (+ grid_nov3 for
# N0 4800), everything else grid_nov3. Writes claude-answers/img/09-grid.png and
# results/stage5/grid_nov3/compare.txt (outcome per point, fusion).

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie                 # 2D charts only
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5  = joinpath(@__DIR__, "..", "results", "stage5")
const IMG = joinpath(@__DIR__, "..", "..", "claude-answers", "img")
const B   = "july_frB4.5_frC6_muB1.9"
const LAYERS = (1.0, 2.0, 3.0, 4.0, 8.0)
const NS = (150, 300, 600, 1200, 2400)
inmap(L, N) = (L in LAYERS && N in NS) || (L in (3.0, 4.0) && N == 4800)

rows_of(dir) = with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, dir))))
all_ = Dict(d => rows_of(d) for d in ("phase_D", "phase_E", "tdiv", "grid_nov3", "tune"))
# (name, f_att BB, tdiv suffix, folders)
SETS = [("6.4", "", ["phase_D", "tune"]), ("6.1", "_faB6.1", ["grid_nov3"]), ("5.8", "_faB5.8", ["phase_E"]),
		("6.4", "", ["tdiv", "grid_nov3"]), ("6.1", "_faB6.1", ["grid_nov3"]), ("5.8", "_faB5.8", ["grid_nov3"])]
TD = ["", "", "", "@tdiv195", "@tdiv195", "@tdiv195"]

function pick_rows(i)
	fa, suf, dirs = SETS[i]; td = TD[i]
	conf = filter(r -> r["kind"] == "confined" && String(r["cand"]) == B * suf * "@lw1.5" * td && r["seed"] == 1 &&
					   inmap(Float64(r["layers"]), Int(r["N0"])), vcat([all_[d] for d in dirs]...))
	free = filter(r -> r["kind"] == "free" && String(r["cand"]) == B * suf * td && r["seed"] == 1 && Int(r["N0"]) in NS,
				  vcat([all_[d] for d in dirs]...))
	# one row per point (a run can sit in two folders)
	u(rs, k) = [first(filter(r -> k(r) == x, rs)) for x in unique(k.(rs))]
	(u(conf, r -> (r["layers"], r["N0"])), u(free, r -> r["N0"]))
end

titles = [L"f_\mathrm{att}^{BB} = 6.4\ \mathrm{(now)},\ \tau^B_\mathrm{div} = 20\,\mathrm{h},\ \tau^C_\mathrm{div} = 50\,\mathrm{h}",
		  L"f_\mathrm{att}^{BB} = 6.1,\ \tau^B_\mathrm{div} = 20\,\mathrm{h},\ \tau^C_\mathrm{div} = 50\,\mathrm{h}",
		  L"f_\mathrm{att}^{BB} = 5.8,\ \tau^B_\mathrm{div} = 20\,\mathrm{h},\ \tau^C_\mathrm{div} = 50\,\mathrm{h}",
		  L"f_\mathrm{att}^{BB} = 6.4,\ \tau_\mathrm{div} = 13\,\mathrm{h}",
		  L"f_\mathrm{att}^{BB} = 6.1,\ \tau_\mathrm{div} = 13\,\mathrm{h}",
		  L"f_\mathrm{att}^{BB} = 5.8,\ \tau_\mathrm{div} = 13\,\mathrm{h}"]
fig = Figure(size = (2400, 1700), figure_padding = 20, fonts = TEX_FONTS)
abbr = Dict("POLARISED" => "o", "CENTRAL" => "c", "PATTERN" => "X", "OTHER" => "n")
lines_ = String["o Polarisation, c B-centroid, X Clusters (seed 1, readout class_read at fB 0.12)", ""]
for i in 1:6
	conf, free = pick_rows(i)
	cr = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), class = String(r["class_read"]), late = String(r["class_read_late"])) for r in conf]
	fr = [(N0 = Int(r["N0"]), class = String(r["class_read"]), late = String(r["class_read_late"])) for r in free]
	phase_space!(fig[1 + (i - 1) ÷ 3, 1 + (i - 1) % 3], cr; free = fr, y = :layers_lin, pies = true, title = titles[i], fontsize = 22)
	push!(lines_, "== f_att BB $(SETS[i][1]) $(TD[i] == "" ? "tau now" : "tau 13 h") ($(length(conf)) confined, $(length(free)) free)")
	for L in sort(unique(Float64(r["layers"]) for r in conf))
		push!(lines_, @sprintf("  %-2g  ", L) * join([@sprintf("%d:%s", r["N0"], abbr[phase_class(String(r["class_read"]))])
												 for r in sort(filter(r -> r["layers"] == L, conf); by = r -> r["N0"])], "  "))
	end
	push!(lines_, "  free " * join([@sprintf("%d:%s", r["N0"], abbr[phase_class(String(r["class_read"]))]) for r in sort(free; by = r -> r["N0"])], "  "))
end
els = [PolyElement(color = PHASE_COLORS[c]) for c in ("POLARISED", "CENTRAL", "PATTERN")]
Legend(fig[3, 1:3], els, [PHASE_NAMES[c] for c in ("POLARISED", "CENTRAL", "PATTERN")]; orientation = :horizontal,
	   framevisible = false, labelsize = 26)
save(joinpath(IMG, "09-grid.png"), fig)
# fusion of the new set
for r in filter(r -> r["kind"] == "fusion", all_["grid_nov3"])
	push!(lines_, @sprintf("fusion %s: BB %.2f CC %.2f BC %.2f", r["cand"], r["bb24"], r["cc24"], r["bc24"]))
end
write(joinpath(S5, "grid_nov3", "compare.txt"), join(lines_, "\n") * "\n")
println(join(lines_, "\n"))
