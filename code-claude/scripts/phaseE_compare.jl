# The phase space of f_att BB 5.8 (results/stage5/phase_E, Nov 2 2026) next to the final set's
# (phase_D), with the figures' readout (class_exp50_02).
#
#   julia --project=. -t 16 code-claude/scripts/phaseE_compare.jl
#
# Writes claude-answers/img/07-faB58-phase-space.png (the two charts, background + area lines),
# claude-answers/img/07-faB58-thick.png (3D readout frames of the thick-slab clusters) and
# results/stage5/phase_E/compare.txt (majority class per point, both sets, every seed).

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie                 # everything here is opaque (confined panels, cut at the mid-plane)
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5   = joinpath(@__DIR__, "..", "results", "stage5")
const IMG  = joinpath(@__DIR__, "..", "..", "claude-answers", "img")
const CCOL = "class_exp50_02"
const SETS = (("now (f_att BB 6.4)", "phase_D", "july_frB4.5_frC6_muB1.9"),
			  ("f_att BB 5.8", "phase_E", "july_frB4.5_frC6_muB1.9_faB5.8"))

rows = Dict{String,Any}()
for (name, dir, base) in SETS
	rs = with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, dir))))
	conf = with_diameter!(filter(r -> r["kind"] == "confined" && String(r["cand"]) == base * "@lw1.5", rs))
	free = filter(r -> r["kind"] == "free" && String(r["cand"]) == base, rs)
	if isempty(free)          # the final set's unconfined runs are in tune/
		free = filter(r -> r["kind"] == "free" && String(r["cand"]) == base,
					  with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, "tune")))))
		Ns = Set(r["N0"] for r in filter(r -> r["kind"] == "free" && String(r["cand"]) == SETS[2][3], rs))
	end
	rows[name] = (conf = conf, free = free)
end
# the same unconfined N0 values in both
Nf = Set(r["N0"] for r in rows[SETS[2][1]].free)
for (name, _, _) in SETS
	rows[name] = (conf = rows[name].conf, free = filter(r -> r["N0"] in Nf, rows[name].free))
end

pt(r) = (Float64(r["layers"]), Int(r["N0"]))
cls(r) = phase_class(String(r[CCOL]))
abbr = Dict("POLARISED" => "o", "CENTRAL" => "c", "PATTERN" => "X", "OTHER" => "n")
majority(v) = (u = unique(v); u[argmax([count(==(c), v) for c in u])])

# ---- the two charts ---------------------------------------------------------------------
fig = Figure(size = (1900, 900), figure_padding = 20, fonts = TEX_FONTS)
for (j, (name, _, _)) in enumerate(SETS)
	R = rows[name]
	cr = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), class = String(r[CCOL]), late = String(r["class_exp50_012"]))
		  for r in R.conf]
	fr = [(N0 = Int(r["N0"]), class = String(r[CCOL]), late = String(r["class_exp50_012"])) for r in R.free]
	cont = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), A = pi * r["D"]^2 * 25) for r in R.conf]
	phase_space!(fig[1, j], cr; free = fr, title = name, contours = cont)
end
els = Any[PolyElement(color = PHASE_COLORS[c]) for c in ("POLARISED", "CENTRAL", "PATTERN", "OTHER")]
push!(els, LineElement(color = (:black, 0.7), linewidth = 3))
Legend(fig[2, 1:2], els, vcat([PHASE_NAMES[c] for c in ("POLARISED", "CENTRAL", "PATTERN", "OTHER")], "pancake area at squeezing (µm²)");
	   orientation = :horizontal, framevisible = false, labelsize = 24)
save(joinpath(IMG, "07-faB58-phase-space.png"), fig)

# ---- per point ---------------------------------------------------------------------------
lines_ = String["o polarised, c B-centroid, X clusters, n other; one letter per seed; * = majority changed", ""]
A, B = rows[SETS[1][1]], rows[SETS[2][1]]
for L in sort(unique(first.(pt.(B.conf))))
	parts = String[]
	for N in sort(unique(last.(filter(p -> p[1] == L, pt.(B.conf)))))
		a = [cls(r) for r in A.conf if pt(r) == (L, N)]; b = [cls(r) for r in B.conf if pt(r) == (L, N)]
		ch = !isempty(a) && majority(a) != majority(b) ? "*" : ""
		push!(parts, @sprintf("%d:%s>%s%s", N, join(abbr[c] for c in a), join(abbr[c] for c in b), ch))
	end
	push!(lines_, @sprintf("%-3g %s", L, join(parts, "  ")))
end
fa = Dict(N => [cls(r) for r in A.free if r["N0"] == N] for N in Nf); fb = Dict(N => [cls(r) for r in B.free if r["N0"] == N] for N in Nf)
push!(lines_, "free " * join([@sprintf("%d:%s>%s", N, join(abbr[c] for c in fa[N]), join(abbr[c] for c in fb[N])) for N in sort(collect(Nf))], "  "))
# readout time, growth, thickness
for (name, _, _) in SETS
	R = rows[name]
	push!(lines_, @sprintf("%s: readout %.1f h on average (%.1f-%.1f); layers_eff/layers %.2f-%.2f",
		name, mean((r["t_02"] - 25) / 15 for r in R.conf), minimum((r["t_02"] - 25) / 15 for r in R.conf),
		maximum((r["t_02"] - 25) / 15 for r in R.conf), minimum(r["layers_eff"] / r["layers"] for r in R.conf),
		maximum(r["layers_eff"] / r["layers"] for r in R.conf)))
end
write(joinpath(S5, "phase_E", "compare.txt"), join(lines_, "\n") * "\n")
println(join(lines_, "\n"))

# ---- 3D readout frames of the thick-slab points (seed 1, both sets) -------------------------
pts = [(3, 1200), (3, 2400), (4, 2400), (4, 4800), (6, 2400)]
g = Figure(size = (360 * length(pts), 860), figure_padding = (15, 15, 15, 60))
for (j, p) in enumerate(pts), (i, (name, _, _)) in enumerate(SETS)
	rs = sort(filter(r -> pt(r) == p, rows[name].conf); by = r -> r["seed"]); isempty(rs) && continue
	r = rs[1]; d = load_frames(joinpath(r["dir"], String(r["file"]))); k = frame_at(d, 0.2)
	c = (x = d.x[k], y = d.y[k], z = d.z[k], r = d.r[k], cell_state = d.cell_state[k])
	draw_cells!(g[i, j], c; view = :top, labelsize = 15,
				label = @sprintf("%s\n%g layers, N0 %d, seed %d\n%s (N_eff %.1f)", name, p[1], p[2], r["seed"],
								 PHASE_NAMES[cls(r)], r["N_eff50"]))
end
save(joinpath(IMG, "07-faB58-thick.png"), g)
println("wrote ", IMG)
