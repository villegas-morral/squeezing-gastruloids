# Nov 4: plate friction on the cohesive sets (claude-answers/10 §2.8). Plate friction acts on confined runs
# only; does more of it bring the 3-4-layer clusters of BB 7.5 / CC 3.6 and BB 7 / CC 3.3 back to smaller N0?
# Seed 1, readout class_read at fB 0.12; friction 1.5 (the current value) from grid_nov3b/c, 3 and 5 from
# wall_nov4. Prints per point the class for each friction, the B looseness S and the time to fB 0.2.
#
#   julia --project=. code-claude/scripts/wall_compare.jl
#
# Writes results/stage5/wall_nov4/compare.txt and claude-answers/img/10-wall.png (the readout frames).

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5  = joinpath(@__DIR__, "..", "results", "stage5")
const IMG = joinpath(@__DIR__, "..", "..", "claude-answers", "img")
const B   = "july_frB4.5_frC6_muB1.9"
const DIRS = ("phase_D", "tune", "grid_nov3b", "grid_nov3c", "wall_nov4")
rows = vcat([with_readout!(filter(r -> r["ok"] == true && r["kind"] == "confined" && r["seed"] == 1, read_rows(joinpath(S5, d))))
			 for d in DIRS if isdir(joinpath(S5, d))]...)
# (name, base candidate, frictions)
SETS = [("now", B, ("1.5",)), ("7 / 3.3", B * "_faB7_faC3.3", ("1.5", "3")), ("7.5 / 3.6", B * "_faB7.5_faC3.6", ("1.5", "3", "5"))]
abbr = Dict("POLARISED" => "o", "CENTRAL" => "c", "PATTERN" => "X", "OTHER" => "n")
find_row(cand, L, N) = (i = findfirst(r -> String(r["cand"]) == cand && Float64(r["layers"]) == L && Int(r["N0"]) == N, rows);
						i === nothing ? nothing : rows[i])
pts = [(4.0, 4800), (4.0, 2400), (4.0, 1200), (3.0, 4800), (3.0, 2400), (3.0, 1200), (3.0, 600), (8.0, 2400), (8.0, 1200),
	   (2.0, 2400), (2.0, 1200), (2.0, 600), (2.0, 300), (1.0, 2400), (1.0, 600), (1.0, 150)]
out = String["o Polarisation, c B-centroid, X Clusters (seed 1, class_read at fB 0.12); S = fraction of B among the",
			 "neighbours of a B cell at fB 0.2 (low = loose); t = hours to fB 0.2 from the start of differentiation", ""]
push!(out, @sprintf("%-12s", "layers, N0") * join([@sprintf("%-22s", "$(s[1]) @lw$f") for s in SETS for f in s[3]]))
for (L, N) in pts
	line = @sprintf("%-12s", @sprintf("%g, %d", L, N))
	for (nm, base, fs) in SETS, f in fs
		r = find_row("$(base)@lw$(f)", L, N)
		line *= r === nothing ? @sprintf("%-22s", "–") :
				@sprintf("%-22s", @sprintf("%s  S %.2f  t %.1f", abbr[phase_class(String(r["class_read"]))], r["S"], (Float64(r["t_02"]) - 25.0) / 15))
	end
	push!(out, line)
end
mkpath(joinpath(S5, "wall_nov4"))
write(joinpath(S5, "wall_nov4", "compare.txt"), join(out, "\n") * "\n")
println(join(out, "\n"))

# readout frames of the thick-slab points, rows = the five set/friction combinations of the two cohesive sets
fr_pts = [(2.0, 1200), (3.0, 1200), (3.0, 2400), (3.0, 4800), (4.0, 2400), (4.0, 4800)]
combos = [(nm, base, f) for (nm, base, fs) in SETS[2:3] for f in fs]
g = Figure(size = (330 * length(fr_pts), 360 * length(combos)), fonts = TEX_FONTS)
for (i, (nm, base, f)) in enumerate(combos), (j, (L, N)) in enumerate(fr_pts)
	r = find_row("$(base)@lw$(f)", L, N)
	r === nothing && continue
	d = load_frames(joinpath(r["dir"], String(r["file"]))); k = frame_at(d, 0.12)
	c = (x = d.x[k], y = d.y[k], z = d.z[k], r = d.r[k], cell_state = d.cell_state[k])
	draw_cells!(g[i, j], c; view = :top, labelsize = 15,
				label = @sprintf("%s, friction %s\n%g layers, N0 %d: %s", nm, f, L, N, PHASE_NAMES[phase_class(String(r["class_read"]))]))
end
save(joinpath(IMG, "10-wall.png"), g)
