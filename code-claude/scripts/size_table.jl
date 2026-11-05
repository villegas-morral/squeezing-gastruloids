# Number of cells against size at the moment of squeezing (Nov 3 2026; your To-Do in the
# November report: "a comparison table with diameter and number of cells accounting for some
# overlapping").
#
#   julia --project=. -t 8 code-claude/scripts/size_table.jl
#
# From the seeding frame (fB 0.6, just before differentiation starts) of the final set's runs in
# results/stage5/phase_D:
#   sphere   unconfined runs: diameter 2 (sqrt(5/3) Rg + r), Rg the radius of gyration of the cell
#            centres (exact for a uniform ball), r = 1 = 5 um
#   pancake  confined runs: diameter 2 sqrt(2 <rho^2>) + 2r (uniform disc; the D column)
# and the packing: the volume fraction phi = N (4/3 pi r^3) / volume. Cells overlap (overlap
# factor 0.3, neighbours sit at ~1.8 r instead of 2 r), so phi can exceed the 0.64 of random close
# packing of hard spheres. Prints a markdown table.

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie                 # the analysis library needs a Makie backend loaded
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5 = joinpath(@__DIR__, "..", "results", "stage5")
const CAND = "july_frB4.5_frC6_muB1.9"
const R_UM = 5.0
rows = filter(r -> r["ok"] == true && r["seed"] == 1 && ((r["kind"] == "free" && r["cand"] == CAND) ||
			  (r["kind"] == "confined" && r["cand"] == CAND * "@lw1.5" && r["layers"] in (1, 2, 3, 4, 8))),
			  read_rows(joinpath(S5, "phase_D")))
res = Vector{Any}(undef, length(rows))
Threads.@threads for i in eachindex(rows)
	r = rows[i]; d = load_frames(joinpath(r["dir"], String(r["file"]))); k = frame_at(d, 0.6)
	x, y, z = d.x[k], d.y[k], d.z[k]; n = length(x)
	X = x .- mean(x); Y = y .- mean(y); Z = z .- mean(z)
	if r["kind"] == "free"
		Rg = sqrt(mean(X.^2 .+ Y.^2 .+ Z.^2)); Rb = sqrt(5/3) * Rg + 1.0
		res[i] = (kind = "sphere", L = 0, N0 = r["N0"], n = n, D = 2Rb, phi = n / Rb^3)       # (4/3 pi r^3)/(4/3 pi R^3) = 1/R^3
	else
		Rd = sqrt(2mean(X.^2 .+ Y.^2)) + 1.0; h = maximum(Z) - minimum(Z) + 2.0
		res[i] = (kind = "pancake", L = r["layers"], N0 = r["N0"], n = n, D = 2Rd, phi = n * (4/3) / (Rd^2 * h))   # disc pi Rd^2 h
	end
end
Ns = sort(unique(r.N0 for r in res if r.N0 in (150, 300, 600, 1200, 2400, 4800)))
get_(kind, L, N) = (m = filter(r -> r.kind == kind && r.L == L && r.N0 == N, res); isempty(m) ? nothing : m[1])
println("| N0 | sphere diameter (µm) | φ sphere | pancake diameter (µm): 1 layer | 2 layers | 3 layers | 4 layers | 8 layers |")
println("|---|---|---|---|---|---|---|---|")
for N in Ns
	s = get_("sphere", 0, N)
	cells = [(p = get_("pancake", L, N); p === nothing ? "–" : @sprintf("%.0f", p.D * R_UM)) for L in (1, 2, 3, 4, 8)]
	@printf("| %d | %s | %s | %s |\n", N, s === nothing ? "–" : @sprintf("%.0f", s.D * R_UM), s === nothing ? "–" : @sprintf("%.2f", s.phi),
			join(cells, " | "))
end
# a fit: sphere diameter against N0 (D ∝ N0^(1/3) if the packing is constant)
sp = sort(filter(r -> r.kind == "sphere", res); by = r -> r.N0)
for r in sp
	@printf("sphere N0 %4d (N %4d at squeezing): D %.1f r = %.0f um, phi %.2f, D/N^(1/3) %.2f r\n", r.N0, r.n, r.D, r.D * R_UM, r.phi, r.D / r.n^(1/3))
end
pc = filter(r -> r.kind == "pancake", res)
@printf("pancakes: phi %.2f-%.2f (median %.2f)\n", minimum(r.phi for r in pc), maximum(r.phi for r in pc), median([r.phi for r in pc]))
