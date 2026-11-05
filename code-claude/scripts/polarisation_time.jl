# Why large unconfined aggregates do not polarise (Nov 3 2026, claude-answers/10): the
# polarisation P(t) and the B fraction of every unconfined run of a set, against the time from
# the start of differentiation, by N0.
#
#   julia --project=. -t 16 code-claude/scripts/polarisation_time.jl [folder cand out.png]
#
# Default: results/stage5/phase_D, the final set's unconfined runs ("july_frB4.5_frC6_muB1.9"),
# -> claude-answers/img/10-polarisation-time.png and the same name .csv.

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5 = joinpath(@__DIR__, "..", "results", "stage5")
const FOLDER = length(ARGS) >= 1 ? ARGS[1] : "phase_D"
const CAND = length(ARGS) >= 2 ? ARGS[2] : "july_frB4.5_frC6_muB1.9"
const OUT = length(ARGS) >= 3 ? ARGS[3] : joinpath(@__DIR__, "..", "..", "claude-answers", "img", "10-polarisation-time.png")
const T_RELAX = 25.0

rows = filter(r -> r["ok"] == true && r["kind"] == "free" && String(r["cand"]) == CAND, read_rows(joinpath(S5, FOLDER)))
sort!(rows; by = r -> (r["N0"], r["seed"]))
println(length(rows), " unconfined runs")
res = Vector{Any}(undef, length(rows))
Threads.@threads for i in eachindex(rows)
	r = rows[i]; d = load_frames(joinpath(r["dir"], String(r["file"])))
	k0 = frame_at(d, 0.6)
	ks = k0:length(d.x)
	t = [(d.t[k] - T_RELAX) / 15 for k in ks]
	P = [exp_polarisation(d.x[k], d.y[k], d.z[k], d.cell_state[k]).P for k in ks]
	fB = [count(==(2), d.cell_state[k]) / length(d.cell_state[k]) for k in ks]
	# radius of the aggregate (98th percentile of the distance from the centre) and depth of the
	# B centroid: how far it still has to go to reach the surface
	R = [quantile(sqrt.((d.x[k] .- mean(d.x[k])).^2 .+ (d.y[k] .- mean(d.y[k])).^2 .+ (d.z[k] .- mean(d.z[k])).^2), 0.98) for k in ks]
	res[i] = (N0 = r["N0"], seed = r["seed"], t = t, P = P, fB = fB, R = R)
end
fig = Figure(size = (1500, 560), figure_padding = 20, fonts = TEX_FONTS)
ax1 = Axis(fig[1, 1]; xlabel = L"t\ \mathrm{(h)}", ylabel = L"P", title = "polarisation over time")
ax2 = Axis(fig[1, 2]; xlabel = L"f_B", ylabel = L"P", xreversed = true, title = L"\mathrm{against\ the\ }B\mathrm{\ fraction}")
Ns = sort(unique(r.N0 for r in res))
cols = cgrad(:viridis, length(Ns); categorical = true)
for r in res
	c = cols[findfirst(==(r.N0), Ns)]
	lines!(ax1, r.t, r.P; color = c, label = Makie.LaTeXString("\$N_0 = $(r.N0)\$"))
	lines!(ax2, r.fB, r.P; color = c)
end
for ax in (ax1, ax2); hlines!(ax, [0.3]; color = :black, linestyle = :dash); end
vlines!(ax2, [0.2, 0.12]; color = (:grey, 0.6), linestyle = :dot)
axislegend(ax1; position = :lt, unique = true, framevisible = false)
save(OUT, fig)
open(replace(OUT, r"\.png$" => ".csv"), "w") do io
	println(io, "N0,seed,t_h,P,fB,R")
	for r in res, i in eachindex(r.t)
		@printf(io, "%d,%d,%.2f,%.4f,%.4f,%.3f\n", r.N0, r.seed, r.t[i], r.P[i], r.fB[i], r.R[i])
	end
end
# summary: time to reach P = 0.3, P at fB 0.2 and 0.12, max P
for r in res
	i3 = findfirst(>=(0.3), r.P)
	j2 = findfirst(<=(0.2), r.fB); j12 = findfirst(<=(0.12), r.fB)
	@printf("N0 %4d s%d: P>=0.3 at %s h | P(fB .2) %.2f  P(fB .12) %.2f  max P %.2f | R0 %.1f r\n", r.N0, r.seed,
			i3 === nothing ? "never" : @sprintf("%.0f", r.t[i3]), j2 === nothing ? NaN : r.P[j2],
			j12 === nothing ? NaN : r.P[j12], maximum(r.P), r.R[1])
end
