# Shapes of the polarisation screen's unconfined runs (Nov 3 2026, claude-answers/10 §2): a
# mid-plane cross-section of each aggregate at its readout frame (fB 0.12), rotated so that the
# long axis (PC1) is horizontal: a pole shows B at one end, a dumbbell B at the waist. CairoMakie,
# no desktop needed.
#
#   julia --project=. -t 8 code-claude/scripts/pol_shapes.jl
#
# -> claude-answers/img/10-pol-shapes.png

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5 = joinpath(@__DIR__, "..", "results", "stage5")
const BASE = "july_frB4.5_frC6_muB1.9"
# (folder, candidate, row label)
const SETS = [("phase_D", BASE, "now"), ("pol_nov3", BASE * "@lam0.5", "\$\\lambda \\times 0.5\$"),
			  ("pol2_nov3", BASE * "_faB7.5_faC3.6", "BB 7.5, CC 3.6"), ("pol_nov3", BASE * "_faC3.6", "CC 3.6"),
			  ("pol2_nov3", BASE * "_faB8.5_faC4.2", "BB 8.5, CC 4.2")]
const NS = (600, 1200, 2400, 4800)

function section(file)
	d = load_frames(file)
	fB = [count(==(2), s) / length(s) for s in d.cell_state]
	k0 = frame_at(d, 0.6)                                                 # differentiation starts here
	k = something(findfirst(k -> k >= k0 && fB[k] <= 0.12, eachindex(fB)), length(d.x))
	X = hcat(d.x[k] .- mean(d.x[k]), d.y[k] .- mean(d.y[k]), d.z[k] .- mean(d.z[k]))
	E = eigen(Symmetric(X' * X ./ size(X, 1))).vectors[:, end:-1:1]      # PC1, PC2, PC3
	Y = X * E
	bi = d.cell_state[k] .== 2
	mean(Y[bi, 1]) < 0 && (Y[:, 1] .*= -1)                                # B to the right
	keep = abs.(Y[:, 3]) .< 1.5
	(Y = Y[keep, :], B = bi[keep], t = (d.t[k] - 25.0) / 15)
end

jobs = [(i, j) for i in eachindex(SETS), j in eachindex(NS)]
res = Dict{Tuple{Int,Int},Any}()
lk = ReentrantLock()
Threads.@threads for (i, j) in jobs
	folder, cand, _ = SETS[i]
	m = filter(r -> r["ok"] == true && r["kind"] == "free" && String(r["cand"]) == cand && r["N0"] == NS[j] && r["seed"] == 1,
			   read_rows(joinpath(S5, folder)))
	isempty(m) && continue
	s = section(joinpath(m[1]["dir"], String(m[1]["file"])))
	lock(lk) do; res[(i, j)] = s; end
end

fig = Figure(size = (1500, 1550), figure_padding = 15, fonts = TEX_FONTS)
for (i, (_, _, name)) in enumerate(SETS)
	Label(fig[i, 0], Makie.LaTeXString(name); rotation = pi/2, fontsize = 20, tellheight = false)
	for (j, N) in enumerate(NS)
		ax = Axis(fig[i, j]; aspect = DataAspect(), title = i == 1 ? Makie.LaTeXString("\$N_0 = $(N)\$") : "", titlesize = 20)
		hidedecorations!(ax); hidespines!(ax)
		haskey(res, (i, j)) || (text!(ax, 0, 0; text = "running", align = (:center, :center), color = :grey); continue)
		s = res[(i, j)]
		R = maximum(abs.(s.Y[:, 1:2])) + 1
		c = .!s.B
		scatter!(ax, s.Y[c, 1], s.Y[c, 2]; color = color_map[3], markersize = 2.0, markerspace = :data, strokewidth = 0.3, strokecolor = (:black, 0.25))
		scatter!(ax, s.Y[s.B, 1], s.Y[s.B, 2]; color = color_map[2], markersize = 2.0, markerspace = :data, strokewidth = 0.3, strokecolor = (:black, 0.25))
		limits!(ax, -R, R, -R, R)
		text!(ax, -R, -R; text = @sprintf("%.0f h", s.t), align = (:left, :bottom), fontsize = 14, color = :grey30)
	end
end
out = joinpath(@__DIR__, "..", "..", "claude-answers", "img", "10-pol-shapes.png")
save(out, fig)
println("wrote ", out)
