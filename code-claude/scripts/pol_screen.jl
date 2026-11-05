# The polarisation screens of Nov 3 (results/stage5/pol_nov3 and pol2_nov3; claude-answers/10 §2.4): for every set,
# the unconfined runs' polarisation over time, against the current set (phase_D).
#
#   julia --project=. -t 16 code-claude/scripts/pol_screen.jl
#
# Per run: when P first reaches 0.3 (hours from the start of differentiation), P at fB 0.12, the
# compression at squeezing (median distance to the nearest neighbour, in r), the class read
# (class_read), and the shape at fB 0.12: where the B centroid sits along the long axis (PC1) of
# the aggregate, "end" = |B centroid - middle| / half the length (about 1 for a B pole at one end,
# about 0 for a belt of B at the waist between two C lobes, which also gives a large P). Per set: the fusion assay. Writes results/stage5/pol_nov3/summary.txt and
# claude-answers/img/10-pol-screen.png (P(t), one panel per N0).

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5 = joinpath(@__DIR__, "..", "results", "stage5")
const T_RELAX = 25.0
const BASE = "july_frB4.5_frC6_muB1.9"

# both screens: pol_nov3 (C side, stiffness, half friction) and pol2_nov3 (B and C both more
# cohesive; weaker BC attraction as a diagnostic), when it exists
scr = filter(r -> r["ok"] == true, vcat([read_rows(joinpath(S5, f)) for f in ("pol_nov3", "pol2_nov3") if isdir(joinpath(S5, f))]...))
base = filter(r -> r["ok"] == true && r["kind"] == "free" && String(r["cand"]) == BASE && r["seed"] == 1 &&
				   r["N0"] in (600, 1200, 2400, 4800), read_rows(joinpath(S5, "phase_D")))
runs = with_readout!(vcat(base, filter(r -> r["kind"] == "free", scr)))     # once, before the threads (one cache file per folder)
res = Vector{Any}(undef, length(runs))
Threads.@threads for i in eachindex(runs)
	r = runs[i]; d = load_frames(joinpath(r["dir"], String(r["file"])))
	k0 = frame_at(d, 0.6); ks = k0:length(d.x)
	t = [(d.t[k] - T_RELAX) / 15 for k in ks]
	P = [exp_polarisation(d.x[k], d.y[k], d.z[k], d.cell_state[k]).P for k in ks]
	fB = [count(==(2), d.cell_state[k]) / length(d.cell_state[k]) for k in ks]
	x, y, z = d.x[k0], d.y[k0], d.z[k0]
	nb = neighbour_lists(x, y, z; rad = 3.0)
	nn = median([minimum([sqrt((x[a]-x[b])^2 + (y[a]-y[b])^2 + (z[a]-z[b])^2) for b in nb[a]]; init = Inf) for a in eachindex(x)])
	j12 = findfirst(<=(0.12), fB)
	kk = j12 === nothing ? length(d.x) : ks[j12]
	X = hcat(d.x[kk] .- mean(d.x[kk]), d.y[kk] .- mean(d.y[kk]), d.z[kk] .- mean(d.z[kk]))
	bi = findall(==(2), d.cell_state[kk])
	e1 = eigen(Symmetric(X' * X ./ size(X, 1))).vectors[:, end]; pr = X * e1
	mid = (maximum(pr) + minimum(pr)) / 2; half = (maximum(pr) - minimum(pr)) / 2
	pole = isempty(bi) ? NaN : abs(mean(pr[bi]) - mid) / half
	res[i] = (cand = String(r["cand"]), N0 = Int(r["N0"]), t = t, P = P, fB = fB, nn = nn,
			  t03 = (i3 = findfirst(>=(0.3), P); i3 === nothing ? NaN : t[i3]), P12 = j12 === nothing ? NaN : P[j12], pole = pole,
			  cls = phase_class(String(r["class_read"])))
end
cands = vcat([BASE], sort(unique(r.cand for r in res if r.cand != BASE)))
short(c) = c == BASE ? "now" : replace(c, BASE * "_" => "", "july_" => "", "_muB1.9" => "", BASE => "now")
lines_ = String["set | N0 600 | 1200 | 2400 | 4800  (each: hours to P >= 0.3 / P at fB 0.12 / class / end = B centroid along the long axis (pole ~1, belt ~0); NN = nearest-neighbour distance at squeezing, r)", ""]
for c in cands
	cells = String[]
	for N in (600, 1200, 2400, 4800)
		m = filter(r -> r.cand == c && r.N0 == N, res)
		push!(cells, isempty(m) ? "–" : @sprintf("%s / %.2f / %s / end %.2f / NN %.2f", isnan(m[1].t03) ? "never" : @sprintf("%.0f h", m[1].t03),
												m[1].P12, lowercase(m[1].cls)[1:3], m[1].pole, m[1].nn))
	end
	push!(lines_, @sprintf("%-22s | ", short(c)) * join(cells, " | "))
end
push!(lines_, "")
for r in filter(r -> r["kind"] == "fusion", scr)
	push!(lines_, @sprintf("fusion %-22s BB %.2f CC %.2f BC %.2f", short(String(r["cand"])), r["bb24"], r["cc24"], r["bc24"]))
end
write(joinpath(S5, "pol_nov3", "summary.txt"), join(lines_, "\n") * "\n")
println(join(lines_, "\n"))

fig = Figure(size = (1900, 520), figure_padding = 20, fonts = TEX_FONTS)
cols = Makie.wong_colors()
pal = Dict(c => (c == BASE ? :black : cols[mod1(i, 7)], i > 8 ? :dash : :solid) for (i, c) in enumerate(cands))
for (j, N) in enumerate((600, 1200, 2400, 4800))
	ax = Axis(fig[1, j]; xlabel = L"t\ \mathrm{(h)}", ylabel = j == 1 ? L"P" : "", title = Makie.LaTeXString("\$N_0 = $(N)\$"))
	hlines!(ax, [0.3]; color = (:black, 0.4), linestyle = :dot)
	for r in filter(r -> r.N0 == N, res)
		lines!(ax, r.t, r.P; color = pal[r.cand][1], linestyle = pal[r.cand][2], linewidth = r.cand == BASE ? 3 : 1.6, label = short(r.cand))
	end
	ylims!(ax, -0.02, 0.85)
end
Legend(fig[1, 5], [LineElement(color = pal[c][1], linestyle = pal[c][2], linewidth = c == BASE ? 3 : 1.6) for c in cands], short.(cands);
	   framevisible = false, labelsize = 16)
save(joinpath(@__DIR__, "..", "..", "claude-answers", "img", "10-pol-screen.png"), fig)
println("wrote 10-pol-screen.png")
