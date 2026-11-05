# Nov 4: the large unconfined aggregates of the candidate sets (claude-answers/10 §2.6).
# Prints the time to P >= 0.3, P and the pole-vs-belt "end" at fB 0.12, and draws the readout frames: a one-cell
# slice through the centre containing the long axis (PC1, horizontal, B end to the right) and the B centroid,
# so an off-axis B lump shows too. Runs whose trajectory is not there yet (N0 9600 + λ × 0.5) are skipped.
#   julia --project=. code-claude/scripts/pol_large.jl
using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra, CairoMakie
const ROOT = joinpath(@__DIR__, "..", "..")
include(joinpath(@__DIR__, "gastruloids_analysis.jl")); include(joinpath(@__DIR__, "stage5_common.jl"))
const S5 = joinpath(ROOT, "code-claude", "results", "stage5")
const T_RELAX = 25.0
B = "july_frB4.5_frC6_muB1.9"
runs = [("now", "phase_D/com_free_N4800_$(B)_w0_s1.jld2"),
		("7 / 3.3, seed 1", "grid_nov3b/com_free_N4800_$(B)_faB7_faC3.3_w0_s1.jld2"),
		("7.5 / 3.6, seed 1", "grid_nov3b/com_free_N4800_$(B)_faB7.5_faC3.6_w0_s1.jld2"),
		("7.5 / 3.6, seed 2", "grid_nov3b/com_free_N4800_$(B)_faB7.5_faC3.6_w0_s2.jld2"),
		("7.5 / 3.6, seed 3", "lam_nov4/com_free_N4800_$(B)_faB7.5_faC3.6_w0_s3.jld2"),
		("7.5 / 3.6 + λ × 0.5, seed 1", "lam_nov4/com_free_N4800_$(B)_faB7.5_faC3.6@lam0.5_w0_lam0.5_s1.jld2"),
		("7.5 / 3.6 + λ × 0.5, seed 2", "lam_nov4/com_free_N4800_$(B)_faB7.5_faC3.6@lam0.5_w0_lam0.5_s2.jld2"),
		("N0 9600: 7.5 / 3.6", "grid_nov3c/com_free_N9600_$(B)_faB7.5_faC3.6_w0_s1.jld2"),
		("N0 9600: 7.5 / 3.6 + λ × 0.5", "lam_nov4/com_free_N9600_$(B)_faB7.5_faC3.6@lam0.5_w0_lam0.5_s1.jld2")]
cs = []
for (nm, f) in runs
	isfile(joinpath(S5, f)) || (println("not there yet: ", f); continue)
	d = load_frames(joinpath(S5, f))
	ks = frame_at(d, 0.6):length(d.x)
	P_t = [exp_polarisation(d.x[k], d.y[k], d.z[k], d.cell_state[k]).P for k in ks]
	i3 = findfirst(>=(0.3), P_t)
	k = frame_at(d, 0.12)
	X = hcat(d.x[k] .- mean(d.x[k]), d.y[k] .- mean(d.y[k]), d.z[k] .- mean(d.z[k]))
	V = eigen(Symmetric(X' * X ./ size(X, 1))).vectors
	e1 = V[:, 3]
	bi = findall(==(2), d.cell_state[k])
	u = vec(mean(X[bi, :]; dims = 1)); u .-= dot(u, e1) * e1   # B offset across the long axis
	e2 = norm(u) > 0.5 ? u / norm(u) : V[:, 2]
	Y = X * hcat(e1, e2, cross(e1, e2))              # long axis, the B side across it, normal of the slice
	Y[:, 1] .*= sign(mean(Y[bi, 1]) + eps())
	pr = Y[:, 1]; mid = (maximum(pr) + minimum(pr)) / 2; half = (maximum(pr) - minimum(pr)) / 2
	en = abs(mean(pr[bi]) - mid) / half
	P = exp_polarisation(d.x[k], d.y[k], d.z[k], d.cell_state[k]).P
	keep = findall(w -> abs(w) < 1.0, Y[:, 3])
	c = (x = -Y[keep, 2], y = Y[keep, 1], z = Y[keep, 3], r = d.r[k][keep], cell_state = d.cell_state[k][keep])
	t = (d.t[k] - T_RELAX) / 15
	t3 = i3 === nothing ? "never" : @sprintf("at %.0f h", (d.t[ks[i3]] - T_RELAX) / 15)
	@printf("%-30s P >= 0.3 %-8s  fB 0.12 at %.1f h: P %.2f, end %.2f, %d cells\n", nm, t3, t, P, en, length(d.x[k]))
	push!(cs, (nm, c, @sprintf("P ≥ 0.3 %s; P %.2f, end %.2f", t3, P, en), maximum(abs, Y[:, 1:2])))
end
bs = maximum(last.(cs)) + 1
nc = 4; nr = cld(length(cs), nc)
fig = Figure(size = (420 * nc, 470 * nr), fonts = TEX_FONTS)
for (i, (nm, c, s, _)) in enumerate(cs)
	draw_cells!(fig[cld(i, nc), mod1(i, nc)], c; view = :top, cut_top = false, label = "$nm\n$s", labelsize = 17, bs = bs)
end
save(joinpath(ROOT, "claude-answers", "img", "10-pol-large.png"), fig)
