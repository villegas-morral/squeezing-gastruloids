# Stage 5: the investigation figures of claude-answers/03 (3D, from the tuning runs).
#
#   julia --project=. code-claude/scripts/stage5_answer_figs.jl
#
# Reads results/stage5/tune/ (and phase_A/ for one run); writes claude-answers/img/03-*.png:
#   03-core-and-friction.png  2 layers, N0 150: sharp vs diluted core, and plate friction
#   03-looseness.png          1 layer N0 1500 and 2 layers N0 2400 at fB = 0.2, per set
#   03-friction-monolayers.png  big monolayers without and with plate friction
# (The final figures, with their data, are made by stage6_figures.jl.)

using Statistics, Printf, DelimitedFiles, JLD2, CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))      # frame_at

const TUNE = joinpath(@__DIR__, "..", "results", "stage5", "tune")
const IMG  = joinpath(@__DIR__, "..", "..", "claude-answers", "img")

"The saved run of a set at a point (tuning runs, then the phase spaces), or nothing."
function run_file(cand, core, L, N0, seed)
	base, _, lw = parse_cand_name(cand)
	f = @sprintf("com_L%g_N%d_%s_rigid_%s%s_s%d.jld2", L, N0, cand, core == "sharp" ? "w0" : core,
				 lw > 0 ? @sprintf("_lw%g", lw) : "", seed)
	for sub in ("tune", "phase_B", "phase_A")
		p = joinpath(TUNE, "..", sub, f)
		isfile(p) && return p
	end
	nothing
end
parse_cand_name(c) = (p = split(c, "@"); (p[1], 0.0, length(p) > 1 ? parse(Float64, p[2][3:end]) : 0.0))

function panel(f, which, label)
	d = load_frames(f); k = frame_at(d, which)
	st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k])
	(d = d, k = k, label = which >= 0.6 ? label : label * @sprintf("\n%s  offset %.2f  S %.2f", lowercase(classify_pattern(st)), st.offset, st.S))
end

# ---- sharp vs diluted core, and friction: 2 layers, N0 150 --------------------------------
pans = NamedTuple[]
for (cand, core, name) in (("july_frB4.5_frC6", "sharp", "sharp core"), ("july_frB4.5_frC6", "fin0.75", "diluted core"),
						   ("july_frB4.5_frC6@lw1.5", "sharp", "sharp core + friction 1.5"))
	f = run_file(cand, core, 2, 150, 1)
	f === nothing || push!(pans, panel(f, 0.6, "$name\nstart (seed 1)"))
	for s in 1:3
		f = run_file(cand, core, 2, 150, s)
		f === nothing || push!(pans, panel(f, 0.2, "$name\nfB = 0.2, seed $s"))
	end
end
fig = plot_gallery(pans; ncol = 4, views = (:top, :side), size_ratio = 0.5, labelsize = 13,
				   title = "2 layers, N0 150, july_frB4.5_frC6: where B is at the readout")
save(joinpath(IMG, "03-core-and-friction.png"), fig); println("wrote 03-core-and-friction.png")

# ---- looseness of the patterns ------------------------------------------------------------
pans = NamedTuple[]
sets = (("july_frB4.5", "july_frB4.5 (your pick)"), ("july_frB4.5_frC6", "july_frB4.5_frC6"),
		("july_frB4.5_frC6@lw1.5", "july_frB4.5_frC6 + friction 1.5 (chosen)"), ("fusion", "FUSION"))
for (L, N0) in ((1, 1500), (2, 2400)), (cand, name) in sets
	f = run_file(cand, "sharp", L, N0, 1)
	f === nothing && (f = run_file(cand, "fin0.75", L, N0, 1))
	f === nothing || push!(pans, panel(f, 0.2, "$name\n$L layer$(L > 1 ? "s" : ""), N0 $N0, fB = 0.2"))
end
fig = plot_gallery(pans; ncol = 4, views = (:top,), size_ratio = 0.6, labelsize = 13,
				   title = "How loose the patterns are (seed 1, top cutaway; S = B fraction around a B cell, lower = looser)")
save(joinpath(IMG, "03-looseness.png"), fig); println("wrote 03-looseness.png")

# ---- plate friction and big monolayers ----------------------------------------------------
pans = NamedTuple[]
for N0 in (4800, 2400), (cand, name) in (("july_frB4.5_frC6", "no friction"), ("july_frB4.5_frC6@lw1.5", "friction 1.5"),
										 ("july_frB4.5_frC6@lw3", "friction 3"))
	f = run_file(cand, "sharp", 1, N0, 1)
	f === nothing || push!(pans, panel(f, 0.2, "1 layer, N0 $N0, $name"))
end
fig = plot_gallery(pans; ncol = 3, views = (:top,), size_ratio = 0.9, labelsize = 15)
save(joinpath(IMG, "03-friction-monolayers.png"), fig); println("wrote 03-friction-monolayers.png")
