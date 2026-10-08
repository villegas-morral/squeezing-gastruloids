# Stage 4: 3D galleries of the screened candidates (reports/REPORT.md).
#
#   julia --project=. code-claude/scripts/stage4_gallery.jl [fB_read] [cand1,cand2,...] [nooverview]
#
# Reads results/stage4/screen.csv and the saved trajectories; draws nothing new
# from simulation. For every candidate:
#   results/stage4/gallery_<cand>.png   top: the 5 screens at the readout
#                                        bottom: the 2-layer run through time
# and one overview, results/stage4/overview.png (top cutaways only).

using Statistics, Printf, DelimitedFiles, JLD2, CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))

const OUT   = joinpath(@__DIR__, "..", "results", "stage4")
const FREAD = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.2
raw, hdr = readdlm(joinpath(OUT, "screen.csv"), ','; header = true)
col = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
rows = [Dict(k => raw[r, i] for (k, i) in col) for r in 1:size(raw, 1)]
cands = length(ARGS) >= 2 ? split(ARGS[2], ",") : unique(String(r["cand"]) for r in rows)

"Frame index of a trajectory where fB first drops to `f` (after the seeding frame)."
function frame_at(d, f)
	fB = [count(==(2), s) / length(s) for s in d.cell_state]
	m1 = something(findlast(iszero, d.t), 1)          # initialize_diff! resets t to 0
	k = findfirst(i -> i >= m1 && fB[i] <= f, eachindex(fB))
	k === nothing ? length(fB) : k
end
"One-frame NamedTuple in the load_frames format, from a single configuration."
single(c, t) = (x = [c.x], y = [c.y], z = [c.z], r = [c.r], cell_state = [c.cell_state], t = [t], N = [length(c.x)])

overview = NamedTuple[]
for cand in cands
	sel = filter(r -> String(r["cand"]) == cand, rows)
	get_row(s) = (i = findfirst(r -> String(r["screen"]) == s, sel); i === nothing ? nothing : sel[i])
	panels = NamedTuple[]
	for (scr, ttl, confined) in [("small_L1", "1 layer, N0 200", true), ("islands_L1", "1 layer, N0 1500", true),
								 ("islands_L2", "2 layers, N0 2500", true), ("polar", "unconfined, N0 500", false)]
		r = get_row(scr)
		(r === nothing || String(r["file"]) == "-") && continue
		d = load_frames(joinpath(OUT, String(r["file"])))
		k = frame_at(d, FREAD)
		st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = confined)
		lbl = @sprintf("%s\n%s\nn=%d  S=%.2f  iso=%.2f", ttl, classify_pattern(st), st.n_isl, st.S, st.f_iso)
		push!(panels, (d = d, k = k, label = lbl))
		scr in ("islands_L1", "islands_L2", "small_L1", "polar") &&
			push!(overview, (d = d, k = k, label = @sprintf("%s\n%s", cand, ttl)))
	end
	# the fusion assay is stored once per base candidate ("name" for "name@w0.1")
	base = String(split(cand, "@")[1])
	fi = findfirst(r -> String(r["cand"]) == base && String(r["screen"]) == "fusion", rows)
	r = fi === nothing ? nothing : rows[fi]
	if r !== nothing && String(r["file"]) != "-"
		f = jldopen(joinpath(OUT, String(r["file"])), "r") do fh; fh["BC"]; end
		push!(panels, (d = single(f, 360.0), k = 1, label = @sprintf("BC doublet, 24 h\nmetric %.2f", Float64(r["bc24"]))))
	end
	# the 2-layer run through time
	r2 = get_row("islands_L2")
	if r2 !== nothing && String(r2["file"]) != "-"
		d = load_frames(joinpath(OUT, String(r2["file"])))
		for f in (0.5, 0.35, 0.2, 0.12)
			k = frame_at(d, f)
			push!(panels, (d = d, k = k, label = @sprintf("2 layers N0 2500\nt = %.0f (%.1f h), fB %.2f", d.t[k], d.t[k]/15,
							count(==(2), d.cell_state[k]) / length(d.cell_state[k]))))
		end
	end
	isempty(panels) && continue
	fig = plot_gallery(panels; ncol = 5, views = (:top, :side), size_ratio = 0.55, labelsize = 17,
					   title = "candidate $cand   (readout at fB = $FREAD)")
	save(joinpath(OUT, "gallery_$(cand).png"), fig)
	println("wrote gallery_$(cand).png")
end
if !isempty(overview) && !(length(ARGS) >= 3 && ARGS[3] == "nooverview")
	fig = plot_gallery(overview; ncol = 4, views = (:top,), size_ratio = 0.5, labelsize = 15,
					   title = "all candidates at fB = $FREAD (top cutaways)")
	save(joinpath(OUT, "overview.png"), fig)
	println("wrote overview.png")
end
