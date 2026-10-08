# Stage 2: what the B-core initial condition looks like (reports/REPORT.md).
#
#   julia -t 24 --project=. code-claude/scripts/stage2_ic.jl
#
# Seeded states (the frame right after initialize_diff!) for a sharp core and for
# diluted edges, confined (2 layers, N0 1500) and unconfined (N0 500). Writes the
# B fraction of every seeded state and its radial profile to results/stage2/ic.csv,
# and a 3D gallery to results/stage2/ic_gallery.png.

using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new

const OUT = joinpath(@__DIR__, "..", "results", "stage2", "ic"); mkpath(OUT)
const WIDTHS = [0.0, 0.05, 0.1, 0.2]

"B fraction in rings of the normalised radius (in-plane if confined)."
function radial_profile(d, k; confined = true, edges = 0:0.1:1.2)
	x, y, z, s = d.x[k], d.y[k], d.z[k], d.cell_state[k]
	dz = confined ? 0.0 : 1.0
	rho = sqrt.((x .- mean(x)).^2 .+ (y .- mean(y)).^2 .+ dz .* (z .- mean(z)).^2)
	rho ./= quantile(rho, 0.98)
	[begin
		sel = findall(v -> edges[i] <= v < edges[i+1], rho)
		isempty(sel) ? NaN : count(==(2), s[sel]) / length(sel)
	 end for i in 1:length(edges)-1]
end

panels = NamedTuple[]
open(joinpath(OUT, "ic.csv"), "w") do io
	println(io, "case,width,fB_seeded,profile_0.0-0.1,...")
	for w in WIDTHS
		r = run_confined(2, 1500; wall = :rigid, width = w, tf = 1.0, outdir = OUT, verbose = false)
		d = load_frames(r.file)
		fB = count(==(2), d.cell_state[r.m1]) / length(d.cell_state[r.m1])
		@printf("confined 2 layers, width %.2f: fB = %.3f\n", w, fB)
		println(io, "confined_L2,", w, ",", round(fB, digits = 4), ",", join(round.(radial_profile(d, r.m1), digits = 3), ","))
		push!(panels, (d = d, k = r.m1, label = @sprintf("2 layers, width %.2f", w)))
	end
	for w in [0.0, 0.1]
		u = run_unconfined(500; width = w, tf = 1.0, outdir = OUT)
		d = u.frames
		fB = count(==(2), d.cell_state[u.m1]) / length(d.cell_state[u.m1])
		@printf("unconfined, width %.2f: fB = %.3f\n", w, fB)
		println(io, "unconfined,", w, ",", round(fB, digits = 4), ",", join(round.(radial_profile(d, u.m1; confined = false), digits = 3), ","))
		push!(panels, (d = d, k = u.m1, label = @sprintf("unconfined, width %.2f", w)))
	end
end

fig = plot_gallery(panels; ncol = length(panels), views = (:top, :side), size_ratio = 0.6,
				   title = "Seeded initial condition: B = 60 % of all cells (green B, orange C)")
save(joinpath(OUT, "ic_gallery.png"), fig)
println("wrote ", joinpath(OUT, "ic_gallery.png"))
