# Quick check that every adapted function runs (small sizes, short times).
#
#   julia -t auto --project=. code-claude/scripts/smoke_test.jl

using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new
out = joinpath(@__DIR__, "..", "results", "smoke"); mkpath(out)

fracB(d, k) = count(==(2), d.cell_state[k]) / length(d.cell_state[k])

r = run_confined(2, 400; tf = 60.0, wall = :hookean, outdir = out)
d = load_frames(r.file)
@printf("seeded fB = %.3f (want 0.600)   last fB = %.3f   frames %d\n", fracB(d, r.m1), fracB(d, length(d.x)), length(d.x))

r = run_confined(1, 300; tf = 60.0, wall = :rigid, width = 0.1, outdir = out)
d = load_frames(r.file)
@printf("rigid 1 layer: seeded fB = %.3f   z in [%.3f, %.3f] (band [1, 1.1])\n", fracB(d, r.m1), minimum(d.z[end]), maximum(d.z[end]))

u = run_unconfined(300; tf = 60.0, outdir = out)
@printf("unconfined: seeded fB = %.3f\n", fracB(u.frames, u.m1))

f = fusion_assay(PSETS["july"].mu, PSETS["july"].f_att, PSETS["july"].f_rep; n = 40, hours = 1.0, pairs_ = [(2,3)])

save(joinpath(out, "smoke_3d.png"), plot_pancake_frames(d, color_map, [r.m1, length(d.x)]))
println("smoke test done")

# --- Oct 7 additions: the diluted core (f_in), plate friction, the final-figure functions ---
r = run_confined(2, 300; tf = 60.0, f_in = 0.75, lambda_wall = 3.0, outdir = out, verbose = false)
d = load_frames(r.file)
x, y, s = d.x[r.m1], d.y[r.m1], d.cell_state[r.m1]
rho = hypot.(x .- mean(x), y .- mean(y)); core = sortperm(rho)[1:count(==(2), s)]
@printf("f_in 0.75 + friction 3: seeded fB = %.3f   B inside the core = %.3f (want ~0.75)\n",
		fracB(d, r.m1), count(==(2), s[core]) / count(==(2), s))

fd = joinpath(out, "figures")
write_panels(joinpath(fd, "data", "smoke_fig"), "smoke figure",
	[(panel = "p$j", row = 1, col = j, views = "top+side", label = "frame $k",
	  cells = (x = d.x[k], y = d.y[k], z = d.z[k], r = d.r[k], cell_state = d.cell_state[k]), meta = (frame = k,))
	 for (j, k) in enumerate((r.m1, length(d.x)))])
replot_figures(fd)
c = load_panel(joinpath(fd, "data", "smoke_fig", "p1.csv"))
@printf("panel round trip: %d cells, B fraction %.3f\n", length(c.x), count(==(2), c.cell_state) / length(c.x))
record_pancake_frames(d, color_map, r.m1:length(d.x); filename = joinpath(out, "smoke.mp4"))
println("smoke test (Oct 7 additions) done")

# --- Oct 8 additions: B fraction, C mixed into the core (c_core), lambda_scale, the final set, the
#     four phase-space versions ---
r = run_confined(2, 300; pset = "final", tf = 60.0, fB = 0.65, c_core = 0.1, lambda_wall = LAMBDA_WALL,
				 outdir = out, verbose = false)
d = load_frames(r.file)
x, y, s = d.x[r.m1], d.y[r.m1], d.cell_state[r.m1]
rho = hypot.(x .- mean(x), y .- mean(y)); nB = count(==(2), s)
core = sortperm(rho)[1:round(Int, nB / 0.9)]
@printf("final set, fB 0.65, c_core 0.1: seeded fB = %.3f (want 0.650)   C in the core = %.3f (want 0.10)   file %s\n",
		fracB(d, r.m1), count(==(3), s[core]) / length(core), basename(r.file))
u = run_unconfined(200; pset = "final", tf = 30.0, lambda_scale = 2.0, outdir = out, verbose = false)
@printf("unconfined, lambda x 2: seeded fB = %.3f   file %s\n", fracB(u.frames, u.m1), basename(u.file))
rows = [(layers = L, N0 = N, class = c, D = sqrt(N / L)) for (L, N, c) in
		((1, 100, "CENTRAL"), (1, 400, "ISLANDS"), (2, 100, "CENTRAL"), (2, 400, "RING"), (4, 100, "POLARISED"), (4, 400, "CENTRAL"))]
free = [(N0 = 100, class = "POLARISED"), (N0 = 400, class = "CENTRAL")]
for (x_, y_) in ((:log, :layers), (:linear, :layers), (:log, :diameter), (:linear, :diameter))
	save(joinpath(out, "smoke_phase_$(x_)_$(y_).png"), plot_phase_space(rows; free = free, x = x_, y = y_, title = "$x_ $y_"))
end
fig = Figure(size = (400, 400))
phase_space!(fig, rows; free = free, region_labels = Dict("PATTERN" => "patterns", "CENTRAL" => "central"),
			 pie_scale = 0.5, axis = (bbox = BBox(40, 380, 40, 380),))
save(joinpath(out, "smoke_phase_bbox.png"), fig)
println("smoke test (Oct 8 additions) done")
