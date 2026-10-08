# Stage 2 figures: how each plate model holds the slab thickness.
#
#   julia --project=. code-claude/scripts/stage2_wall_figs.jl
#
# Reads results/stage2/walls.csv and the saved runs (no simulation). Writes
#   results/stage2/walls_thickness.png   effective layers vs time, per wall model
#   results/stage2/walls_side3d.png      3D side views at the end, 2 and 4 layers

using Statistics, Printf, DelimitedFiles, JLD2, CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))

const OUT = joinpath(@__DIR__, "..", "results", "stage2")
raw, hdr = readdlm(joinpath(OUT, "walls.csv"), ','; header = true)
col = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
rows = [(wall = String(raw[r, col["wall"]]), layers = Int(raw[r, col["layers"]]),
		 N0 = Int(raw[r, col["N0"]]), secs = Float64(raw[r, col["secs"]]),
		 file = String(raw[r, col["file"]])) for r in 1:size(raw, 1) if raw[r, col["ok"]] == true]

walls  = ["constant", "hookean", "hookean_stiff", "rigid"]
names  = Dict("constant" => "constant force (original)", "hookean" => "Hookean (July)",
			  "hookean_stiff" => "Hookean, stiff (dt 0.005)", "rigid" => "rigid")
colors = Dict("constant" => :gray40, "hookean" => :steelblue, "hookean_stiff" => :darkorange, "rigid" => :seagreen)

fig = Figure(size = (1500, 430))
for (j, L) in enumerate((1, 2, 4))
	ax = Axis(fig[1, j], title = "$L layer$(L > 1 ? "s" : "") (gap $(round(2.1L, digits = 1)))",
			  xlabel = "time (h)", ylabel = j == 1 ? "effective layers" : "",
			  yticks = 0:0.5:5)
	hlines!(ax, [L]; color = :black, linestyle = :dash)
	for w in walls
		i = findfirst(r -> r.wall == w && r.layers == L, rows)
		i === nothing && continue
		d = load_frames(joinpath(OUT, rows[i].file))
		m1 = something(findlast(iszero, d.t), 1)
		th = thickness_timeline(d, 2.1L; from = m1)
		lines!(ax, [x.t / 15 for x in th], [x.layers_eff for x in th]; color = colors[w], linewidth = 3,
			   label = @sprintf("%s (%.0f s)", names[w], rows[i].secs))
	end
	axislegend(ax; position = :lt, labelsize = 11)
end
save(joinpath(OUT, "walls_thickness.png"), fig)
println("wrote walls_thickness.png")

panels = NamedTuple[]
for L in (2, 4), w in walls
	i = findfirst(r -> r.wall == w && r.layers == L, rows)
	i === nothing && continue
	d = load_frames(joinpath(OUT, rows[i].file))
	k = length(d.x)
	th = thickness_timeline(d, 2.1L; from = k)[1]
	push!(panels, (d = d, k = k, label = @sprintf("%s, %d layers\n%.2f effective layers", names[w], L, th.layers_eff)))
end
fig = plot_gallery(panels; ncol = 4, views = (:side,), size_ratio = 0.6, labelsize = 15,
				   title = "end of the run (t = 325), side views: green B, orange C")
save(joinpath(OUT, "walls_side3d.png"), fig)
println("wrote walls_side3d.png")
