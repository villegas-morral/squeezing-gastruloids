# Cell diffusion in a monolayer, for the Damkohler number (November 2026).
#
#   julia --project=. code-claude/scripts/msd_monolayer.jl run <N0> <seed>
#   julia --project=. code-claude/scripts/msd_monolayer.jl analyse
#
# run      one monolayer (1 layer, the final set "final", rigid plates, plate friction
#          LAMBDA_WALL) from seeding to fB = 0.15, one frame per time unit (4 min) instead
#          of one per hour, with the cell ids (save_frames) -> results/msd/
# analyse  msd_2d of every run in results/msd/: the in-plane MSD of B and C cells,
#          followed by id between divisions, with the growth flow removed; D from a
#          straight line at long lag times, and the Damkohler numbers
#          -> results/msd/msd.csv, msd_fit.csv, msd.png
#
# There is no noise in the model: cells move because neighbours divide and push, and
# because differentiation changes the adhesion. So D is measured in growing
# monolayers, not in a resting one (that would give D = 0).

using Random, Statistics, Printf, DelimitedFiles, LinearAlgebra, JLD2, CairoMakie

const MODE = ARGS[1]
const OUT  = joinpath(@__DIR__, "..", "results", "msd")
mkpath(OUT)

if MODE == "run"
	using Distributions, ProgressMeter, CellBasedModels
	import CellBasedModels: update!
	include(joinpath(@__DIR__, "gastruloids_lib.jl"))
	model = model_new
	N0, seed = parse(Int, ARGS[2]), parse(Int, ARGS[3])
	r = run_confined(1, N0; pset = "final", wall = :rigid, lambda_wall = LAMBDA_WALL, tf = 900.0, fB_stop = 0.15,
					 seed = seed, outdir = OUT, save_every = 1.0, verbose = false)
	println("done ", basename(r.file), "  (", round(Int, r.secs), " s)")
else
	include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
	include(joinpath(@__DIR__, "stage5_common.jl"))      # frame_at
	const Q = 0.004667                  # B -> C rate (per time unit)
	const UM2_PER_H = 25.0 * 15.0       # r^2 per time unit -> um^2 per hour (r = 5 um, 15 units = 1 h)
	files = sort(filter(f -> startswith(f, "com_") && endswith(f, ".jld2"), readdir(OUT)))
	fits = NamedTuple[]; curves = Dict{Tuple{String,Symbol},Any}()
	for f in files
		m = match(r"_N(\d+)_.*_s(\d+)\.jld2", f); N0, seed = parse(Int, m[1]), parse(Int, m[2])
		d = load_frames(joinpath(OUT, f))
		m1 = frame_at(d, 0.6)                                   # seeding frame (t = 0), old ids
		from = findfirst(k -> k > m1 && d.t[k] >= 26.0, eachindex(d.t))   # after the 25-unit relaxation: cells divide
		x, y = d.x[m1], d.y[m1]
		R0 = sqrt(2mean((x .- mean(x)).^2 .+ (y .- mean(y)).^2)) + 1.0    # pancake radius at squeezing (units of r)
		for drift in (:affine, :local)
			c = msd_2d(d; from = from, drift = drift)
			curves[(f, drift)] = c
			for (typ, ms) in (("B", c.msd_B), ("C", c.msd_C), ("all", c.msd))
				ft = fit_diffusion(c.tau, ms; tmin = 30, tmax = 120)
				push!(fits, (file = f, N0 = N0, seed = seed, drift = drift, type = typ, D = ft.D, slope = ft.slope,
							 D_um2h = ft.D * UM2_PER_H, R0 = R0, Da_cell = Q * 4 / ft.D, Da_R = Q * R0^2 / ft.D,
							 hours = (d.t[end] - d.t[from]) / 15))
			end
		end
		@printf("%s: window %.1f h from t = %.0f\n", f, (d.t[end] - d.t[from]) / 15, d.t[from])
	end
	open(joinpath(OUT, "msd_fit.csv"), "w") do io
		println(io, "file,N0,seed,drift,type,D_r2_per_unit,D_um2_per_h,slope_loglog,R0_r,Da_cell,Da_R,window_h")
		for r in fits
			@printf(io, "%s,%d,%d,%s,%s,%.5f,%.2f,%.3f,%.2f,%.4f,%.3f,%.1f\n", r.file, r.N0, r.seed, r.drift, r.type,
					r.D, r.D_um2h, r.slope, r.R0, r.Da_cell, r.Da_R, r.hours)
		end
	end
	open(joinpath(OUT, "msd.csv"), "w") do io
		println(io, "file,drift,tau,msd_B,msd_C,msd_all,n_B,n_C")
		for ((f, drift), c) in sort(collect(curves); by = first), i in eachindex(c.tau)
			@printf(io, "%s,%s,%.3f,%.5f,%.5f,%.5f,%d,%d\n", f, drift, c.tau[i], c.msd_B[i], c.msd_C[i], c.msd[i], c.n_B[i], c.n_C[i])
		end
	end
	for r in fits
		r.type == "all" || @printf("N0 %4d s%d %-7s %s: D = %.4f r^2/unit = %5.1f um^2/h  slope %.2f  Da_cell %.3f  Da_R %.2f\n",
								   r.N0, r.seed, r.drift, r.type, r.D, r.D_um2h, r.slope, r.Da_cell, r.Da_R)
	end
	# the figure: MSD against lag time (log-log), B and C, both sizes, affine-corrected; the fits
	fig = Figure(size = (1100, 480), fontsize = 15, fonts = TEX_FONTS)
	for (j, drift) in enumerate((:affine, :local))
		ax = Axis(fig[1, j]; xscale = log10, yscale = log10, xlabel = L"\mathrm{lag\ time}\ \tau\ \mathrm{(h)}",
				  ylabel = L"\langle |\Delta \mathbf{r}|^2 \rangle\ (r^2)",
				  title = drift == :affine ? "growth flow removed (affine fit)" : L"\mathrm{relative\ to\ the\ neighbours}\ (5r)")
		for f in files, (typ, col, key) in (("B", Makie.wong_colors()[3], :msd_B), ("C", Makie.wong_colors()[2], :msd_C))
			c = curves[(f, drift)]; N0 = parse(Int, match(r"_N(\d+)_", f)[1])
			ok = findall(>(0), getfield(c, key))
			lines!(ax, c.tau[ok] ./ 15, getfield(c, key)[ok]; color = col, linestyle = N0 == 600 ? :dash : :solid,
				   label = Makie.LaTeXString("\$$(typ)\$, \$N_0 = $(N0)\$"))
		end
		t = [2.0, 8.0]; lines!(ax, t, 0.3 .* t ./ 2; color = :black, linewidth = 1)   # slope 1 guide
		text!(ax, 8.0, 1.2; text = "slope 1", fontsize = 12)
		vspan!(ax, 2.0, 8.0; color = (:grey, 0.12))
		j == 2 && axislegend(ax; position = :lt, unique = true, framevisible = false)
	end
	save(joinpath(OUT, "msd.png"), fig)
	println("wrote ", joinpath(OUT, "msd_fit.csv"), ", msd.csv, msd.png")
end
