# =============================================================================
# gastruloids_lib.jl -- model, parameters and analysis for the squeezing project
#
# EXPORTED FROM squeezing-gastruloids-claude.ipynb. Do not edit by hand; edit the
# notebook and re-export with scripts/export_lib.ps1 so the two cannot drift.
#
# Usage:
#   using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
#   using CellBasedModels; import CellBasedModels: update!
#   using JLD2, CairoMakie
#   include("gastruloids_lib.jl")
#   model = model_new
#
# `differentiate!` calls IJulia.clear_output only when show_progress = true, so
# batch scripts must pass show_progress = false (run_confined already does).
# =============================================================================
# ---- notebook cell 51154174 ----------------------------------------------------
color_map = Dict(
    1 => Makie.wong_colors()[1],
    2 => Makie.wong_colors()[3],
    3 => Makie.wong_colors()[2]
)
color_map_alpha = Dict(
    1 => Makie.wong_colors()[1],
    2 => Makie.wong_colors()[3],
    3 => RGBAf(Makie.wong_colors()[2], 0.2f0)
);

# ---- notebook cell c17fa1b1 ----------------------------------------------------
function initialize_growth(parameters; dt, simBox = [-60.0 60; -60 60; -60 60])

	com = Community(
		model,
		N = 1,
		dt = dt,
		simBox = simBox,
	)

	# Prameters
	for (par, val) in pairs(parameters)
		com[par] = val
	end

	# Initialization
	# The aggregate grows as C cells; initialize_diff! seeds the B core later.
	# (The original grew type 1 = A, which no longer exists. The growth mechanics,
	# parameters_tfm, are the same for every type, so this only matters for the
	# relaxation right after switching to parameters_do: it now runs under CC
	# mechanics, the tissue that surrounds the core.)
	com.cell_state = 3
	com.x = 0.0
	com.y = 0.0
	com.z = 0.0
	com.vx = 0.0
	com.vy = 0.0
	com.vz = 0.0
	com.t_div = 1

	com.g_on = true
	com.d_on = false

	return com

end;

# ---- notebook cell f2f50d7d ----------------------------------------------------
function grow_size!(com, save_each, n_cells;
	n_control=100, show_progress = true)

	com.abm.neighbors.cellEdge = (2 * maximum(com.r)) * maximum(com.mu)

	loadToPlatform!(com, preallocateAgents = round(Int, n_cells * 1.1))
	i = 0
	n_total = 0

	while (com.N .< n_cells)
		if show_progress && com.N > n_total
			println("N > $(n_total)")
			n_total += n_control
		end
		i += 1
		agentStepDE!(com)
		agentStepRule!(com)
		update!(com)
		computeNeighbors!(com)
		if i % save_each == 0
			saveRAM!(com)
		end
	end
	for i in 1:save_each
		agentStepDE!(com)
		update!(com)
		computeNeighbors!(com)
	end
	saveRAM!(com)
	bringFromPlatform!(com)
end;

# ---- notebook cell f818a62d ----------------------------------------------------
"""
	initialize_diff!(com; g_on = true, t_reset = true, fB = 0.6, width = 0.0, f_in = nothing,
	                 confined = true, c_core = 0.0)

Seed the B core and switch differentiation on.

* `fB`        fraction of ALL cells that start as B (meeting: about 60 %). Exactly
              round(fB * N) cells become B: the innermost ones.
* `width`     0 = sharp core: the innermost cells are B and all the others C.
              width > 0 dilutes the edge: a cell at distance rho from the centre
              (in units of the aggregate radius) is B with probability
              1 / (1 + exp((rho - rho0) / width)), with rho0 such that exactly
              round(fB * N) cells are still B. 0.05 is a soft edge, 0.2 a strong
              dilution.
* `f_in`      another way to ask for the same dilution: the fraction of the B
              cells that stay inside the sharp core (the round(fB * N) innermost
              cells). 1 = sharp; 0.75 = a quarter of the B cells outside it, the
              rest of the core C (meeting, Oct 2026). The width that gives it is
              found by bisection; `width` is then ignored.
* `confined`  true: rho is the IN-PLANE distance, so the core is a cylinder
              through the whole slab (squeezing a sphere that has a core gives a
              core spanning the thickness). false: the 3D distance, a spherical core.
* `c_core`    fraction of the core left C, mixed in at random (the July seeding: 75 %
              B + 15 % A + 10 % C shuffled inside the core). The core is then the
              round(fB * N / (1 - c_core)) innermost cells, so exactly round(fB * N)
              are still B. 0 = a pure B core.

The aggregate radius is the 98th percentile of rho, so a few stray cells cannot
inflate it. Returns the width used.
"""
function initialize_diff!(com;
    g_on = true,
    t_reset = true,
    fB = 0.6,
    width = 0.0,
    f_in = nothing,
    confined = true,
    c_core = 0.0
    )

    # distance of every cell from the aggregate centre, in units of its radius
    xc, yc, zc = mean(com.x), mean(com.y), mean(com.z)
    if confined
        rho = sqrt.((com.x .- xc).^2 .+ (com.y .- yc).^2)
    else
        rho = sqrt.((com.x .- xc).^2 .+ (com.y .- yc).^2 .+ (com.z .- zc).^2)
    end
    rho = rho ./ quantile(rho, 0.98)
    n_B = round(Int, fB * length(rho))

    # The n_B cells with the smallest score become B. With width > 0 the score is
    # the radius blurred by logistic noise: that turns the sharp cut into the
    # smooth profile above while keeping exactly n_B cells B.
    score = copy(rho)
    if (f_in !== nothing && f_in < 1) || width > 0
        u = rand(length(rho))
        noise = log.(u ./ (1 .- u))
        if f_in !== nothing && f_in < 1
            core = falses(length(rho))
            core[sortperm(rho)[1:n_B]] .= true
            inside(w) = count(core[sortperm(rho .+ w .* noise)[1:n_B]]) / n_B
            lo, hi = 0.0, 0.1
            while inside(hi) > f_in; hi *= 2; end
            for _ in 1:30
                mid = (lo + hi) / 2
                inside(mid) > f_in ? (lo = mid) : (hi = mid)
            end
            width = hi
        end
        score .+= width .* noise
    end
    com.cell_state .= 3
    com.cell_state[sortperm(score)[1:n_B]] .= 2

    # C cells mixed into the core (salt and pepper inside it), as seeded in July
    if c_core > 0
        core = sortperm(score)[1:min(length(score), round(Int, n_B / (1 - c_core)))]
        com.cell_state .= 3
        com.cell_state[shuffle(core)[1:n_B]] .= 2
    end

    dt = com.dt
	com.g_on = g_on
    com.d_on = true

    if t_reset
        setfield!(com,:t, 0);
    end

    if g_on
        for i in 1:com.N
            sample_right = 2 * com.tau_div[com.cell_state[i]] * com.sigma_div[] + dt
            com.t_div[i] = com.t + CBMDistributions.uniform(dt, sample_right)
        end
    end

	saveRAM!(com)

    return width
end;

# ---- notebook cell 9154a7b7 ----------------------------------------------------
"Fraction of B cells (cell_state == 2) among the current cells. Works while loaded on the platform."
function fraction_B(com)
	n = com.N[1]
	cs = com.cell_state
	nB = 0
	@inbounds for i in 1:n
		cs[i] == 2 && (nB += 1)
	end
	return nB / n
end;

"""
	differentiate!(com, save_each, tf; fB_stop = nothing, show_progress = true, ...)

Proliferate and differentiate (B -> C) for a time `tf`, saving a frame every
`save_each` steps. By default the whole `tf` runs. With `fB_stop` set, the run
stops at the first saved frame where the fraction of B cells is <= fB_stop
(handy in sweeps; 0.2 is the experimental 96 h reference).
"""
function differentiate!(com, save_each, tf;
    preallocate = 2 * com.N * 2 ^ (tf / minimum(com.tau_div)),
    fB_stop = nothing,
    show_progress = true      # false in batch runs: skips the inline plot
    )

	com.abm.neighbors.cellEdge = (2 * maximum(com.r)) * maximum(com.mu)

    message1 = "Possible numerical instabilities.
            You might have to send argument preallocate = highnumber
            Default is preallocate = 2 * com.N * 2 ^ (tf / minimum(com.tau_div))"
    message2 = "Numerical instabilities.
            You might have to decrease the timestep"

	steps = round(Int64, tf / com.dt)
    step_control = steps / 10
    progress = step_control
    i_end = steps

    if com.g_on[]
        show_progress && println("Initial N: $(com.N)")
        loadToPlatform!(com, preallocateAgents = round(Int, preallocate))
    else
        loadToPlatform!(com)
    end

    for i in 1:steps
		if i > progress
            if show_progress
                IJulia.clear_output(true)
			    println("$(round(100 * i / steps))% \t N=$(com.N) \t fB=$(round(fraction_B(com), digits = 3))")
                plot_pancake(com, color_map_alpha, length(com), length(com); savefig=false, n_views=1, n=1)
            end
			progress += step_control
		end
        agentStepDE!(com)
        agentStepRule!(com)
        update!(com)
        computeNeighbors!(com)
        if i % save_each == 0
            saveRAM!(com)
            if fB_stop !== nothing && fraction_B(com) <= fB_stop
                i_end = i
                break
            end
        end
        if all(com.N .> preallocate)
            println(message1)
            break
        end
    end

	if i_end % save_each != 0
		saveRAM!(com)
    end

	bringFromPlatform!(com)

    if show_progress
        IJulia.clear_output(true)
        println("Final N: $(com.N) \t fB: $(round(fraction_B(com), digits = 3))")
    end

end;

# ---- notebook cell f06d6c5e ----------------------------------------------------
function stabilize!(com, save_each, tf; show_progress = true)

	com.abm.neighbors.cellEdge = (2 * maximum(com.r)) * maximum(com.mu)
	
	steps = round(Int64, tf / com.dt)
	i = 0

	prog_bar = show_progress ? Progress(steps, desc="Stabilizing") : nothing

	loadToPlatform!(com)

	while (i<steps)
		i += 1
		agentStepDE!(com)
		update!(com)
		computeNeighbors!(com)
		if i % save_each == 0
			saveRAM!(com)
		end
		show_progress && next!(prog_bar)
	end
	if i % save_each != 0
		saveRAM!(com)
	end
	bringFromPlatform!(com)

end;

# ---- notebook cell 462afe86 ----------------------------------------------------
function plot_pancake(com, color_map, mstart, mstop;
    boxsize = maximum([
        maximum(com.x) - minimum(com.x),
        maximum(com[mstart].x) - minimum(com[mstart].x)
    ]) / 1.3,
    n = 4,
    showtime = false,
    shownumbers = true,
    filename = "pancake",
    n_views = 3,
    savefig = false
)
    filename = "$(filename)_$(layers)layers_n$(com[mstart].N)"
    views = [
        (perspective = pi/2,  show_z = false, show_x = true, name = "top"),
        (perspective = 2*pi+pi/32,  show_z = true,  show_x = false, name = "side"),
        (perspective = -pi/2,  show_z = false, show_x = true, name = "bottom"),
    ]
    views = views[1:n_views]
    
    d = getParameter(com, [:x, :y, :z, :r, :cell_state])

    for (j, view) in enumerate(views)
        fig = Figure(size = (n * 640, 600), figure_padding = 40)
        labelsize = 40
        
        for (i, pos) in enumerate(range(start = mstart, length = n, stop = mstop))
            pos = floor(Int, pos)
            if j == 1
                println("Column $i: timestamp $pos")
            end
            t = round(com[pos].t, digits = 2)

            ax = Axis3(
                fig[1, i],
                aspect = :data,
                azimuth = 0,
                elevation = view.perspective,
                xlabel = "",
                ylabel = "",
                zlabel = "",
                xticklabelsize = labelsize,
                yticklabelsize = labelsize,
                zticklabelsize = labelsize,
                titlevisible = showtime,
                titlealign = :center,
                titlegap = 12,
                titlesize = labelsize,
                title = L"t=%$(t)"
            )

            color = [color_map[j] for j in d[:cell_state][pos]]
            meshscatter!(
                ax,
                d[:x][pos],
                d[:y][pos],
                d[:z][pos],
                markersize = d[:r][pos],
                color = color
            )

            xlims!(ax, -boxsize, boxsize)
            ylims!(ax, -boxsize, boxsize)
            zlims!(ax, -boxsize, boxsize)

            if !shownumbers
                ax.xticklabelsize = 0
                ax.yticklabelsize = 0
                ax.zticklabelsize = 0
            end

            ax.zticklabelsvisible = view.show_z
            ax.xticklabelsvisible = view.show_x
            ax.yticklabelsvisible = false
        end
        display(fig)
        if savefig
            save("$(filename)_$(view.name).png", fig)
        end
    end
end;

# ---- notebook cell ba36859e ----------------------------------------------------
function record_pancake(com, color_map, mstart, mstop, layers;
    boxsize = maximum([
        maximum(com.x) - minimum(com.x),
        maximum(com[mstart].x) - minimum(com[mstart].x)
    ]) / 1.5,
    n = (mstop-mstart),
    fps = 30,
    shownumbers = true,
    filename = "pancake",
    n_views = 3,
    size_ratio=1
)
    filename = "$(filename)_$(layers)layers_n$(com[mstart].N)_$(fps)fps.mp4"

    views = [
        (perspective = pi/2,  show_z = false, show_x = true,  name = "top"),
        (perspective = 2*pi+pi/32,  show_z = true,  show_x = false, name = "side"),
        (perspective = -pi/2,  show_z = false, show_x = true, name = "bottom"),
    ]
    views = views[1:n_views]

    d = getParameter(com, [:x, :y, :z, :r, :cell_state])
    frames = floor.(Int, range(start = mstart, length = n, stop = mstop))

    labelsize = 40*size_ratio
    fig = Figure(size = (n_views * 640*size_ratio, 600*size_ratio), figure_padding = 40*size_ratio)

    # Build all axes and scatter plots once
    axes = []
    scatter_plots = []
    pos0 = frames[1]

    for (i, view) in enumerate(views)
        ax = Axis3(
            fig[1, i],
            aspect = :data,
            azimuth = 0,
            elevation = view.perspective,
            xlabel = "", ylabel = "", zlabel = "",
            xticklabelsize = labelsize,
            yticklabelsize = labelsize,
            zticklabelsize = labelsize,
            titlevisible = true,
            titlealign = :center,
            titlegap = 12,
            titlesize = labelsize,
            title = view.name,
        )

        if !shownumbers
            ax.xticklabelsize = 0
            ax.yticklabelsize = 0
            ax.zticklabelsize = 0
        end
        ax.zticklabelsvisible = view.show_z
        ax.xticklabelsvisible = view.show_x
        ax.yticklabelsvisible = false

        xlims!(ax, -boxsize, boxsize)
        ylims!(ax, -boxsize, boxsize)
        zlims!(ax, -boxsize, boxsize)

        color0 = [color_map[cs] for cs in d[:cell_state][pos0]]
        sp = meshscatter!(
            ax,
            d[:x][pos0], d[:y][pos0], d[:z][pos0],
            markersize = d[:r][pos0],
            color = color0,
        )

        push!(axes, ax)
        push!(scatter_plots, sp)
    end

    println("Recording $filename")
    record(fig, filename, frames; framerate = fps) do pos
        # t = round(com[pos].t, digits = 2)
        # println("  frame t=$t")

        color = [color_map[cs] for cs in d[:cell_state][pos]]

        for (ax, sp) in zip(axes, scatter_plots)
            sp[1][] = Point3f.(d[:x][pos], d[:y][pos], d[:z][pos])
            sp.markersize[] = d[:r][pos]
            sp.color[] = color
            # ax.title = L"t=%$(t)"
        end
    end

    println("Saved $filename")
end;

# ---- notebook cell 4cd7bbef ----------------------------------------------------
# --- saving and reloading runs ------------------------------------------------
#
# `jldsave("com.jld2"; com)` does NOT round-trip. ABM() eval's the generated step
# functions as ANONYMOUS closures in Main, and the pickled Community stores a
# reference to them by their gensym name. Reloading needs a session whose gensym
# counter produced the same names, which a notebook will not reproduce once cells
# are re-run in a different order or the kernel restarts, so you get
#     Warning: type Main.#19#20 does not exist in workspace; reconstructing
#     MethodError: Cannot `convert` an object of type
#                  JLD2.ReconstructedSingleton{Symbol("#19#20")} to ... Function
# Plain arrays carry no such dependency, load in any session with the model not
# even defined, and are ~7x smaller.

"""
Save the trajectory of `com` as plain arrays (session-independent). Also saves
each cell's `id` (CellBasedModels numbers every cell; the two daughters of a
division get new ids), so cells can be followed between divisions (msd_2d).
"""
function save_frames(file, com)
	d  = getParameter(com, [:x, :y, :z, :r, :cell_state])
	nf = length(com)
	isempty(dirname(file)) || mkpath(dirname(file))
	jldsave(file;
		x = d[:x], y = d[:y], z = d[:z], r = d[:r], cell_state = d[:cell_state],
		id = [Vector{Int}(v) for v in getParameter(com, :id)],
		t = [com[k].t for k in 1:nf], N = [com[k].N for k in 1:nf])
	return file
end;

"""
Load what `save_frames` wrote. Returns a NamedTuple of per-frame vectors; `id` is
`nothing` for files saved before November 2026 (no ids).
"""
load_frames(file) = jldopen(file, "r") do fh
	(x = fh["x"], y = fh["y"], z = fh["z"], r = fh["r"],
	 cell_state = fh["cell_state"], t = fh["t"], N = fh["N"],
	 id = haskey(fh, "id") ? fh["id"] : nothing)
end;

"""
	plot_frames(d, color_map, ks; slice=true, halfwidth=1.2, ...)

Flat 2D version, one column per frame index in `ks`.

`slice=true` keeps only cells within `halfwidth` of the mid-plane. Use it: a top
projection through the whole slab lets the C roof hide the islands underneath.
"""
function plot_frames(d, color_map, ks;
		slice = true, halfwidth = 1.2, markersize = 4,
		filename = "frames", savefig = false)
	fig = Figure(size = (320*length(ks), 340), fonts = TEX_FONTS)
	for (i, k) in enumerate(ks)
		x, y, z, s = d.x[k], d.y[k], d.z[k], d.cell_state[k]
		keep = slice ? findall(v -> abs(v - mean(z)) < halfwidth, z) : collect(eachindex(z))
		fAB  = count(v -> v == 1 || v == 2, s) / length(s)
		ax = Axis(fig[1, i], aspect = DataAspect(), titlesize = 13,
			title = "t=$(round(d.t[k], digits=0))  N=$(length(s))  fAB=$(round(fAB, digits=2))")
		hidedecorations!(ax); hidespines!(ax)
		o = sort(keep; by = j -> s[j], rev = true)      # draw C first
		scatter!(ax, x[o] .- mean(x), y[o] .- mean(y);
			color = [color_map[s[j]] for j in o], markersize = markersize)
	end
	savefig && save("$(filename).png", fig)
	fig
end;

"""
	plot_pancake_frames(d, color_map, ks; kw...)

`plot_pancake` for the array format from `load_frames`: the same Axis3 +
meshscatter rendering, one column per frame index in `ks`, one row per view.

Kept separate from `plot_pancake` because that one needs a live `Community`
(`com[k].t`, `com[k].N`), which the session-independent format deliberately does
not carry.

    ghost = nothing   opaque cells (default). ghost = GHOST: every view shows the
                      whole aggregate, C translucent, so the B cells show through
                      the C cells around them (unconfined aggregates; needs GLMakie)
    cut_top = true    with ghost = nothing: drop the cells above the mid-plane in
                      the TOP view instead. It is a cutaway, not a slice -- the
                      side view is left intact, so the constant slab thickness
                      is still checkable by eye.
"""
function plot_pancake_frames(d, color_map, ks;
		n_views = 2, showtime = true, shownumbers = false, cut_top = true,
		boxsize = nothing, filename = "pancake_frames", savefig = false,
		labels = nothing, size_ratio = 1.0, ghost = nothing)

	views = [
		(perspective = pi/2,        show_z = false, show_x = true,  name = "top"),
		(perspective = 2*pi+pi/32,  show_z = true,  show_x = false, name = "side"),
		(perspective = -pi/2,       show_z = false, show_x = true,  name = "bottom"),
	][1:n_views]

	bs = boxsize === nothing ?
		maximum(maximum(d.x[k]) - minimum(d.x[k]) for k in ks) / 1.9 : boxsize
	labelsize = 26*size_ratio
	fig = Figure(fonts = TEX_FONTS, size = (round(Int, 520*size_ratio*length(ks)),
						 round(Int, 470*size_ratio*n_views)), figure_padding = 20)

	for (j, view) in enumerate(views), (i, k) in enumerate(ks)
		x, y, z, r, s = d.x[k], d.y[k], d.z[k], d.r[k], d.cell_state[k]
		zc = mean(z)
		gh = ghost
		keep = (view.name == "top" && cut_top && gh === nothing) ? findall(v -> v <= zc + 0.2, z) :
			   collect(eachindex(z))
		fAB = count(v -> v == 1 || v == 2, s) / length(s)
		ttl = j == 1 ?
			(labels === nothing ?
				"t=$(round(d.t[k], digits=0))  N=$(length(s))  fAB=$(round(fAB, digits=2))" :
				labels[i]) : ""
		ax = Axis3(fig[j, i], aspect = :data, azimuth = 0,
			elevation = view.perspective,
			xlabel = "", ylabel = "", zlabel = "",
			xticklabelsize = labelsize, yticklabelsize = labelsize,
			zticklabelsize = labelsize,
			titlevisible = showtime && j == 1, titlealign = :center,
			titlegap = 6, titlesize = labelsize, title = ttl)
		o = sort(keep; by = m -> s[m], rev = true)     # draw C first
		cells_mesh!(ax, x[o] .- mean(x), y[o] .- mean(y), z[o] .- zc, r[o], s[o]; ghost = gh, cmap = color_map)
		xlims!(ax, -bs, bs); ylims!(ax, -bs, bs); zlims!(ax, -bs, bs)
		if !shownumbers
			ax.xticklabelsize = 0; ax.yticklabelsize = 0; ax.zticklabelsize = 0
		end
		ax.zticklabelsvisible = view.show_z
		ax.xticklabelsvisible = view.show_x
		ax.yticklabelsvisible = false
	end
	savefig && save("$(filename).png", fig)
	fig
end;

# ---- notebook cell 6b8e0c35 ----------------------------------------------------
# --- 3D gallery: several runs and/or frames side by side ----------------------------

"""
	plot_gallery(panels; ncol = length(panels), views = (:top, :side), cut_top = true,
	             size_ratio = 0.7, labelsize = 22, title = nothing,
	             filename = "gallery", savefig = false, ghost = nothing)

3D panels in the plot_pancake look (Axis3 + meshscatter, C drawn first), for
several runs and/or frames at once. `panels` is a vector of NamedTuples
(d, k, label): d from load_frames, k the frame to draw, label the panel title.
Each panel gets one row per entry of `views`:

	:top    seen from above
	:side   seen from the side (shows the slab thickness)

With ghost = nothing (default) everything is opaque and, with cut_top, the top view
drops the cells above the mid-plane (a cutaway); ghost = GHOST draws the whole
aggregate with C translucent (unconfined aggregates, see draw_cells!).
"""
function plot_gallery(panels; ncol = length(panels), views = (:top, :side), cut_top = true,
		size_ratio = 0.7, labelsize = 22, title = nothing,
		filename = "gallery", savefig = false, ghost = nothing)
	nv = length(views)
	nrow = cld(length(panels), ncol)
	toff = title === nothing ? 0 : 1
	fig = Figure(fonts = TEX_FONTS, size = (round(Int, 520*size_ratio*ncol),
						 round(Int, 470*size_ratio*nv*nrow) + 70*toff), figure_padding = 20)
	title === nothing || Label(fig[1, 1:ncol], title; fontsize = labelsize + 6, font = :bold)
	for (p, pan) in enumerate(panels)
		rb, cc = divrem(p - 1, ncol)
		x, y, z, r, s = pan.d.x[pan.k], pan.d.y[pan.k], pan.d.z[pan.k], pan.d.r[pan.k], pan.d.cell_state[pan.k]
		xc, yc, zc = mean(x), mean(y), mean(z)
		bs = max(maximum(x) - minimum(x), maximum(y) - minimum(y)) / 2 + 2.0   # + margin: cells are not clipped
		for (jv, v) in enumerate(views)
			gh = ghost
			keep = (v == :top && cut_top && gh === nothing) ? findall(w -> w <= zc + 0.2, z) : collect(eachindex(z))
			ax = Axis3(fig[toff + rb*nv + jv, cc + 1], aspect = :data, azimuth = 0,
				elevation = v == :top ? pi/2 : 2*pi + pi/32,
				xlabel = "", ylabel = "", zlabel = "",
				xgridvisible = false, ygridvisible = false, zgridvisible = false,
				xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
				xticksvisible = false, yticksvisible = false, zticksvisible = false,
				xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false,
				titlevisible = jv == 1, title = pan.label, titlesize = labelsize, titlegap = 4)
			o = sort(keep; by = m -> s[m], rev = true)     # draw C first
			cells_mesh!(ax, x[o] .- xc, y[o] .- yc, z[o] .- zc, r[o], s[o]; ghost = gh)
			xlims!(ax, -bs, bs); ylims!(ax, -bs, bs); zlims!(ax, -bs, bs)
		end
	end
	title === nothing || rowgap!(fig.layout, 1, 60)     # room between the title and the panel labels
	savefig && save("$(filename).png", fig)
	fig
end;

# ---- notebook cell a7f3c2d1 ----------------------------------------------------
# --- final figures: 3D panels, the data behind them, the phase space, movies -------
#
# Every final figure can be redrawn from code-claude/figures/data alone, without
# re-running anything (replot_figures). Per figure, data/<figure>/ holds
#
#   title.txt       the figure title
#   panels.tsv      one row per panel: panel, row, col, views, label, file
#                   (tab-separated; "|" in a label is a line break)
#   <panel>.csv     one row per cell: x, y, z, r, cell_state -- the frame shown
#   <panel>.jld2    the same arrays plus where they come from (meta: run, frame, t, fB)
#
# and the phase space has data/phase_space/phase_space.csv (one row per run) plus
# its insets in panels.tsv, where the row and col columns hold the layers and N0 of
# the run the inset comes from. reports/FIGURES.md says how each figure was made.

PHASE_COLORS = Dict(                   # outcome classes of the phase space: the colours of your sketch (November 2026)
	"CENTRAL"   => Makie.to_color("#1fb3c8"),  # cyan
	"PATTERN"   => Makie.to_color("#ee82b4"),  # pink (ring or islands)
	"POLARISED" => Makie.to_color("#a9a9a9"),  # grey
	"OTHER"     => Makie.to_color("#e6e6e6"),  # light grey
);

"A darker shade of a colour, for text in the colour of a region."
darker(c) = RGBf(0.55c.r, 0.55c.g, 0.55c.b);

"The phase-space class of a classify_pattern outcome: rings and islands are both PATTERN."
phase_class(c) = c in ("RING", "ISLANDS") ? "PATTERN" : c in ("CENTRAL", "POLARISED") ? c : "OTHER";

"""
	save_panel(dir, panel, x, y, z, r, s; meta...)

Write the cells of one 3D panel as <panel>.csv (x, y, z, r, cell_state) and
<panel>.jld2 (the same arrays plus `meta`, e.g. run = file, frame = k, t, fB).
"""
function save_panel(dir, panel, x, y, z, r, s; meta...)
	mkpath(dir)
	open(joinpath(dir, "$(panel).csv"), "w") do io
		println(io, "x,y,z,r,cell_state")
		for i in eachindex(x)
			@printf(io, "%.4f,%.4f,%.4f,%.4f,%d\n", x[i], y[i], z[i], r[i], s[i])
		end
	end
	jldsave(joinpath(dir, "$(panel).jld2"); x = x, y = y, z = z, r = r, cell_state = s,
			meta = Dict(String(k) => v for (k, v) in meta))
	"$(panel).csv"
end;

"Read a panel written by save_panel: (x, y, z, r, cell_state)."
function load_panel(file)
	raw = readdlm(file, ','; header = true)[1]
	(x = Float64.(raw[:, 1]), y = Float64.(raw[:, 2]), z = Float64.(raw[:, 3]),
	 r = Float64.(raw[:, 4]), cell_state = Int.(raw[:, 5]))
end;

"""
	GHOST

Opacity of the C cells in the views of unconfined (polarising) aggregates (November
2026, your choice: 0.2, then 0.3, then 0.5 on Nov 3): the whole ball is drawn, B opaque and C translucent, so the B core
shows through the C cells around it. As in your original plot_pancake: one
meshscatter with RGBA colours and no `transparency = true` (whose order-independent
blending turns many translucent layers dark brown), B drawn first so the C spheres
in front of it blend over it. Needs GLMakie; CairoMakie draws each sphere as many
small triangles whose edges add up. Confined aggregates are drawn opaque with the
mid-plane cutaway (`ghost = nothing`, the default of draw_cells!).
"""
GHOST = 0.5;

"Colours of the cells with states `s`, all with opacity `a`."
cell_colors(s; a = 1f0) = [RGBAf(color_map[v], a) for v in s];

"""
	cells_mesh!(ax, x, y, z, r, s; ghost = nothing, cmap = color_map)

The cells as spheres in `ax`, one meshscatter. ghost = nothing: opaque, C drawn first.
ghost = a number: B drawn first and opaque, then C with that opacity (see GHOST).
"""
function cells_mesh!(ax, x, y, z, r, s; ghost = nothing, cmap = color_map)
	o = sortperm(s; rev = ghost === nothing)                   # opaque: C first; ghost: B first
	a = ghost === nothing ? 1f0 : Float32(ghost)
	isempty(o) || meshscatter!(ax, x[o], y[o], z[o]; markersize = r[o],
							   color = [RGBAf(cmap[v], v == 2 ? 1f0 : a) for v in s[o]])
	ax
end;

"Half-width of a box that fits the aggregate in the plane, with a margin."
box_halfwidth(c) = max(maximum(c.x) - minimum(c.x), maximum(c.y) - minimum(c.y)) / 2 + 2.0;

"A panel title line: plain text, or LaTeX (text with inline \$...\$ math) when it has a \$."
panel_text(s) = occursin('$', s) ? Makie.LaTeXString(s) : s;

"Camera elevation of the :side view (20 degrees: some of the top shows) and of the :oblique view (45)."
SIDE_ELEVATION = pi/9;
OBLIQUE_ELEVATION = pi/4;

"""
	draw_cells!(pos, c; view = :top, ghost = nothing, cut_top = true, label = "", labelsize = 18, bs = nothing)

One 3D panel in the plot_pancake look (Axis3 + meshscatter) at the figure position
`pos`, for the cells `c` (anything with x, y, z, r, cell_state):

	:top      seen from above
	:side     seen from the side, tilted by SIDE_ELEVATION so some of the top shows
	:front    seen from the side, along y (x runs left to right), level
	:oblique  a three-quarter view from above (OBLIQUE_ELEVATION)

With ghost = nothing (the default) the cells are opaque, and with cut_top the top
view leaves out the cells above the mid-plane, so the B layout shows on the cut face
(confined aggregates). With ghost = GHOST every view shows the whole aggregate with C
translucent (unconfined aggregates, where B polarises; November 2026).

`bs` is the half-width of the box (default: fits the aggregate); give every panel
of a figure the same `bs` to draw them at the same scale (the box is then a cube, drawn at
a scale that does not depend on the view). Without `bs`, a flat aggregate seen from the
side or obliquely gets a flat box, so it fills the panel.

`label`: the title, one line per "\n", each plain text or LaTeX (`panel_text`), above the axis.
"""
function draw_cells!(pos, c; view = :top, ghost = nothing, cut_top = true, label = "", labelsize = 18, labelcolor = :black, bs = nothing)
	x, y, z, r, s = c.x, c.y, c.z, c.r, c.cell_state
	xc, yc, zc = mean(x), mean(y), mean(z)
	same_scale = bs !== nothing
	bs = same_scale ? bs : box_halfwidth(c)
	gh = ghost
	keep = (view == :top && cut_top && gh === nothing) ? findall(w -> w <= zc + 0.2, z) : collect(eachindex(z))
	flat = !same_scale && view != :top && maximum(z) - minimum(z) < 1.2bs
	lines_ = isempty(label) ? String[] : String.(split(label, "\n"))
	g = isempty(lines_) ? nothing : GridLayout(pos)
	for (i, l) in enumerate(lines_)
		Label(g[i, 1], panel_text(l); fontsize = labelsize, color = labelcolor, tellwidth = false)
	end
	ax = Axis3(g === nothing ? pos : g[length(lines_) + 1, 1], aspect = :data,
		azimuth = view == :oblique ? pi/4 : view == :front ? pi/2 : 0,
		elevation = view == :top ? pi/2 : view == :side ? SIDE_ELEVATION : view == :front ? 2*pi + pi/32 : OBLIQUE_ELEVATION,
		viewmode = same_scale ? :fit : :fitzoom,
		xlabel = "", ylabel = "", zlabel = "",
		xgridvisible = false, ygridvisible = false, zgridvisible = false,
		xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
		xticksvisible = false, yticksvisible = false, zticksvisible = false,
		xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false)
	g === nothing || rowgap!(g, 2)
	cells_mesh!(ax, x[keep] .- xc, y[keep] .- yc, z[keep] .- zc, r[keep], s[keep]; ghost = gh)
	xlims!(ax, -bs, bs); ylims!(ax, -bs, bs)
	flat ? zlims!(ax, minimum(z) - zc - 1.2, maximum(z) - zc + 1.2) : zlims!(ax, -bs, bs)
	ax
end;

"""
	plot_panels(dir; size_ratio = 0.6, labelsize = 18)

A final 3D figure from its data folder (title.txt, panels.tsv, <panel>.csv). Each
panel sits at (row, col) of panels.tsv and gets one axis row per view
("top", "side", "front", "oblique", or several joined by "+"; "arrow" draws a down
arrow instead of cells, e.g. between the two rows of the fusion figure). An optional
7th column `bs` fixes the box half-width of that panel (the same value in every panel
= the same scale); an 8th column `ghost` = 1 draws the panel whole with C translucent
(GHOST: unconfined aggregates), else opaque with the mid-plane cutaway. An empty title
or label draws none; labels can hold LaTeX (\$...\$). All text in LaTeX's fonts.
"""
function plot_panels(dir; size_ratio = 0.6, labelsize = 18)
	lines_ = filter(l -> !isempty(strip(l)), readlines(joinpath(dir, "panels.tsv")))
	pans = [split(l, '\t') for l in lines_[2:end]]       # panel, row, col, views, label, file[, bs]
	pbs(p) = length(p) >= 7 && !isempty(strip(p[7])) ? parse(Float64, p[7]) : nothing
	pghost(p) = length(p) >= 8 && strip(p[8]) == "1" ? GHOST : nothing
	views(p) = Symbol.(split(p[4], "+"))
	nrow = maximum(parse(Int, p[2]) for p in pans)
	ncol = maximum(parse(Int, p[3]) for p in pans)
	# every figure row is as tall as its panel with the most views
	nv = [maximum([length(views(p)) for p in pans if parse(Int, p[2]) == i]; init = 1) for i in 1:nrow]
	start = cumsum(vcat(1, nv[1:end-1]))
	title = isfile(joinpath(dir, "title.txt")) ? strip(read(joinpath(dir, "title.txt"), String)) : ""
	toff = isempty(title) ? 0 : 1
	# a side view is flat: an axis row with only side views gets 45 % of the height, one
	# with only arrows 30 %
	h = round(Int, 470*size_ratio)
	rowkind(i, jv) = (vs = [views(p)[jv] for p in pans if parse(Int, p[2]) == i && length(views(p)) >= jv];
					  all(==(:arrow), vs) ? 0.3 : all(in((:side, :arrow)), vs) ? 0.45 : 1.0)
	frac = [rowkind(i, jv) for i in 1:nrow for jv in 1:nv[i]]
	fig = Figure(size = (round(Int, 520*size_ratio*ncol), sum(round(Int, f*h) for f in frac) + 70*toff),
				 figure_padding = 20, fonts = TEX_FONTS)
	# the title shrinks to fit a narrow figure
	W = round(Int, 520*size_ratio*ncol)
	isempty(title) || Label(fig[1, 1:ncol], title; fontsize = min(labelsize + 6, round(Int, 1.8W / length(title))), font = :bold)
	for p in pans
		row, col = parse(Int, p[2]), parse(Int, p[3])
		if views(p) == [:arrow]
			Label(fig[toff + start[row], col], L"\downarrow"; fontsize = 3labelsize, tellwidth = false)
			continue
		end
		c = load_panel(joinpath(dir, p[6]))
		for (jv, v) in enumerate(views(p))
			draw_cells!(fig[toff + start[row] + jv - 1, col], c; view = v, bs = pbs(p), ghost = pghost(p),
						label = jv == 1 ? replace(p[5], "|" => "\n") : "", labelsize = labelsize)
		end
	end
	for (j, f) in enumerate(frac)
		f < 1 && rowsize!(fig.layout, toff + j, Fixed(round(Int, f*h)))
	end
	isempty(title) || rowgap!(fig.layout, 1, 50)
	fig
end;

"""
	write_panels(dir, title, panels; same_scale = false)

Write a final figure's data folder. `panels` is a vector of NamedTuples
(panel, row, col, views, label, cells, meta[, ghost]): cells has x, y, z, r, cell_state
(e.g. one frame of load_frames), meta a NamedTuple saved in the .jld2, ghost = true
to draw the panel whole with C translucent (unconfined aggregates). A panel with
views = "arrow" has no cells (cells = nothing).
`same_scale` writes one box half-width (that of the largest panel) for every
panel, so plot_panels draws them at the same scale.
"""
function write_panels(dir, title, panels; same_scale = false)
	mkpath(dir)
	write(joinpath(dir, "title.txt"), title)
	bs = same_scale ? @sprintf("%.3f", maximum(box_halfwidth(p.cells) for p in panels if p.cells !== nothing)) : ""
	open(joinpath(dir, "panels.tsv"), "w") do io
		println(io, "panel\trow\tcol\tviews\tlabel\tfile\tbs\tghost")
		for p in panels
			c = p.cells
			f = c === nothing ? "" : save_panel(dir, p.panel, c.x, c.y, c.z, c.r, c.cell_state; p.meta...)
			gh = hasproperty(p, :ghost) && p.ghost ? "1" : ""
			println(io, join((p.panel, p.row, p.col, p.views, replace(p.label, "\n" => "|"), f, bs, gh), '\t'))
		end
	end
	dir
end;

PHASE_LABELS = Dict("CENTRAL" => "central cluster", "PATTERN" => "pattern (ring or islands)",
					"POLARISED" => "polarised", "OTHER" => "other");

"The names of the outcomes in the figures since November 2026 (legend and inset titles)."
PHASE_NAMES = Dict("CENTRAL" => "B-centroid", "PATTERN" => "Clusters", "POLARISED" => "Polarisation", "OTHER" => "other");

"Figure fonts: LaTeX's Computer Modern for all the text (November 2026)."
TEX_FONTS = (; regular = Makie.texfont(), bold = Makie.texfont(:bold), italic = Makie.texfont(:italic));

"""
	phase_space!(pos, rows; insets = NamedTuple[], free = NamedTuple[], title = "",
	             x = :log, y = :layers, colors = PHASE_COLORS, background = 0.5, fontsize = 28,
	             pie_scale = 1.0, region_labels = nothing, contours = nothing,
	             contour_levels = AREA_LEVELS, axis = (;))

The chart of plot_phase_space, drawn at the figure position `pos` (so it can be
composed into a bigger figure). Returns the axis.

* the background is the majority outcome of the nearest point; a tie is broken by the
  majority at the later readout (the rows' `late` class, if given): what it becomes
* `pies`, `thresholds`, `letters`: the seeds' pies, the threshold lines and the inset
  letters (all off since November 2026)
* `x`  :log         log10 N0, as in your sketch;  :linear  N0
* `y`  :layers      the plate gap in cell layers, on a log scale (log2: each doubling
                    of the gap is one step; the ticks show the number of layers);
                    :layers_lin  the same on a linear scale;
                    :diameter  the pancake diameter at squeezing, in cell diameters,
                    log scale (rows need D)
* `background`      opacity of the region colours
* `pie_scale`       scales the pies (smaller for a sketch-like chart)
* `region_labels`   class => text, written in the middle of that class's region
* `contours`        rows (layers, N0, A) with the pancake area at squeezing A (um^2):
                    labelled iso-area lines at `contour_levels` (width `contour_width`,
                    labels `contour_labelsize`),
                    interpolated in log A along each row of layers and then between rows
* `axis`            more Axis attributes (e.g. bbox, to place it by hand)
"""
AREA_LEVELS = [2.5e3, 5e3, 1e4, 2e4, 4e4, 8e4, 1.6e5];

"An area in um^2 as a short contour label: 2.5k, 10k, 160k."
area_label(A) = A >= 1e3 ? @sprintf("%gk", round(A / 1e3; sigdigits = 2)) : @sprintf("%g", round(A; sigdigits = 2));

function phase_space!(pos, rows; insets = NamedTuple[], free = NamedTuple[], title = "",
		x = :log, y = :layers, colors = PHASE_COLORS, background = 0.5, fontsize = 28,
		pie_scale = 1.0, region_labels = nothing, contours = nothing, contour_levels = AREA_LEVELS,
		contour_width = 3.0, contour_labelsize = fontsize - 8, pies = false, thresholds = false, letters = false, axis = (;))
	order = ["CENTRAL", "PATTERN", "POLARISED", "OTHER"]
	pts = unique([(Float64(r.layers), Int(r.N0)) for r in rows])
	outc = Dict(p => [phase_class(r.class) for r in rows if (r.layers, r.N0) == p] for p in pts)
	late = Dict(p => [phase_class(r.late) for r in rows if (r.layers, r.N0) == p && hasproperty(r, :late)] for p in pts)
	# a point sits at the mean pancake diameter of its seeds
	Dm = y == :diameter ? Dict(p => mean(r.D for r in rows if (r.layers, r.N0) == p) for p in pts) : Dict()
	# the unconfined runs: one more row above the thickest gap (layers = Inf)
	for N in unique([r.N0 for r in free])
		push!(pts, (Inf, N)); outc[(Inf, N)] = [phase_class(r.class) for r in free if r.N0 == N]
		late[(Inf, N)] = [phase_class(r.late) for r in free if r.N0 == N && hasproperty(r, :late)]
	end
	conf = filter(p -> isfinite(p[1]), pts)
	# the majority outcome of a point; a tie goes to the tied class most frequent later
	function major(p)
		n = [count(==(c), outc[p]) for c in order]
		tied = findall(==(maximum(n)), n)
		length(tied) == 1 && return order[tied[1]]
		l = [count(==(order[t]), get(late, p, String[])) for t in tied]
		order[tied[argmax(l)]]
	end
	majo = Dict(p => major(p) for p in pts)

	# plot coordinates: N0 across 10 units, the confinement up 8, the free row at 9.4
	fx(N) = x == :log ? log10(N) : float(N)
	fy(p) = y == :layers ? log2(p[1]) : y == :layers_lin ? float(p[1]) : log2(Dm[p])
	x0, x1 = x == :log ? (log10(minimum(last.(pts))), log10(maximum(last.(pts)))) : (0.0, float(maximum(last.(pts))))
	y0, y1 = extrema(fy.(conf))
	Yv(v) = 8 * (v - y0) / (y1 - y0)
	X(N) = 10 * (fx(N) - x0) / (x1 - x0)
	Y(p) = isfinite(p[1]) ? Yv(fy(p)) : 9.4
	top = any(p -> !isfinite(p[1]), pts) ? 10.1 : 8.7

	if x == :log
		xt = filter(t -> x0 - 0.05 <= t <= x1 + 0.05, 1.5:0.5:4.5)
		xticks = (10 .* (xt .- x0) ./ (x1 - x0), [@sprintf("%g", t) for t in xt]); xlabel = L"\log_{10} N_0"
	else
		xt = 0:1000:maximum(last.(pts))
		xticks = (10 .* xt ./ x1, string.(xt)); xlabel = L"N_0"
	end
	if y in (:layers, :layers_lin)
		# log: ticks at 1, 2, 3, 4, 6, 8, 12 layers, labelled with the number of layers (not
		# its log); linear: at every sampled number of layers
		Ls = y == :layers ? filter(L -> y0 - 1e-9 <= log2(L) <= y1 + 1e-9, [1, 2, 3, 4, 6, 8, 12, 16]) : sort(unique(first.(conf)))
		yticks = (vcat(Yv.(y == :layers ? log2.(Ls) : float.(Ls)), top > 9 ? [9.4] : Float64[]),
				  vcat([@sprintf("%g", L) for L in Ls], top > 9 ? ["unconfined"] : String[]))
		ylabel = y == :layers ? L"n_\mathrm{layers}\;\mathrm{(log\ scale)}" : L"n_\mathrm{layers}"
	else
		yt = filter(t -> 2.0^y0 / 1.1 <= t <= 2.0^y1 * 1.1, [4, 5, 7, 10, 15, 20, 30, 40, 60, 80, 120])
		yticks = (vcat(Yv.(log2.(yt)), top > 9 ? [9.4] : Float64[]), vcat(string.(yt), top > 9 ? ["unconfined"] : String[]))
		ylabel = L"D \;\mathrm{(cell\ diameters,\ log\ scale)}"
	end
	ax = Axis(pos; aspect = DataAspect(), title = title, titlesize = fontsize + 2,
		xlabel = xlabel, ylabel = ylabel, xlabelsize = fontsize, ylabelsize = fontsize,
		xticks = xticks, yticks = yticks, xgridvisible = false, ygridvisible = false,
		xticklabelsize = fontsize - 4, yticklabelsize = fontsize - 4, axis...)
	limits!(ax, -0.7, 10.7, -0.7, top)

	# background: majority outcome of the nearest point
	xs = range(-0.7, 10.7, length = 460); ys = range(-0.7, top, length = 400)
	P = [(X(p[2]), Y(p)) for p in pts]
	grid = [findfirst(==(majo[pts[argmin([(q[1] - u)^2 + (q[2] - v)^2 for q in P])]]), order) for u in xs, v in ys]
	fade(c, a) = RGBAf(c.r, c.g, c.b, a)
	heatmap!(ax, xs, ys, grid; colormap = cgrad([fade(colors[c], background) for c in order], length(order); categorical = true),
			 colorrange = (0.5, length(order) + 0.5))
	# the axis break between the plates and the unconfined row
	top > 9 && poly!(ax, Rect2f(-0.7, 8.55, 11.4, 0.3); color = :white)

	# thresholds, per number of layers (the unconfined row included)
	thresholds && for (which, style) in (("PATTERN", :dash), ("POLARISED", :dot))
		tx = Float64[]; ty = Float64[]
		for L in sort(unique(first.(pts)))
			row = sort(filter(p -> p[1] == L, pts); by = last)
			m = [majo[p] == which for p in row]
			if which == "PATTERN"     # patterns from some N0 on, up to the largest
				i = findlast(!, m)                     # the last N0 without a pattern
				i == length(row) && continue
				i === nothing ? (push!(tx, max(X(row[1][2]) - 0.5, -0.3)); push!(ty, Y(row[1]))) :
					(push!(tx, (X(row[i][2]) + X(row[i+1][2])) / 2); push!(ty, (Y(row[i]) + Y(row[i+1])) / 2))
			else                      # polarised up to some N0
				i = findfirst(!, m)                    # the first N0 that does not polarise
				i == 1 && continue
				i === nothing ? (push!(tx, min(X(row[end][2]) + 0.5, 10.3)); push!(ty, Y(row[end]))) :
					(push!(tx, (X(row[i-1][2]) + X(row[i][2])) / 2); push!(ty, (Y(row[i-1]) + Y(row[i])) / 2))
			end
		end
		length(tx) >= 1 && lines!(ax, tx, ty; color = colors[which], linewidth = 3, linestyle = style)
		length(tx) >= 1 && scatter!(ax, tx, ty; color = colors[which], markersize = 9)
	end

	# one pie per point: the seeds' outcomes; smaller where points crowd
	rad = Dict(p => pie_scale * clamp(0.45 * minimum([hypot(X(p[2]) - X(q[2]), Y(p) - Y(q)) for q in pts if q != p]; init = 1.0), 0.13, 0.42)
			   for p in pts)
	pies && for p in pts
		o = outc[p]; n = length(o); a0 = pi/2; cx, cy = X(p[2]), Y(p)
		for c in order
			k = count(==(c), o); k == 0 && continue
			a1 = a0 - 2pi * k / n
			th = range(a0, a1, length = max(2, ceil(Int, 40 * k / n)))
			poly!(ax, Point2f[(cx, cy); [(cx + rad[p]*cos(t), cy + rad[p]*sin(t)) for t in th]];
				  color = colors[c], strokecolor = :white, strokewidth = k == n ? 0 : 1.5)
			a0 = a1
		end
	end

	# iso-area lines of the pancake at squeezing, over the pies so the labels show: log A
	# interpolated (and extrapolated to the edges) along each row of layers, then linearly
	# between the rows
	if contours !== nothing && y in (:layers, :layers_lin)
		cpts = unique([(Float64(c.layers), Int(c.N0)) for c in contours])
		lA = Dict(p => mean(log10(c.A) for c in contours if (c.layers, c.N0) == p) for p in cpts)
		crow = Tuple{Float64,Vector{Float64},Vector{Float64}}[]          # (Y, X's, log A's) per row
		for L in sort(unique(first.(cpts)))
			ps = sort(filter(p -> p[1] == L, cpts); by = last)
			length(ps) >= 2 && push!(crow, (Yv(fy(ps[1])), [X(p[2]) for p in ps], [lA[p] for p in ps]))
		end
		along(u, xs_, vs_) = (i = clamp(searchsortedlast(xs_, u), 1, length(xs_) - 1);
							  vs_[i] + (vs_[i+1] - vs_[i]) * (u - xs_[i]) / (xs_[i+1] - xs_[i]))
		cx = range(-0.7, 10.7, length = 230); cy = range(crow[1][1], crow[end][1], length = 200)   # only between sampled rows
		Z = fill(NaN, length(cx), length(cy))
		for (i, u) in enumerate(cx), (j, v) in enumerate(cy)
			# below the first row and above the last one: that row's values (no NaN,
			# which breaks the contour labels)
			k = something(findlast(r -> r[1] <= v + 1e-9, crow), 1)
			a = along(u, crow[k][2], crow[k][3])
			Z[i, j] = (k == length(crow) || v < crow[1][1]) ? a :
				a + (along(u, crow[k+1][2], crow[k+1][3]) - a) * (v - crow[k][1]) / (crow[k+1][1] - crow[k][1])
		end
		contour!(ax, cx, cy, Z; levels = log10.(contour_levels), color = (:black, 0.7), linewidth = contour_width,
				 labels = true, labelsize = contour_labelsize, labelcolor = (:black, 0.85),
				 labelformatter = v -> area_label(10^v))
	end

	# region names, below the unconfined row: inside each region (the whole text width
	# within it), at the spot farthest from any pie and from the names already placed,
	# pulled a little towards the middle of the region
	if region_labels !== nothing
		placed = Tuple{Float64,Float64}[]
		for c in order
			haskey(region_labels, c) || continue
			txt = region_labels[c]; k = findfirst(==(c), order)
			hw = 0.12 * length(txt)              # about half the text width
			di = ceil(Int, hw / step(xs)); dj = ceil(Int, 0.3 / step(ys))
			inside(i, j) = all(grid[clamp(i + a, 1, length(xs)), clamp(j + b, 1, length(ys))] == k
							   for a in (-di, 0, di), b in (-dj, 0, dj))
			in_ = [(u, v) for (i, u) in enumerate(xs), (j, v) in enumerate(ys)
				   if grid[i, j] == k && hw - 0.5 < u < 10.5 - hw && 0 < v < 8.2 && inside(i, j)]
			isempty(in_) && continue
			m = (mean(first.(in_)), mean(last.(in_)))
			score(q) = min(minimum(hypot(q[1] - P_[1], q[2] - P_[2]) for P_ in P),
						   minimum((hypot((q[1] - p_[1]) / 3, q[2] - p_[2]) for p_ in placed); init = Inf)) -
					   0.15hypot(q[1] - m[1], q[2] - m[2])
			at = in_[argmax(score.(in_))]
			push!(placed, at)
			cc = colors[c]
			text!(ax, at...; text = txt, align = (:center, :center), fontsize = fontsize + 1,
				  font = :bold, color = RGBf(0.5cc.r, 0.5cc.g, 0.5cc.b))
		end
	end

	# the insets' letters at their points
	letters && for (i, ins) in enumerate(insets)
		p = (ins.layers, ins.N0); haskey(outc, p) || continue
		text!(ax, X(p[2]) + rad[p] + 0.06, Y(p) + rad[p]; text = string('A' + i - 1), fontsize = fontsize + 2,
			  font = :bold, align = (:left, :bottom))
	end
	ax
end;

"""
	plot_phase_space(rows; insets = NamedTuple[], free = NamedTuple[], title = "",
	                 x = :log, y = :layers, colors = PHASE_COLORS, labels = PHASE_LABELS)

The phase space. `rows` has one entry per run with fields layers, N0, class (the
class at the readout, classify_pattern) and, for y = :diameter, D.

* x = :log (log10 N0) or :linear (N0); y = :layers (plate gap) or :diameter (the
  pancake diameter at squeezing); see phase_space!
* every (layers, N0) point is a pie of the outcomes of its seeds
* the background takes the colour of the majority outcome of the nearest point
* two thresholds per number of layers, where they exist: the N0 above which
  patterns form (dashed), and the N0 below which the aggregate polarises (dotted),
  half-way between the two points that bracket the change; half a step outside
  the grid when every point is past it
* `insets`: (cells, label, layers, N0, cls) drawn in 3D on the right (top and side
  view), titled in the colour of their region; the polarised one whole with C
  translucent (GHOST), the others with the mid-plane cutaway
* `free`: unconfined runs (N0, class), drawn as an extra top row "unconfined" --
  what a gap much wider than the aggregate tends to
"""
function plot_phase_space(rows; insets = NamedTuple[], free = NamedTuple[], title = "",
		x = :log, y = :layers, colors = PHASE_COLORS, labels = PHASE_NAMES, contours = nothing)
	order = ["CENTRAL", "PATTERN", "POLARISED", "OTHER"]
	fig = Figure(size = (isempty(insets) ? 1000 : 1600, 900), figure_padding = 20, fonts = TEX_FONTS)
	phase_space!(fig[1, 1], rows; insets = insets, free = free, title = title, x = x, y = y, colors = colors,
				 contours = contours)
	if !isempty(insets)
		# one example per outcome, titled with its name in the colour of its region
		g = fig[1, 2] = GridLayout()
		for (i, ins) in enumerate(insets)
			cls = hasproperty(ins, :cls) ? ins.cls : ""
			col = haskey(colors, cls) ? darker(colors[cls]) : :black
			gh = cls == "POLARISED" ? GHOST : nothing
			draw_cells!(g[i, 1], ins.cells; view = :top, ghost = gh, label = ins.label, labelsize = 30, labelcolor = col)
			draw_cells!(g[i, 2], ins.cells; view = :side, ghost = gh)
		end
		colsize!(fig.layout, 2, Relative(0.36))
	end
	# the legend: only the outcomes that occur
	seen = Set(phase_class(r.class) for r in vcat(collect(rows), collect(free)))
	shown = filter(c -> c in seen, ["POLARISED", "CENTRAL", "PATTERN", "OTHER"])   # the order of the insets
	els = Any[PolyElement(color = colors[c]) for c in shown]; lbs = [labels[c] for c in shown]
	if contours !== nothing && y in (:layers, :layers_lin)
		push!(els, LineElement(color = (:black, 0.7), linewidth = 3)); push!(lbs, "pancake area at squeezing (µm²)")
	end
	Legend(fig[2, 1], els, lbs; orientation = :horizontal, framevisible = false, labelsize = 24)
	fig
end;

"""
	record_pancake_frames(d, color_map, ks; views = (:top, :side), fps = 8,
	                      size_ratio = 0.8, filename = "movie.mp4")

record_pancake for the array format of load_frames: one column per view, one
movie frame per index in `ks`, the time (hours after `t0`, e.g. the start of
differentiation) and B fraction in the title.
ghost = nothing (default): opaque cells and the mid-plane cutaway in the top view
(confined runs); ghost = GHOST: the whole aggregate with C translucent (unconfined
runs; needs GLMakie). Cells are padded to
the largest frame (unused ones get size 0) so the plot can be updated in place
while cells divide.
"""
function record_pancake_frames(d, color_map, ks; views = (:top, :side), fps = 8,
		size_ratio = 0.8, filename = "movie.mp4", ghost = nothing, t0 = 0.0)
	nmax = maximum(length(d.x[k]) for k in ks)
	bs = maximum(max(maximum(d.x[k]) - minimum(d.x[k]), maximum(d.y[k]) - minimum(d.y[k])) for k in ks) / 2 + 2.0
	fig = Figure(size = (round(Int, 600*size_ratio*length(views)), round(Int, 600*size_ratio)), figure_padding = 20,
				 fonts = TEX_FONTS)
	ttl = Observable("")
	Label(fig[0, 1:length(views)], ttl; fontsize = 24)
	# one sphere plot per part: :all (opaque; in the top view without ghost, the cutaway),
	# or with ghost :B (opaque, drawn first) and :C (translucent, blended over B; see GHOST)
	parts(v) = ghost !== nothing ? (:B, :C) : (:all,)
	frame_data(k, v, part) = begin
		x, y, z, r, s = d.x[k], d.y[k], d.z[k], d.r[k], d.cell_state[k]
		zc = mean(z)
		show_ = [part == :all ? !(v == :top && z[i] > zc + 0.2) : part == :B ? s[i] == 2 : s[i] != 2
				 for i in eachindex(z)]
		pad = nmax - length(x)
		P = [Point3f(x[i] - mean(x), y[i] - mean(y), z[i] - zc) for i in eachindex(x)]
		C = [RGBAf(color_map[v_], part == :C ? ghost : 1f0) for v_ in s]
		(vcat(P, fill(Point3f(0, 0, 0), pad)),
		 vcat([show_[i] ? Float32(r[i]) : 0f0 for i in eachindex(r)], zeros(Float32, pad)),
		 vcat(C, fill(RGBAf(color_map[3], 0f0), pad)))
	end
	obs = []
	for (j, v) in enumerate(views)
		ax = Axis3(fig[1, j], aspect = :data, azimuth = 0,
			elevation = v == :top ? pi/2 : SIDE_ELEVATION,
			xgridvisible = false, ygridvisible = false, zgridvisible = false,
			xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
			xticksvisible = false, yticksvisible = false, zticksvisible = false,
			xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false,
			xlabel = "", ylabel = "", zlabel = "", title = String(v), titlesize = 20)
		for part in parts(v)
			P, R, C = frame_data(ks[1], v, part)
			Po = Observable(P); Ro = Observable(R); Co = Observable(C)
			meshscatter!(ax, Po; markersize = Ro, color = Co)
			push!(obs, (v, part, Po, Ro, Co))
		end
		xlims!(ax, -bs, bs); ylims!(ax, -bs, bs); zlims!(ax, -bs, bs)
	end
	record(fig, filename, ks; framerate = fps) do k
		for (v, part, Po, Ro, Co) in obs
			P, R, C = frame_data(k, v, part)
			Po[] = P; Ro[] = R; Co[] = C        # padded: the lengths never change
		end
		fB = count(==(2), d.cell_state[k]) / length(d.cell_state[k])
		ttl[] = @sprintf("t = %.1f h   fB = %.2f   N = %d", max(0, d.t[k] - t0) / 15,   # 15 time units = 1 h, from t0
			 fB, length(d.x[k]))
	end
	filename
end;

"""
	replot_figures(figdir; save_png = true)

Redraw every final figure from figdir/data alone: each data/<figure>/ with a
panels.tsv becomes figdir/<figure>.png; each data/phase_space*/ (phase_space.csv,
free.csv, panels.tsv for the insets, readout.txt naming the class column to draw and,
on its second line, the column that breaks ties between seeds) becomes three charts:
figdir/phase_space*.png (x = log10 N0, y = number of layers, linear), *_lin.png (both
axes linear) and *_diameter.png (y = pancake diameter at squeezing, when phase_space.csv
has the column D).
"""
function replot_figures(figdir; save_png = true)
	figs = Dict{String,Any}()
	for name in sort(readdir(joinpath(figdir, "data")))
		dir = joinpath(figdir, "data", name)
		isdir(dir) || continue
		if startswith(name, "phase_space") && isfile(joinpath(dir, "phase_space.csv"))
			# which readout to draw: readout.txt names the class column (default class_02) and the
			# tie-breaking one (default: the same readout at fB 0.12, class_..._012)
			rl = isfile(joinpath(dir, "readout.txt")) ? filter(!isempty, strip.(readlines(joinpath(dir, "readout.txt")))) : String[]
			ccol = isempty(rl) ? "class_02" : String(rl[1])
			raw, hdr = readdlm(joinpath(dir, "phase_space.csv"), ','; header = true)
			col = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
			lcol = length(rl) >= 2 ? String(rl[2]) : replace(ccol, r"_02$" => "_012")
			rows = [(layers = Float64(raw[i, col["layers"]]), N0 = Int(raw[i, col["N0"]]),
					 class = String(raw[i, col[ccol]]), D = haskey(col, "D") ? Float64(raw[i, col["D"]]) : NaN,
					 late = haskey(col, lcol) ? String(raw[i, col[lcol]]) : String(raw[i, col[ccol]]))
					for i in 1:size(raw, 1)]
			# iso-area lines of the pancake at squeezing, when phase_space.csv has them
			cont = haskey(col, "A_um2") ?
				[(layers = r.layers, N0 = r.N0, A = Float64(raw[i, col["A_um2"]])) for (i, r) in enumerate(rows)] : nothing
			ins = NamedTuple[]
			if isfile(joinpath(dir, "panels.tsv"))
				for l in readlines(joinpath(dir, "panels.tsv"))[2:end]
					isempty(strip(l)) && continue
					p = split(l, '\t')       # panel, layers, N0, views, label, file
					push!(ins, (cells = load_panel(joinpath(dir, p[6])), label = replace(p[5], "|" => "\n"),
								layers = parse(Float64, p[2]), N0 = parse(Int, p[3]), cls = uppercase(String(p[1]))))
				end
			end
			free = NamedTuple[]
			if isfile(joinpath(dir, "free.csv"))      # unconfined runs: N0, seed, class_02, ...
				raw, hdr = readdlm(joinpath(dir, "free.csv"), ','; header = true)
				fc = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
				free = [(N0 = Int(raw[i, fc["N0"]]), class = String(raw[i, fc[ccol]]),
						 late = haskey(fc, lcol) ? String(raw[i, fc[lcol]]) : String(raw[i, fc[ccol]])) for i in 1:size(raw, 1)]
			end
			title = isfile(joinpath(dir, "title.txt")) ? strip(read(joinpath(dir, "title.txt"), String)) : ""
			# the number of layers on a linear axis (Nov 3), and a version with both axes linear
			for (suffix, x, y) in (("", :log, :layers_lin), ("_lin", :linear, :layers_lin), ("_diameter", :log, :diameter))
				y == :diameter && !haskey(col, "D") && continue
				# the outcomes by the names of the figures (PHASE_NAMES); the October readout counts rings as clusters too
				lab = PHASE_NAMES
				figs[name * suffix] = plot_phase_space(rows; insets = ins, free = free, title = title, x = x, y = y, labels = lab,
													   contours = cont)
				save_png && save(joinpath(figdir, "$(name)$(suffix).png"), figs[name * suffix])
			end
		elseif isfile(joinpath(dir, "panels.tsv"))
			figs[name] = plot_panels(dir)
			save_png && save(joinpath(figdir, "$(name).png"), figs[name])
		end
		haskey(figs, name) && println("redrew $(name)")
	end
	figs
end;

# ---- notebook cell 9557e7da ----------------------------------------------------
# --- island detection ---------------------------------------------------------

"""
	cluster_sizes(px, py, pz, cutoff)

Sizes of the single-linkage clusters of the given points, descending. Two cells
are in the same cluster if they are within `cutoff`. O(N) via a uniform grid.
"""
function cluster_sizes(px, py, pz, cutoff)
	n = length(px)
	n == 0 && return Int[]
	lo = (minimum(px), minimum(py), minimum(pz))
	cell(i) = (floor(Int,(px[i]-lo[1])/cutoff), floor(Int,(py[i]-lo[2])/cutoff),
			   floor(Int,(pz[i]-lo[3])/cutoff))
	bucket = Dict{NTuple{3,Int},Vector{Int}}()
	for i in 1:n; push!(get!(bucket, cell(i), Int[]), i); end
	parent = collect(1:n)
	function root(a)
		while parent[a] != a; parent[a] = parent[parent[a]]; a = parent[a]; end
		a
	end
	c2 = cutoff^2
	for i in 1:n
		ci = cell(i)
		for dx in -1:1, dy in -1:1, dz in -1:1
			b = get(bucket, (ci[1]+dx, ci[2]+dy, ci[3]+dz), nothing)
			b === nothing && continue
			for j in b
				j <= i && continue
				if (px[i]-px[j])^2 + (py[i]-py[j])^2 + (pz[i]-pz[j])^2 < c2
					a1, a2 = root(i), root(j)
					a1 != a2 && (parent[a1] = a2)
				end
			end
		end
	end
	cnt = Dict{Int,Int}()
	for i in 1:n; rt = root(i); cnt[rt] = get(cnt, rt, 0) + 1; end
	sort(collect(values(cnt)); rev = true)
end;

"Minimum cells for a cluster to count as an island at `layers` layers."
nmin_island(layers) = max(8, round(Int, 10*layers));

"""
	island_stats(x, y, z, s, layers; cut = 2.6)

One frame. Returns `(n_isl, frac_max, clustered, nAB, top)`:

* `n_isl`     number of A,B clusters of at least `nmin_island(layers)` cells
* `frac_max`  largest cluster / all A,B cells   (1 => still one blob)
* `clustered` fraction of A,B cells living in real islands at all
              (low => salt-and-pepper dispersion, which is NOT a pattern even
              though `frac_max` is then tiny -- this is the discriminator)
* `top`       the eight largest cluster sizes

The island threshold scales with `layers` so a monolayer island is not
penalised for containing fewer cells at the same in-plane footprint.
"""
function island_stats(x, y, z, s, layers; cut = 2.6)
	sel = findall(i -> s[i] == 1 || s[i] == 2, eachindex(s))
	nAB = length(sel)
	nAB == 0 && return (n_isl = 0, frac_max = NaN, clustered = 0.0, nAB = 0, top = Int[])
	cs  = cluster_sizes(x[sel], y[sel], z[sel], cut)
	big = filter(>=(nmin_island(layers)), cs)
	(n_isl = length(big), frac_max = cs[1]/nAB,
	 clustered = sum(big; init = 0)/nAB, nAB = nAB, top = cs[1:min(end,8)])
end;

"""
	classify_islands(st)

`"gone"` / `"dispersed"` / `"islands"` / `"single"` from `island_stats`.

Validated against ~60 runs by eye: `dispersed` is the salt-and-pepper texture,
`single` is the central blob that never breaks up, `islands` is the phenomenon.
"""
function classify_islands(st)
	st.nAB == 0        && return "gone"
	st.clustered < 0.6 && return "dispersed"
	(st.n_isl >= 2 && st.frac_max <= 0.7) && return "islands"
	return "single"
end;


# ---- notebook cell 756cfff2 ----------------------------------------------------
# --- whole-trajectory readouts -------------------------------------------------

"""
	island_timeline(d, layers; from = 1, verbose = true)

Score EVERY frame of `d` (from `load_frames`, or `getParameter`-style output).
Returns a vector of NamedTuples and prints a table.

Scoring a single late frame is how the monolayer was mis-read: the A,B phase is
a transient, so a readout keyed to a low fAB lands after the islands have both
formed and started to erode. Always look at the whole thing.
"""
function island_timeline(d, layers; from = 1, verbose = true)
	out = NamedTuple[]
	verbose && @printf("%7s %6s %6s %6s | %5s %6s %6s | %-9s %s\n",
		"t", "N", "fAB", "nAB", "n_isl", "fmax", "clust", "class", "top sizes")
	for k in from:length(d.x)
		x, y, z, s = d.x[k], d.y[k], d.z[k], d.cell_state[k]
		st  = island_stats(x, y, z, s, layers)
		cls = classify_islands(st)
		push!(out, (t = d.t[k], N = length(s), fAB = st.nAB/length(s),
					n_isl = st.n_isl, frac_max = st.frac_max,
					clustered = st.clustered, class = cls))
		verbose && @printf("%7.0f %6d %6.3f %6d | %5d %6.3f %6.3f | %-9s %s\n",
			d.t[k], length(s), st.nAB/length(s), st.nAB, st.n_isl,
			st.frac_max, st.clustered, cls, join(st.top, ","))
	end
	out
end;

"""
	island_verdict(tl)

Collapse an `island_timeline` into one word.

`"ISLANDS"` means islands are present at the END of the run, not merely at some
point: at low adhesion the monolayer fragments and then dissolves again, which
scores `"transient"` and is NOT the phenomenon.
"""
function island_verdict(tl)
	isempty(tl) && return "empty"
	nis = count(r -> r.class == "islands", tl)
	tl[end].class == "islands" ? "ISLANDS" :
		(nis >= 3 ? "transient" : uppercase(tl[end].class))
end;

"""
	polar_stats(x, y, z, s)

Unconfined polarisation. Returns `(offset, r_AB, r_C, fAB, R)` where lengths are
in units of the aggregate radius `R`.

* `offset` distance of the A,B centroid from the aggregate centre.
  ~0 = A,B still a central core; ->1 = fully polarised to one pole.
* `r_AB`, `r_C` mean radial position of each population; `r_AB > r_C` means the
  A,B cells have moved out to the surface.
"""
function polar_stats(x, y, z, s)
	cx, cy, cz = mean(x), mean(y), mean(z)
	rho = sqrt.((x.-cx).^2 .+ (y.-cy).^2 .+ (z.-cz).^2)
	R  = quantile(rho, 0.98)
	ab = findall(v -> v == 1 || v == 2, s)
	cc = findall(==(3), s)
	isempty(ab) && return (offset = NaN, r_AB = NaN, r_C = NaN, fAB = 0.0, R = R)
	dx = mean(x[ab]) - cx; dy = mean(y[ab]) - cy; dz = mean(z[ab]) - cz
	(offset = sqrt(dx^2+dy^2+dz^2)/R,
	 r_AB = mean(rho[ab])/R,
	 r_C  = isempty(cc) ? NaN : mean(rho[cc])/R,
	 fAB  = length(ab)/length(s), R = R)
end;

"""
	confinement_report(com_or_z, L; r = 1.0)

Are the plates actually confining? Cell centres should sit in `[r, L-r]`.

This is worth running on every confined simulation. With the old constant-force
plates at layers = 1, 93% of cells sat outside that band and the "monolayer" was
1.95 cells thick -- silently, because nothing else in the pipeline looks at z.
"""
function confinement_report(z::AbstractVector, L; r = 1.0, label = "")
	band = L - 2r
	sp   = maximum(z) - minimum(z)
	out  = count(v -> v < r - 1e-9 || v > L - r + 1e-9, z) / length(z)
	@printf("%s z-spread %.3f (band %.3f)  outside %.1f%%  effective layers %.2f\n",
			label, sp, band, 100out, (sp + 2r)/(2.1r))
	(spread = sp, band = band, outside = out, layers_eff = (sp + 2r)/(2.1r))
end;
confinement_report(com, L; kw...) = confinement_report(com.z, L; kw...);

"""
	msd_2d(d; from, lags = MSD_LAGS, every = 2, rim = 0.8, drift = :affine, rloc = 5.0)

In-plane mean squared displacement of the cells of a run saved with ids (save_frames,
November 2026), for the lag times `lags` (in frames; one frame per time unit in
msd_monolayer.jl). A cell is followed from frame k to k + lag only if its id is in
both, i.e. it did not divide in between. Start frames k run from `from` (the first
frame after seeding, where the ids are fresh) every `every` frames.

	drift  the growth flow (the slab spreads as cells divide) is not diffusion and is
	       removed from every displacement: :affine fits u = A r + b over the followed
	       cells (least squares) and keeps the rest; :local subtracts the mean
	       displacement of the followed cells within `rloc` (relative motion);
	       :none keeps it
	rim    only cells closer to the centre than rim x R at frame k count (R the
	       98th-percentile radius): the free edge moves differently

The cell type is the one at frame k. Returns per lag: tau (time units), msd_B,
msd_C, msd (all), n_B, n_C (displacements averaged).
"""
const MSD_LAGS = [1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 60, 75, 90, 105, 120, 150]

function msd_2d(d; from, lags = MSD_LAGS, every = 2, rim = 0.8, drift = :affine, rloc = 5.0)
	K = length(d.x)
	pos = [Dict(id => i for (i, id) in enumerate(d.id[k])) for k in 1:K]
	sB = zeros(length(lags)); sC = zeros(length(lags)); nB = zeros(Int, length(lags)); nC = zeros(Int, length(lags))
	for (il, lag) in enumerate(lags), k in from:every:(K - lag)
		x0, y0, s0 = d.x[k], d.y[k], d.cell_state[k]
		x1, y1 = d.x[k + lag], d.y[k + lag]
		P1 = pos[k + lag]
		i0 = Int[]; i1 = Int[]
		for (i, id) in enumerate(d.id[k])
			j = get(P1, id, 0)
			j > 0 && (push!(i0, i); push!(i1, j))
		end
		length(i0) < 10 && continue
		ux = x1[i1] .- x0[i0]; uy = y1[i1] .- y0[i0]
		xc, yc = mean(x0), mean(y0)
		px = x0[i0] .- xc; py = y0[i0] .- yc
		if drift == :affine
			M = hcat(px, py, ones(length(px)))
			ux = ux .- M * (M \ ux); uy = uy .- M * (M \ uy)
		elseif drift == :local
			cellk(a) = (floor(Int, px[a] / rloc), floor(Int, py[a] / rloc))
			bk = Dict{Tuple{Int,Int},Vector{Int}}()
			for a in eachindex(px); push!(get!(bk, cellk(a), Int[]), a); end
			bx = similar(ux); by = similar(uy)
			for a in eachindex(px)
				ca = cellk(a); sx = 0.0; sy = 0.0; m = 0
				for ddx in -1:1, ddy in -1:1, b in get(bk, (ca[1] + ddx, ca[2] + ddy), Int[])
					if (px[b] - px[a])^2 + (py[b] - py[a])^2 < rloc^2
						sx += ux[b]; sy += uy[b]; m += 1
					end
				end
				bx[a] = ux[a] - sx / m; by[a] = uy[a] - sy / m
			end
			ux, uy = bx, by
		end
		rho = sqrt.((x0 .- xc).^2 .+ (y0 .- yc).^2)
		R = quantile(rho, 0.98)
		for a in eachindex(i0)
			rho[i0[a]] < rim * R || continue
			u2 = ux[a]^2 + uy[a]^2
			if s0[i0[a]] == 2
				sB[il] += u2; nB[il] += 1
			else
				sC[il] += u2; nC[il] += 1
			end
		end
	end
	dt = d.t[from + 1] - d.t[from]
	(tau = lags .* dt, msd_B = sB ./ max.(nB, 1), msd_C = sC ./ max.(nC, 1),
	 msd = (sB .+ sC) ./ max.(nB .+ nC, 1), n_B = nB, n_C = nC)
end;

"""
	fit_diffusion(tau, msd; tmin = 30, tmax = 120)

D from MSD = 4 D tau + c, a straight line through the lags in [tmin, tmax] (2D, so 4D).
Also the log-log slope over the same lags (1 = diffusive, 2 = ballistic, < 1 caged).
"""
function fit_diffusion(tau, msd; tmin = 30, tmax = 120)
	w = findall(t -> tmin <= t <= tmax, tau)
	length(w) < 2 && return (D = NaN, c = NaN, slope = NaN)
	A = hcat(tau[w], ones(length(w)))
	a, c = A \ msd[w]
	L = hcat(log.(tau[w]), ones(length(w)))
	sl = (L \ log.(msd[w]))[1]
	(D = a / 4, c = c, slope = sl)
end;

# ---- notebook cell 7c1d9e42 ----------------------------------------------------
# --- slab thickness over time ----------------------------------------------------

"""
	thickness_timeline(d, L; from = 1, r = 1.0)

How thick the slab stays, frame by frame. For a plate gap L the cell centres
should sit in [r, L - r]. Returns, per frame:

* `spread`      full z-spread of the centres
* `spread98`    z-spread of the central 98 % (robust to a few stray cells)
* `layers_eff`  effective number of layers, (spread98 + 2r) / (2.1 r): the gap
                for `layers` layers is 2.1 r * layers, so 1.0 = one cell thick
* `outside`     fraction of centres outside [r, L - r]
"""
function thickness_timeline(d, L; from = 1, r = 1.0)
	out = NamedTuple[]
	for k in from:length(d.z)
		z = d.z[k]
		sp   = maximum(z) - minimum(z)
		sp98 = quantile(z, 0.99) - quantile(z, 0.01)
		push!(out, (t = d.t[k], spread = sp, spread98 = sp98,
					layers_eff = (sp98 + 2r) / (2.1r),
					outside = count(v -> v < r - 1e-6 || v > L - r + 1e-6, z) / length(z)))
	end
	out
end;

# ---- notebook cell 2d4f6a81 ----------------------------------------------------
# --- pattern readout for the B/C model --------------------------------------------
#
# Mirrors how the experimental images are segmented: B cells sitting in a dense B
# neighbourhood form clusters (the coloured islands in experiments/pattern-examples),
# the others are isolated cells (the white points there).
#
#   core B cell   at least `kmin` other B cells within `cut` (1.3 cell diameters)
#   cluster       core cells linked within `cut`, plus the B cells touching them
#   island        a cluster of at least `nmin` cells
#   isolated      any other B cell
#
# The July readout (island_stats, cell above) used plain single linkage and called
# a loose texture "dispersed" -- the very kind of pattern we now want -- so the B/C
# readout is built on this instead.

"""
	b_clusters(x, y, z, s; cut = 2.6, kmin = 3)

Cluster label of every B cell (0 = isolated). Returns (idx, label, nbrs): `idx`
are the indices of the B cells, `label[i]` the cluster of idx[i], `nbrs[i]` its
B neighbours (positions in idx).
"""
function b_clusters(x, y, z, s; cut = 2.6, kmin = 3)
	idx = findall(==(2), s)
	n = length(idx)
	label = zeros(Int, n)
	nbrs = [Int[] for _ in 1:n]
	n == 0 && return (idx = idx, label = label, nbrs = nbrs)
	px, py, pz = x[idx], y[idx], z[idx]
	lo = (minimum(px), minimum(py), minimum(pz))
	cell(i) = (floor(Int, (px[i]-lo[1])/cut), floor(Int, (py[i]-lo[2])/cut),
			   floor(Int, (pz[i]-lo[3])/cut))
	bucket = Dict{NTuple{3,Int},Vector{Int}}()
	for i in 1:n; push!(get!(bucket, cell(i), Int[]), i); end
	c2 = cut^2
	for i in 1:n
		ci = cell(i)
		for dx in -1:1, dy in -1:1, dz in -1:1
			b = get(bucket, (ci[1]+dx, ci[2]+dy, ci[3]+dz), nothing)
			b === nothing && continue
			for j in b
				j == i && continue
				if (px[i]-px[j])^2 + (py[i]-py[j])^2 + (pz[i]-pz[j])^2 < c2
					push!(nbrs[i], j)
				end
			end
		end
	end
	core = [length(nbrs[i]) >= kmin for i in 1:n]
	nlab = 0
	for i in 1:n
		(core[i] && label[i] == 0) || continue
		nlab += 1
		stack = [i]; label[i] = nlab
		while !isempty(stack)
			a = pop!(stack)
			for b in nbrs[a]
				if core[b] && label[b] == 0
					label[b] = nlab; push!(stack, b)
				end
			end
		end
	end
	for i in 1:n
		core[i] && continue
		for b in nbrs[i]
			if core[b]; label[i] = label[b]; break; end
		end
	end
	(idx = idx, label = label, nbrs = nbrs)
end;

"Mean, over the B cells `idx`, of the fraction of B among each cell's neighbours within `cut`."
function looseness(x, y, z, s, idx; cut = 2.6)
	n = length(x)
	lo = (minimum(x), minimum(y), minimum(z))
	cell(i) = (floor(Int, (x[i]-lo[1])/cut), floor(Int, (y[i]-lo[2])/cut),
			   floor(Int, (z[i]-lo[3])/cut))
	bucket = Dict{NTuple{3,Int},Vector{Int}}()
	for i in 1:n; push!(get!(bucket, cell(i), Int[]), i); end
	c2 = cut^2
	acc = 0.0; m = 0
	for i in idx
		ci = cell(i); nb = 0; nbB = 0
		for dx in -1:1, dy in -1:1, dz in -1:1
			b = get(bucket, (ci[1]+dx, ci[2]+dy, ci[3]+dz), nothing)
			b === nothing && continue
			for j in b
				j == i && continue
				if (x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2 < c2
					nb += 1
					s[j] == 2 && (nbB += 1)
				end
			end
		end
		if nb > 0; acc += nbB / nb; m += 1; end
	end
	m == 0 ? NaN : acc / m
end;

"""
	pattern_stats(x, y, z, s; confined = true, cut = 2.6, kmin = 3, nmin = 10)

One frame of the B/C model. Returns:

	nB, fB      number and fraction of B cells
	n_isl       number of islands (clusters of >= nmin cells)
	f_isl       fraction of B cells in islands
	f_iso       fraction of B cells that are isolated
	f_max       fraction of B cells in the largest island
	offset      3D shift of the B centroid from the aggregate centre, / R
	            (a polarised aggregate can shift in-plane or towards a plate)
	offset_max  the same for the largest island
	spread      mean distance of the B cells from the centre, / R: how far B is
	            spread out (a compact central blob ~0.1, rings and islands ~0.4-0.6)
	fill        B fraction of the cells in the middle (closer to the centre than
	            0.7 x the radius a compact central blob of the current B would
	            have): ~0.7-0.9 for a central cluster, even a loose one; ~0.2-0.4
	            for a ring (C in the middle) or islands
	S           looseness: mean fraction of B among a B cell's neighbours
	            (compact islands -> ~0.8-0.9, random mixing -> ~fB)
	sizes       island sizes, largest first
	confined    as passed in (classify_pattern uses it)

R and `spread` use in-plane distances when `confined` (a slab), 3D otherwise; R is
the 98th percentile distance of all cells from the aggregate centre.
"""
function pattern_stats(x, y, z, s; confined = true, cut = 2.6, kmin = 3, nmin = 10)
	N = length(s)
	cl = b_clusters(x, y, z, s; cut = cut, kmin = kmin)
	idx, lab = cl.idx, cl.label
	nB = length(idx)
	xc, yc, zc = mean(x), mean(y), mean(z)
	dz = confined ? 0.0 : 1.0
	rho = sqrt.((x .- xc).^2 .+ (y .- yc).^2 .+ dz .* (z .- zc).^2)
	R = quantile(rho, 0.98)
	nB == 0 && return (nB = 0, fB = 0.0, n_isl = 0, f_isl = 0.0, f_iso = 0.0, f_max = 0.0,
					   offset = NaN, offset_max = NaN, spread = NaN, fill = NaN, S = NaN, sizes = Int[],
					   confined = confined)
	counts = Dict{Int,Int}()
	for l in lab; l > 0 && (counts[l] = get(counts, l, 0) + 1); end
	isl = [l for (l, c) in counts if c >= nmin]
	sizes = sort([counts[l] for l in isl]; rev = true)
	shift(sel) = sqrt((mean(x[sel]) - xc)^2 + (mean(y[sel]) - yc)^2 + (mean(z[sel]) - zc)^2) / R
	offset = shift(idx)
	offset_max = isempty(isl) ? NaN : shift(idx[lab .== isl[argmax([counts[l] for l in isl])]])
	inner = findall(<(0.7 * (confined ? sqrt(nB / N) : cbrt(nB / N)) * R), rho)
	fill = isempty(inner) ? NaN : count(i -> s[i] == 2, inner) / length(inner)
	(nB = nB, fB = nB / N, n_isl = length(isl), f_isl = sum(sizes; init = 0) / nB,
	 f_iso = count(==(0), lab) / nB, f_max = isempty(sizes) ? 0.0 : sizes[1] / nB,
	 offset = offset, offset_max = offset_max, spread = mean(rho[idx]) / R, fill = fill,
	 S = looseness(x, y, z, s, idx; cut = cut), sizes = sizes, confined = confined)
end;

"""
	classify_pattern(st; pol = 0.4, cen = nothing, iso = 0.7)

From pattern_stats, in this order:

	"NONE"       no B left, or most B cells isolated (f_iso >= iso)
	"POLARISED"  B as a whole has moved to one side (offset >= pol)
	"CENTRAL"    B gathered in the middle: spread < cen, or a cluster holding
	             half of B that fills the middle (fill >= 0.55)
	"RING"       B spread out, but still one connected network centred on the
	             middle (largest island >= half of B, its centre within 0.2 R)
	"ISLANDS"    B spread out in several islands
	"CENTRAL"    also: a single loose cluster that has not moved (at most one island,
	             offset < 0.25) -- too little B left for a pattern; typical of very
	             small aggregates (Oct 7)
	"OTHER"      anything else -- look at it (e.g. a cluster drifting off-centre)

`cen` defaults to 0.30 in a slab and 0.45 in 3D: B gathered into a central blob
holding 20 % of the cells has spread ~0.30 in a disc and ~0.44 in a ball -- so
at the readout a slightly loose central cluster misses it; the `fill` test
catches those (added Oct 7, checked by eye: central clusters 0.65-0.85, rings
0.24-0.41; reports/REPORT.md, stage 5).
"""
function classify_pattern(st; pol = 0.4, cen = nothing, iso = 0.7)
	(st.nB == 0 || st.f_iso >= iso) && return "NONE"
	st.offset >= pol && return "POLARISED"
	c = cen === nothing ? (st.confined ? 0.30 : 0.45) : cen
	st.spread < c && return "CENTRAL"
	(st.f_max >= 0.5 && st.fill >= 0.55) && return "CENTRAL"
	(st.f_max >= 0.5 && st.offset_max < 0.2) && return "RING"
	st.n_isl >= 2 && return "ISLANDS"
	(st.n_isl <= 1 && st.offset < 0.25) && return "CENTRAL"
	return "OTHER"
end;

"""
	pattern_timeline(d; from = 1, confined = true, kw...)

pattern_stats and classify_pattern for every frame of `d` (from load_frames),
starting at frame `from`. Each row also carries t and the frame index k.
"""
function pattern_timeline(d; from = 1, confined = true, kw...)
	out = NamedTuple[]
	for k in from:length(d.x)
		st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = confined, kw...)
		push!(out, merge((k = k, t = d.t[k], class = classify_pattern(st)), st))
	end
	out
end;

"""
	readout_row(tl; fB_read = 0.2)

The row of a pattern_timeline at the readout: the first frame where the B fraction
has dropped to `fB_read` (the experimental 96 h reference), or the last frame if
it never does.
"""
function readout_row(tl; fB_read = 0.2)
	i = findfirst(r -> r.fB <= fB_read, tl)
	tl[i === nothing ? length(tl) : i]
end;
# --- experiment-like readout (November 2026) ------------------------------------------
#
# How the experimentalists quantify their 3D images (reports/2026-11 november, "The
# experiment"), in model units: r = R0 = 5 um, so their 13 um is 2.6 r.
#
#   core cell   a B cell with at least `frac` (40 %) B among all its neighbours
#               within `rad` (13 um)
#   cluster     core cells linked within `rad`, directly or through other core cells
#   border      a non-core B cell joins the cluster of its nearest core cell, if
#               that core cell is within `rad`
#   kept        clusters of at least `nmin` (20) cells, core and border together
#   N_eff       1 / sum p_i^2, p_i = share of the clustered cells in cluster i
#   P           |r - r_B| / R: shift of the B centre of mass from the aggregate's,
#               over half the length of the aggregate along its longest axis (PC1);
#               polarised if P >= 0.3

using LinearAlgebra: eigen, Symmetric

"For every cell, the indices of the other cells within `rad` (cell-linked buckets)."
function neighbour_lists(x, y, z; rad = 2.6)
	n = length(x)
	lo = (minimum(x), minimum(y), minimum(z))
	cell(i) = (floor(Int, (x[i]-lo[1])/rad), floor(Int, (y[i]-lo[2])/rad),
			   floor(Int, (z[i]-lo[3])/rad))
	bucket = Dict{NTuple{3,Int},Vector{Int}}()
	for i in 1:n; push!(get!(bucket, cell(i), Int[]), i); end
	r2 = rad^2
	nb = [Int[] for _ in 1:n]
	for i in 1:n
		ci = cell(i)
		for dx in -1:1, dy in -1:1, dz in -1:1
			b = get(bucket, (ci[1]+dx, ci[2]+dy, ci[3]+dz), nothing)
			b === nothing && continue
			for j in b
				if j != i && (x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2 < r2
					push!(nb[i], j)
				end
			end
		end
	end
	nb
end;

"""
	exp_clusters(x, y, z, s; rad = 2.6, frac = 0.4, nmin = 20)

The experiments' DBSCAN-like clustering of the B cells (see above). Returns
(n, f_clustered, sizes, N_eff, label): `label` per cell (0 = in no kept cluster),
`sizes` largest first, `f_clustered` the share of the B cells in kept clusters.
"""
function exp_clusters(x, y, z, s; rad = 2.6, frac = 0.4, nmin = 20)
	n = length(x)
	nb = neighbour_lists(x, y, z; rad = rad)
	isB = s .== 2
	core = [isB[i] && !isempty(nb[i]) && count(j -> isB[j], nb[i]) >= frac * length(nb[i]) for i in 1:n]
	lab = zeros(Int, n); nl = 0
	for i in 1:n
		(core[i] && lab[i] == 0) || continue
		nl += 1; lab[i] = nl; stack = [i]
		while !isempty(stack)
			a = pop!(stack)
			for b in nb[a]
				if core[b] && lab[b] == 0
					lab[b] = nl; push!(stack, b)
				end
			end
		end
	end
	for i in 1:n                     # border cells: cluster of the nearest core cell
		(isB[i] && !core[i]) || continue
		best = 0; dbest = Inf
		for j in nb[i]
			core[j] || continue
			d2 = (x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2
			if d2 < dbest
				dbest = d2; best = j
			end
		end
		best > 0 && (lab[i] = lab[best])
	end
	counts = zeros(Int, nl)
	for l in lab; l > 0 && (counts[l] += 1); end
	keep = findall(>=(nmin), counts)
	newl = zeros(Int, nl); newl[keep] = 1:length(keep)
	lab = [l > 0 ? newl[l] : 0 for l in lab]
	sizes = sort(counts[keep]; rev = true)
	tot = sum(sizes; init = 0); nB = count(isB)
	(n = length(sizes), f_clustered = nB == 0 ? 0.0 : tot / nB, sizes = sizes,
	 N_eff = tot == 0 ? 0.0 : 1 / sum((sizes ./ tot) .^ 2), label = lab)
end;

"""
	exp_polarisation(x, y, z, s)

The experiments' polarisation: P = |r - r_B| / R, with r and r_B the centres of mass
of all cells and of the B cells, and R half the extent of all cells along their first
principal axis (PC1). Returns (P, R, shift).
"""
function exp_polarisation(x, y, z, s)
	b = findall(==(2), s)
	isempty(b) && return (P = NaN, R = NaN, shift = NaN)
	X = hcat(x .- mean(x), y .- mean(y), z .- mean(z))
	v = eigen(Symmetric(X' * X ./ length(x))).vectors[:, end]     # PC1
	p = X * v
	R = (maximum(p) - minimum(p)) / 2
	sh = sqrt(mean(X[b, 1])^2 + mean(X[b, 2])^2 + mean(X[b, 3])^2)
	(P = sh / R, R = R, shift = sh)
end;

"""
	exp_stats(x, y, z, s; rad = 2.6, frac = 0.4, nmin = 20)

The experiment-like readout of one frame, as one NamedTuple: n_cl (kept clusters),
f_cl (share of B in them), N_eff, size_max (largest cluster) and P.
"""
function exp_stats(x, y, z, s; rad = 2.6, frac = 0.4, nmin = 20)
	cl = exp_clusters(x, y, z, s; rad = rad, frac = frac, nmin = nmin)
	po = exp_polarisation(x, y, z, s)
	(n_cl = cl.n, f_cl = cl.f_clustered, N_eff = cl.N_eff,
	 size_max = isempty(cl.sizes) ? 0 : cl.sizes[1], P = po.P)
end;

"""
	classify_exp(st; pol = 0.3, neff = 1.5)

From exp_stats, in this order:

	"NONE"       no B cluster of `nmin` cells (B scattered, or too little B left)
	"POLARISED"  P >= pol (the experiments' cutoff, 0.3)
	"ISLANDS"    N_eff >= neff: B split into clusters of comparable size
	"CENTRAL"    one dominant cluster that has not moved. An unbroken ring is
	             CENTRAL here, as in the experiments' medium confinement ("may form
	             a ring-like pattern but without breaking").
"""
function classify_exp(st; pol = 0.3, neff = 1.5)
	st.n_cl == 0 && return "NONE"
	st.P >= pol && return "POLARISED"
	st.N_eff >= neff && return "ISLANDS"
	return "CENTRAL"
end;

# ---- notebook cell 047d722c ----------------------------------------------------
# --- adhesion energetics (Differential Adhesion Hypothesis) --------------------
#
# The pair force is  F(d) = f * (rij/d - 1) * (mu*rij/d - 1),  rij = 2r,
# attractive for rij < d < mu*rij. Integrating it over the attractive branch
# gives the work to separate one bonded pair:
#
#     W_ij = f_att[i,j] * rij * A(mu_ij),   A(mu) = (1+mu)*ln(mu) - 2*(mu-1)
#
# f_rep does NOT appear: it only acts below contact, where the force is
# repulsive, so it sets compressibility, not adhesion. (It still matters for the
# patterns -- it decides whether a tissue can round up against its own
# repulsion -- but it does not enter W.)
#
# Surface tension of a pure aggregate of type i, and interfacial tension between
# i and j, up to the same geometric prefactor (coordination number / 2):
#     sigma_i   = W_ii / 2
#     gamma_ij  = (W_ii + W_jj)/2 - W_ij
#
# Steinberg's classification of a doublet, from these:
#     gamma_ij > sigma_i + sigma_j    -> the two aggregates separate
#     |sigma_i - sigma_j| > gamma_ij  -> one completely engulfs the other
#     otherwise                       -> partial engulfment, a stable doublet
#                                        with a finite contact angle: "glue"
# BC needs to land in the third case.

A_mu(mu) = (1 + mu)*log(mu) - 2*(mu - 1);

"Pair separation work matrix W[i,j] = f_att[i,j] * 2r * A(mu[i,j])."
adhesion_work(f_att, mu; r = 1.0) = f_att .* (2r) .* A_mu.(mu);

"""
	tension_table(f_att, mu; r = 1.0, labels = ["A","B","C"])

Print sigma_i, gamma_ij and the Steinberg class of every pair. Use it to check a
candidate parameter set BEFORE running any simulation.
"""
function tension_table(f_att, mu; r = 1.0, labels = ["A","B","C"])
	W = adhesion_work(f_att, mu; r = r)
	sigma = [W[i,i]/2 for i in 1:size(W,1)]
	println("W (pair separation work)");  display(round.(W, digits = 4))
	@printf("\nsigma: %s\n\n", join((@sprintf("%s=%.4f", labels[i], sigma[i]) for i in eachindex(sigma)), "  "))
	@printf("%-5s %9s %9s  %s\n", "pair", "gamma", "s_i+s_j", "class")
	for i in 1:size(W,1), j in i+1:size(W,1)
		g = (W[i,i] + W[j,j])/2 - W[i,j]
		cls = g > sigma[i] + sigma[j] ? "separate" :
			  abs(sigma[i] - sigma[j]) > g ? "engulfment" : "partial (glue)"
		@printf("%-5s %9.4f %9.4f  %s\n", labels[i]*labels[j], g, sigma[i]+sigma[j], cls)
	end
	(W = W, sigma = sigma)
end;

# The set used below, for reference:
#   tension_table(4*[0.9 1 0.6; 1.0 1.6 0.4; 0.6 0.4 0.75],
#                 2*[1 0.91 1; 0.91 1 1; 1 1 1.07])


# ---- notebook cell cad66be3 ----------------------------------------------------
# --- two-aggregate fusion assay ------------------------------------------------
#
# Time calibration: p = 0.009333 <-> 7.1 h and q = 0.004667 <-> 14.1 h both give
#     1 hour = 15 model time units,      so 12 h = 180, 24 h = 360.
const HOUR = 15.0;

"""
	make_aggregate(n, state, mu_, fatt_, frep_; seed = 99, relax = 200.0)

Grow `n` cells and relax them AS A PURE AGGREGATE OF TYPE `state` under the
assay mechanics. Returns centred `(x, y, z)`.

Relaxing under the right type's mechanics matters: growth uses `parameters_tfm`
(olap = 0.75, i.e. heavily overlapping), and dropping that straight into a stiff
`f_rep` makes the aggregates explode apart on contact.
"""
function make_aggregate(n, state, mu_, fatt_, frep_; seed = 99, relax = 200.0)
	Random.seed!(seed)
	c = Community(model, N = 1, dt = 0.002, simBox = [-80.0 80; -80 80; -80 80])
	for (par, val) in pairs(parameters_tfm); c[par] = val; end
	c.cell_state = 1
	c.x = 0.0; c.y = 0.0; c.z = 0.0
	c.vx = 0.0; c.vy = 0.0; c.vz = 0.0
	c.t_div = 1
	c.g_on = true; c.d_on = false
	c.confinement = false; c.L = 4.2; c.rep_wall = 60; c.alpha_ecm = 0.9
	c.abm.neighbors.cellEdge = (2*maximum(c.r))*maximum(c.mu)
	loadToPlatform!(c, preallocateAgents = round(Int, 1.5n))
	while c.N[1] < n
		agentStepDE!(c); agentStepRule!(c); update!(c); computeNeighbors!(c)
	end
	bringFromPlatform!(c)

	for (par, val) in pairs(parameters_do); c[par] = val; end
	setfield!(c, :dt, 0.015)
	c.mu = mu_; c.f_att = fatt_; c.f_rep = frep_
	c.cell_state .= state
	c.g_on = false; c.d_on = false
	c.vx = 0.0; c.vy = 0.0; c.vz = 0.0
	c.confinement = false
	c.alpha_ecm = 0.9 * c.lambda
	c.abm.neighbors.cellEdge = (2*maximum(c.r))*maximum(c.mu)
	loadToPlatform!(c)
	for _ in 1:round(Int, relax/0.015)
		agentStepDE!(c); update!(c); computeNeighbors!(c)
	end
	bringFromPlatform!(c)
	(copy(c.x .- mean(c.x)), copy(c.y .- mean(c.y)), copy(c.z .- mean(c.z)))
end;

"""
	fusion_assay(mu_, fatt_, frep_; n = 150, hours = 24.0, pairs_ = ...)

Reproduce the doublet-fusion figure. Two aggregates of `n` cells each are placed
just touching, with no division and no differentiation, and the metric is

    COM separation / (Rg1 + Rg2)

~1.6 at contact, ~0.5 once fully rounded into one sphere. LOWER = MORE FUSED.

The signature to look for is not an absolute number but the ORDERING: BC should
be the least-fused pair by a clear margin while the others round up, i.e. two
lobes that adhere without becoming a sphere.
"""
function fusion_assay(mu_, fatt_, frep_; n = 150, hours = 24.0,
		pairs_ = [(2,2),(3,3),(2,3)], labels = ["A","B","C"],
		verbose = true, positions = nothing)
	# only the types that appear in pairs_ (no time spent on A when it is not asked for)
	agg = Dict(st => make_aggregate(n, st, mu_, fatt_, frep_) for st in unique([p[k] for p in pairs_ for k in 1:2]))
	results = Dict{String,NTuple{3,Float64}}()

	for (i, j) in pairs_
		x1, y1, z1 = agg[i]; x2, y2, z2 = agg[j]
		n1, n2 = length(x1), length(x2)
		R1 = maximum(sqrt.(x1.^2 .+ y1.^2 .+ z1.^2))
		R2 = maximum(sqrt.(x2.^2 .+ y2.^2 .+ z2.^2))
		dsep = (R1 + R2)/2 + 1.0

		c = Community(model, N = n1+n2, dt = 0.015, simBox = [-80.0 80; -80 80; -80 80])
		for (par, val) in pairs(parameters_do); c[par] = val; end
		c.mu = mu_; c.f_att = fatt_; c.f_rep = frep_
		c.x = vcat(x1 .- dsep, x2 .+ dsep)
		c.y = vcat(y1, y2); c.z = vcat(z1, z2)
		c.vx = 0.0; c.vy = 0.0; c.vz = 0.0
		c.cell_state = vcat(fill(i, n1), fill(j, n2))
		c.t_div = 1e9
		c.g_on = false; c.d_on = false
		c.confinement = false
		c.alpha_ecm = 0.9 * c.lambda
		setfield!(c, :t, 0.0)
		# the starting configuration too (positions[<pair>_0]), for the before/after figure
		positions === nothing || (positions[labels[i]*labels[j]*"_0"] =
			(x = copy(c.x), y = copy(c.y), z = copy(c.z), r = copy(c.r), cell_state = copy(c.cell_state)))
		c.abm.neighbors.cellEdge = (2*maximum(c.r))*maximum(c.mu)
		loadToPlatform!(c)

		g1 = 1:n1; g2 = n1+1:n1+n2
		function metric()
			x, y, z = c.x, c.y, c.z
			c1 = (mean(x[g1]), mean(y[g1]), mean(z[g1]))
			c2 = (mean(x[g2]), mean(y[g2]), mean(z[g2]))
			sep = sqrt((c1[1]-c2[1])^2 + (c1[2]-c2[2])^2 + (c1[3]-c2[3])^2)
			rg(g, cc) = sqrt(mean((x[g].-cc[1]).^2 .+ (y[g].-cc[2]).^2 .+ (z[g].-cc[3]).^2))
			sep / (rg(g1, c1) + rg(g2, c2))
		end

		m0 = metric(); m12 = NaN
		k12 = round(Int, 12.0*HOUR/0.015)
		for step in 1:round(Int, hours*HOUR/0.015)
			agentStepDE!(c); update!(c); computeNeighbors!(c)
			step == k12 && (m12 = metric())
		end
		m24 = metric()
		bringFromPlatform!(c)
		results[labels[i]*labels[j]] = (m0, m12, m24)
		# final configuration, for 3D plots of the doublet (pass positions = Dict())
		positions === nothing || (positions[labels[i]*labels[j]] =
			(x = copy(c.x), y = copy(c.y), z = copy(c.z), r = copy(c.r), cell_state = copy(c.cell_state)))
	end

	if verbose
		println("\n--- fusion assay (lower = more fused) ---")
		@printf("%-4s %8s %8s %8s\n", "pair", "t0", "12h", "24h")
		for (i, j) in pairs_
			tag = labels[i]*labels[j]
			@printf("%-4s %8.3f %8.3f %8.3f\n", tag, results[tag]...)
		end
	end
	results
end;

# ---- notebook cell 8a8800de ----------------------------------------------------
# --- parameter sets, plates, and one confined run end to end ---------------------
#
# PSETS holds the three starting points for the B/C mechanics and the set chosen
# from them in stage 5 ("final", reports/REPORT.md). In all of them the
# BC entries ([2,3] and [3,2]) are the FUSION values -- mu 2, f_att 1.6, f_rep 3 --
# which stay frozen ("mainly do not touch the [2,3] and [3,2] elements"); only
# the BB and CC entries differ. Row/column 1 (A) is kept so the matrices still
# match the original ones; no A cell is ever created.
PSETS = Dict(
	# your CONFINED PATTERNS set (= 0.5 * parameters_do f_att, mu = 2), with the
	# FUSION BC entry: f_att[2,3] 0.8 -> 1.6
	"confined_patterns" => (
		mu    = 2 * [1 1 1; 1 1 1; 1 1 1],
		f_att =     [3 3 1.2; 3 3 1.6; 1.2 1.6 2.4],
		f_rep = 3 * [1 1 1; 1 1 1; 1 1 1.7]),
	# the July set: FUSION mu and f_att, BB and CC repulsion softened
	"july" => (
		mu    = 2 * [1 0.91 1; 0.91 1 1; 1 1 1.07],
		f_att = 4 * [0.9 1.0 0.6; 1.0 1.6 0.4; 0.6 0.4 0.75],
		f_rep = 3 * [3 1 1; 1 1 1; 1 1 1.7]),
	# Sotiris' FUSION set, unchanged
	"fusion" => (
		mu    = 2 * [1 0.91 1; 0.91 1 1; 1 1 1.07],
		f_att = 4 * [0.9 1.0 0.6; 1.0 1.6 0.4; 0.6 0.4 0.75],
		f_rep = 3 * [3 1 1; 1 2 1; 1 1 2]),
	# the set chosen on Oct 7: B as stiff as in july_frB4.5 (f_rep BB 4.5, the
	# looseness you picked), C as stiff as in FUSION (f_rep CC 6, keeps a small
	# squeezed B core from sliding to the rim). Patterns only in 1-2 layers.
	# In the stage 5 scripts it is "july_frB4.5_frC6@lw1.5".
	"final_oct7" => (
		mu    = 2 * [1 0.91 1; 0.91 1 1; 1 1 1.07],
		f_att = 4 * [0.9 1.0 0.6; 1.0 1.6 0.4; 0.6 0.4 0.75],
		f_rep = 3 * [3 1 1; 1 1.5 1; 1 1 2]),
	# THE CHOSEN SET (Oct 8): final_oct7 with a shorter B-B attraction range, mu BB
	# 2 -> 1.9. B is less cohesive (adhesion energy W_BB 1.02 -> 0.79, still above
	# C's 0.65), so a big B core in a slab of 3+ layers opens into a ring and
	# islands instead of staying one blob. Used with plate friction LAMBDA_WALL.
	# In the stage 5 scripts it is "july_frB4.5_frC6_muB1.9@lw1.5".
	"final" => (
		mu    = 2 * [1 0.91 1; 0.91 0.95 1; 1 1 1.07],
		f_att = 4 * [0.9 1.0 0.6; 1.0 1.6 0.4; 0.6 0.4 0.75],
		f_rep = 3 * [3 1 1; 1 1.5 1; 1 1 2]),
);

LAMBDA_WALL = 1.5;   # friction with the plates used with "final" (stage 5; MODEL.md section 4)

TAU_WALL = 0.15;     # rigid plates: an overshoot is undone in 10 steps of dt = 0.015

"""
	wall_settings(layers, wall = :rigid) -> (wall_type, rep_wall, dt)

The plate models compared in stage 2 (reports/REPORT.md):

	:constant       the original: constant force 30 while a cell overlaps a plate
	:hookean        July: spring; (150, dt 0.005) at 1 layer, (60, dt 0.015) otherwise
	:hookean_stiff  spring 150 at every thickness, dt 0.005
	:rigid          no motion out through a plate; dt 0.015

The Hookean spring is only stable for rep_wall * dt <= ~0.75 (measured in July),
hence the smaller time step whenever it is stiff.
"""
function wall_settings(layers, wall = :rigid)
	wall == :constant      && return (1, 30.0, 0.015)
	wall == :hookean       && return layers <= 1 ? (2, 150.0, 0.005) : (2, 60.0, 0.015)
	wall == :hookean_stiff && return (2, 150.0, 0.005)
	wall == :rigid         && return (3, 0.0, 0.015)
	error("unknown wall $wall")
end;

"""
	run_confined(layers, N0; pset = "july", wall = :rigid, fB = 0.6, width = 0.0, f_in = nothing, c_core = 0.0, lambda_scale = 1.0,
	             lambda_wall = 0.0, tau_div = nothing, tf = 625.0, fB_stop = nothing, seed = 2345, outdir = "coms",
	             save_every = 15.0, verbose = true)

Your Initialization + Differentiation cells as one function, so a sweep does not
depend on the notebook's global state: grow N0 cells between the plates, relax
under the B/C mechanics, seed the B core (initialize_diff!), differentiate for tf
(or until fB_stop), and save the trajectory with save_frames -- one frame every
`save_every` time units (15 = one hour).

`lambda_scale` multiplies the friction lambda (and alpha_ecm = 0.9 lambda with it) after
growth: the mechanics run slower against division and differentiation (Oct 8 test).
`width` / `f_in` dilute the B core and `c_core` mixes C cells into it (see
initialize_diff!).
`lambda_wall` > 0 adds friction with the plates (extra in-plane drag on cells
touching a plate) from the relaxation onwards; growth is always frictionless.
`tau_div` replaces the mean division times [A, B, C] of parameters_do (Nov 2026
test: [195, 195, 195] = 13 h for every cell); nothing = keep 300, 300, 750.

Returns (file, m1, timeline, verdict, secs); m1 is the frame right after seeding.
"""
function run_confined(layers, N0; pset = "july", wall = :rigid, fB = 0.6, width = 0.0, f_in = nothing, c_core = 0.0, lambda_scale = 1.0,
		lambda_wall = 0.0, tau_div = nothing, tf = 625.0, fB_stop = nothing, seed = 2345, outdir = "coms",
		save_every = 15.0, verbose = true)

	wt, rw, dtd = wall_settings(layers, wall)
	Lgap = (1.05*2.0)*layers
	P = PSETS[pset]
	core = f_in === nothing ? @sprintf("w%g", width) : @sprintf("fin%g", f_in)
	fB == 0.6 || (core *= @sprintf("_fB%g", fB))
	c_core > 0 && (core *= @sprintf("_cc%g", c_core))
	lambda_scale == 1 || (core *= @sprintf("_lam%g", lambda_scale))
	tau_div === nothing || (core *= @sprintf("_tdiv%g", tau_div[2]))
	save_every == 15 || (core *= @sprintf("_se%g", save_every))
	name = @sprintf("L%g_N%d_%s_%s_%s%s_s%d", layers, N0, pset, wall, core,
				   lambda_wall > 0 ? @sprintf("_lw%g", lambda_wall) : "", seed)

	Random.seed!(seed)
	c = Community(model, N = 1, dt = dt_tfm, simBox = [-200.0 200; -200 200; -200 200])
	for (par, val) in pairs(parameters_tfm); c[par] = val; end
	c.cell_state = 3
	c.x = 0.0; c.y = 0.0; c.z = 0.0
	c.vx = 0.0; c.vy = 0.0; c.vz = 0.0
	c.t_div = 1
	c.g_on = true; c.d_on = false
	c.confinement = true; c.L = Lgap
	c.rep_wall = rw; c.wall_type = wt; c.tau_wall = TAU_WALL; c.lambda_wall = 0.0
	c.alpha_ecm = 0.9*c.lambda[1]
	grow_size!(c, save_each_tfm, N0; show_progress = false)

	for (par, val) in pairs(parameters_do); c[par] = val; end
	lambda_scale == 1 || (c.lambda = lambda_scale * c.lambda[1])   # slower mechanics (Oct 8 test)
	tau_div === nothing || (c.tau_div = Float64.(tau_div))         # other division times (Nov test)
	setfield!(c, :dt, dtd)
	c.mu = P.mu; c.f_att = P.f_att; c.f_rep = P.f_rep
	c.confinement = true; c.L = Lgap
	c.rep_wall = rw; c.wall_type = wt; c.tau_wall = TAU_WALL; c.lambda_wall = lambda_wall
	c.alpha_ecm = 0.9*c.lambda[1]

	se = round(Int, save_every/dtd)
	stabilize!(c, se, 100; show_progress = false)
	initialize_diff!(c; g_on = true, fB = fB, width = width, f_in = f_in, confined = true, c_core = c_core)
	# fresh ids 1..N: CellBasedModels restarts its id counter at N every time a run is
	# loaded, so older ids could repeat; from here on they stay unique (msd_2d)
	setfield!(c, :id, collect(1:length(getfield(c, :id))))
	m1 = length(c)
	stabilize!(c, se, 25; show_progress = false)
	secs = @elapsed differentiate!(c, se, tf; fB_stop = fB_stop, show_progress = false)

	file = save_frames(joinpath(outdir, "com_$(name).jld2"), c)

	d  = load_frames(file)
	tl = island_timeline(d, layers; from = m1, verbose = false)
	v  = island_verdict(tl)
	if verbose
		e = tl[end]
		@printf("%-40s -> %-9s  n_isl=%d fmax=%.2f clust=%.2f fB=%.3f  (%.0fs)  ",
				name, v, e.n_isl, e.frac_max, e.clustered, e.fAB, secs)
		confinement_report(d.z[end], Lgap)
	end
	(file = file, m1 = m1, timeline = tl, verdict = v, secs = secs)
end;

"""
	phase_scan(layers_v, N0_v; kw...)

`run_confined` over a grid of layers and N0, tabulated (rough version; the
phase space of stage 5 replaces it).
"""
function phase_scan(layers_v, N0_v; kw...)
	rows = NamedTuple[]
	for L in layers_v, N0 in N0_v
		r = run_confined(L, N0; kw...)
		e = r.timeline[end]
		push!(rows, (layers = L, N0 = N0, verdict = r.verdict,
					 n_isl = e.n_isl, frac_max = e.frac_max,
					 clustered = e.clustered, file = r.file))
	end
	println()
	@printf("%7s %7s  %-10s %6s %8s %8s\n", "layers", "N0", "verdict", "n_isl", "frac_max", "clust")
	for r in rows
		@printf("%7.1f %7d  %-10s %6d %8.3f %8.3f\n",
				r.layers, r.N0, r.verdict, r.n_isl, r.frac_max, r.clustered)
	end
	rows
end;

# ---- notebook cell 5e7c1a90 ----------------------------------------------------
# --- one unconfined run end to end: the polarisation check -------------------------

"""
	run_unconfined(N0 = 500; pset = "july", fB = 0.6, width = 0.0, f_in = nothing, c_core = 0.0, lambda_scale = 1.0,
	               tau_div = nothing, tf = 700.0, fB_stop = nothing, seed = 2345, outdir = "coms",
	               save_every = 15.0, show_progress = false, verbose = true)

Grow a free aggregate, seed a SPHERICAL B core (initialize_diff! with
confined = false), differentiate, save, and print how far the B centroid moves
from the centre (polar_stats: offset ~0 = central core, -> 1 = at one pole).

`lambda_scale` multiplies the friction lambda (and alpha_ecm = 0.9 lambda with it) after
growth: the mechanics run slower against division and differentiation (Oct 8 test).
`width` / `f_in` dilute the B core, `c_core` mixes C into it (see initialize_diff!);
`tau_div` replaces the mean division times [A, B, C] (as in run_confined);
verbose = false skips the table.

show_progress = false by default: differentiate!'s progress display calls
IJulia.clear_output and plot_pancake, which only exist inside Jupyter.
"""
function run_unconfined(N0 = 500; pset = "july", fB = 0.6, width = 0.0, f_in = nothing, c_core = 0.0, lambda_scale = 1.0,
		tau_div = nothing, tf = 700.0, fB_stop = nothing, seed = 2345, outdir = "coms", save_every = 15.0,
		show_progress = false, verbose = true)

	P = PSETS[pset]
	core = f_in === nothing ? @sprintf("w%g", width) : @sprintf("fin%g", f_in)
	fB == 0.6 || (core *= @sprintf("_fB%g", fB))
	c_core > 0 && (core *= @sprintf("_cc%g", c_core))
	lambda_scale == 1 || (core *= @sprintf("_lam%g", lambda_scale))
	tau_div === nothing || (core *= @sprintf("_tdiv%g", tau_div[2]))
	name = @sprintf("free_N%d_%s_%s_s%d", N0, pset, core, seed)

	Random.seed!(seed)
	c = Community(model, N = 1, dt = dt_tfm, simBox = [-150.0 150; -150 150; -150 150])
	for (par, val) in pairs(parameters_tfm); c[par] = val; end
	c.cell_state = 3
	c.x = 0.0; c.y = 0.0; c.z = 0.0
	c.vx = 0.0; c.vy = 0.0; c.vz = 0.0
	c.t_div = 1
	c.g_on = true; c.d_on = false
	c.confinement = false; c.L = 4.2; c.rep_wall = 60; c.wall_type = 2; c.tau_wall = TAU_WALL
	c.alpha_ecm = 0.9*c.lambda[1]
	grow_size!(c, save_each_tfm, N0; show_progress = show_progress)

	for (par, val) in pairs(parameters_do); c[par] = val; end
	lambda_scale == 1 || (c.lambda = lambda_scale * c.lambda[1])   # slower mechanics (Oct 8 test)
	tau_div === nothing || (c.tau_div = Float64.(tau_div))         # other division times (Nov test)
	setfield!(c, :dt, 0.015)
	c.mu = P.mu; c.f_att = P.f_att; c.f_rep = P.f_rep
	c.confinement = false
	c.alpha_ecm = 0.9*c.lambda[1]

	se = round(Int, save_every/0.015)
	stabilize!(c, se, 100; show_progress = show_progress)
	initialize_diff!(c; g_on = true, fB = fB, width = width, f_in = f_in, confined = false, c_core = c_core)
	# fresh ids 1..N, unique from here on (see run_confined)
	setfield!(c, :id, collect(1:length(getfield(c, :id))))
	m1 = length(c)
	secs = @elapsed differentiate!(c, se, tf; fB_stop = fB_stop, show_progress = show_progress)

	file = save_frames(joinpath(outdir, "com_$(name).jld2"), c)
	d = load_frames(file)
	if verbose
		@printf("%7s %6s %7s %7s %7s\n", "t", "fB", "offset", "r_B", "r_C")
		for k in m1:length(d.x)
			st = polar_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k])
			@printf("%7.0f %6.3f %7.3f %7.3f %7.3f\n", d.t[k], st.fAB, st.offset, st.r_AB, st.r_C)
		end
	end
	(file = file, frames = d, m1 = m1, secs = secs)
end;

# ---- notebook cell ff3415f8 ----------------------------------------------------
# --- phase diagram ------------------------------------------------------------

"""
	plot_phase(rows; filename = "phase_diagram", savefig = false)

Plot the output of `phase_scan`: layers vs N0, one marker per run.

Two categories are drawn:

    patterns          A,B islands present at the END of the run
    central cluster   the A,B mass stays as one cluster

`transient` runs (islands form and then dissolve) are drawn as *central
cluster*, which is what they actually end as -- the one such run in the map,
layers = 1 at N0 = 200, finishes with a single cluster holding 67% of the A,B
cells. See `island_verdict` if you need to tell the two apart.

A `DISPERSED` ("salt and pepper") run would be neither, so it is drawn in grey
and warned about rather than silently folded into one of the two categories.
There are none in the current map.
"""
function plot_phase(rows; filename = "phase_diagram", savefig = false)
	ispattern(v) = v == "ISLANDS"
	iscluster(v) = v in ("SINGLE", "transient")

	odd = unique(r.verdict for r in rows if !ispattern(r.verdict) && !iscluster(r.verdict))
	isempty(odd) || @warn "verdicts outside the two plotted categories, drawn grey" odd

	fig = Figure(size = (720, 420))
	ax = Axis(fig[1,1], xlabel = "N0 at squeezing", ylabel = "layers",
			  xscale = log10, yticks = sort(unique(r.layers for r in rows)),
			  title = "pattern phase diagram")
	for r in rows
		col, mrk = ispattern(r.verdict) ? (:seagreen, :circle) :
				   iscluster(r.verdict) ? (:steelblue, :rect) : (:gray, :xcross)
		scatter!(ax, [Float64(r.N0)], [Float64(r.layers)];
				 color = col, marker = mrk, markersize = 18, strokewidth = 0.5)
	end
	lv = sort(unique(r.layers for r in rows))
	ylims!(ax, minimum(lv) - 0.5, maximum(lv) + 0.5)
	Legend(fig[1,2],
		[MarkerElement(color = :seagreen, marker = :circle, markersize = 14),
		 MarkerElement(color = :steelblue, marker = :rect,   markersize = 14)],
		["patterns", "central cluster"])
	savefig && save("$(filename).png", fig)
	fig
end;

"""
	critical_N0(rows)

Bracket the pattern threshold per layer count: the largest tested N0 that did
NOT give persistent islands, and the smallest that did.

Experimentally the transition is sharp -- once patterns appear at some N0 they
appear for every larger one -- so a bracket is the whole answer.
"""
function critical_N0(rows)
	@printf("%7s %12s %12s\n", "layers", "largest no", "smallest yes")
	for L in sort(unique(r.layers for r in rows))
		sub  = filter(r -> r.layers == L, rows)
		yes  = [r.N0 for r in sub if r.verdict == "ISLANDS"]
		no   = [r.N0 for r in sub if r.verdict != "ISLANDS"]
		ymin = isempty(yes) ? nothing : minimum(yes)
		below = ymin === nothing ? no : filter(n -> n < ymin, no)
		nmax = isempty(below) ? nothing : maximum(below)
		@printf("%7.1f %12s %12s\n", L,
				nmax === nothing ? "-" : string(nmax),
				ymin === nothing ? "> $(isempty(no) ? "?" : maximum(no))" : string(ymin))
	end
end;

# ---- notebook cell 99c84e68 ----------------------------------------------------
model_new = ABM(3,

	# Global parameters
	model = Dict(
		# Mechanics
		:range => Float64,
		:mu => Array{Float64},
		:lambda => Float64,
		:f_rep => Array{Float64},
		:f_att => Array{Float64},
		:alpha_ecm => Float64,
		# Division
		:tau_div => Array{Float64},
		:sigma_div => Float64,
		:olap => Float64,
		:g_on => Bool,
		# Differentiation: B -> C only. There are no A cells, so the A -> B rate p,
		# the feedback strength k and the A fraction it needed are gone.
		:d_on => Bool,
		:q => Float64,
		# Reference values
		:t0 => Float64,
		:r0 => Float64,
		:f0 => Float64,
		# Confinement
		:confinement => Bool,
		:L => Float64,
		:rep_wall => Float64,
		:wall_type => Int64,		# 1 constant force, 2 Hookean, 3 rigid -- see wall_settings
		:tau_wall => Float64,		# rigid plates only: time to undo an overshoot
		:lambda_wall => Float64,	# extra in-plane drag on cells touching a plate (0 = frictionless)
	),


	# Agent parameters
	agent = Dict(
		:t_div => Float64,
		:ni => Int64,
		:cell_state => Int64,		# 2 = B (T/Bra+), 3 = C (T/Bra-); 1 (A) is never created
		:r => Float64,
		# Mechanics
		:vx => Float64,
		:vy => Float64,
		:vz => Float64,
		:fx => Float64,
		:fy => Float64,
		:fz => Float64,
	),


	# Mechanics
	agentODE = quote

		fx = 0
		fy = 0
		fz = 0
		vsum_x = 0
		vsum_y = 0
		vsum_z = 0
		ni = 0
		state_i = cell_state
		rij = 2 * r

		@loopOverNeighbors it2 begin
			state_j = cell_state[it2]

			dij = CBMMetrics.euclidean(x, x[it2], y, y[it2], z, z[it2])
			mu_ij = mu[state_i, state_j]

			# Cellular forces
			if dij <= mu_ij*rij && dij > 0
				if dij < rij
					f_factor = f_rep[state_i, state_j]
				else
					f_factor = f_att[state_i, state_j]
				end
				inv_dij = 1 / dij
				ratio = rij * inv_dij
				f0_ij = (ratio - 1) * (mu_ij * ratio - 1) * inv_dij
				fx += f_factor * f0_ij * (x - x[it2])
				fy += f_factor * f0_ij * (y - y[it2])
				fz += f_factor * f0_ij * (z - z[it2])
			end

			# Counting neighbours
			if dij < range * rij
				ni += 1
				vsum_x += vx[it2]
				vsum_y += vy[it2]
				vsum_z += vz[it2]
			end
		end

		# Confinement between two plates at z = 0 and z = L. Three models, compared
		# in reports/REPORT.md (stage 2):
		#   wall_type 1  constant force while the cell overlaps a plate (the original)
		#   wall_type 2  Hookean: force grows with the overlap (July)
		#   wall_type 3  rigid: no motion out through a plate (applied to vz below)
		if confinement
			if wall_type == 1
				if z < r
					fz += rep_wall
				end
				diz_top = L - z
				if diz_top < r
					fz -= rep_wall
				end
			elseif wall_type == 2
				if z < r
					fz += rep_wall * (r - z) / r
				end
				diz_top = L - z
				if diz_top < r
					fz -= rep_wall * (r - diz_top) / r
				end
			end
		end

		# Plate friction (optional, off when lambda_wall = 0): a cell touching a
		# plate feels an extra drag lambda_wall on its in-plane motion.
		wall_drag = 0.0
		if confinement && (z < r + 0.05 || z > L - r - 0.05)
			wall_drag = lambda_wall
		end

		# Equations of motion
		ni_relv_min = 2
		if ni < ni_relv_min
			inv = 1 / lambda
			inv_xy = 1 / (lambda + wall_drag)
			vx = inv_xy * fx
			vy = inv_xy * fy
			vz = inv * fz
		else
			inv = 1 / (lambda * ni)
			inv_xy = 1 / (lambda * ni + wall_drag)
			vx = inv_xy * (fx + vsum_x*alpha_ecm)
			vy = inv_xy * (fy + vsum_y*alpha_ecm)
			vz = inv * (fz + vsum_z*alpha_ecm)
		end

		# Rigid plates: a cell touching a plate cannot move further into it, and a
		# small overshoot (at most one time step) is pulled back within tau_wall.
		if confinement && wall_type == 3
			if z < r
				vz = max(vz, (r - z) / tau_wall)
			end
			if z > L - r
				vz = min(vz, (L - r - z) / tau_wall)
			end
		end

		dt(x) = vx
		dt(y) = vy
		dt(z) = vz

	end,


	# Growth and differentiation
	agentRule = quote
		# Growth
		if g_on
			if t > t_div
				x_div = CBMDistributions.normal(0, 1)
				y_div = CBMDistributions.normal(0, 1)
				z_div = CBMDistributions.normal(0, 1)
				norm_div = sqrt(x_div^2 + y_div^2 + z_div^2)
				x_div /= norm_div
				y_div /= norm_div
				z_div /= norm_div

				r_sep = r * olap
				@addAgent(
					x = x + r_sep * x_div,
					y = y + r_sep * y_div,
					z = z + r_sep * z_div,
					vx = vx / 2,
					vy = vy / 2,
					vz = vz / 2,
					t_div = t + CBMDistributions.uniform(tau_div[cell_state] * (1 - sigma_div), tau_div[cell_state] * (1 + sigma_div))
				)
				@addAgent(
					x = x - r_sep * x_div,
					y = y - r_sep * y_div,
					z = z - r_sep * z_div,
					vx = vx / 2,
					vy = vy / 2,
					vz = vz / 2,
					t_div = t + CBMDistributions.uniform(tau_div[cell_state] * (1 - sigma_div), tau_div[cell_state] * (1 + sigma_div))
				)
				@removeAgent()
			end
		end

		# Differentiation: B -> C at the constant rate q (no feedback)
		if d_on == true && cell_state == 2
			ran = CBMDistributions.uniform(0, 1)
			if ran < q * dt
				cell_state = 3
			end
		end

	end,


	agentAlg = CBMIntegrators.Heun(),
	neighborsAlg = CBMNeighbors.CellLinked(cellEdge = 10.0),
);

# ---- notebook cell 6ff203b2 ----------------------------------------------------
# ORIGINAL MINE

parameters_tfm = Dict(
	:range => 1.2,
	:mu => 2 * [1 1 1; 1 1 1 ; 1 1 1],
	:r => 1,
	:lambda => 1,
	:tau_div => 5 * [1, 1, 1],
	:sigma_div => 0.5,
	:olap => 0.75,
	# :p => 0.25,		# removed: no A cells in the B/C model
	:q => 0.125,
	# :k => 4.76,		# removed: no differentiation feedback
	:f_rep => 2.5 * [1 1 1; 1 1 1 ; 1 1 1],
	:f_att => 		[1 1 1; 1 1 1 ; 1 1 1],
	:t0 => 2,
	:r0 => 5,
	:f0 => 20
);

dt_tfm = 0.002
save_each_tfm = round(Int64, 0.25 / dt_tfm);

# ---- notebook cell 93111bf3 ----------------------------------------------------
# ORIOLA PAPER

parameters_do = Dict(
	:range => 1,		# changes	
	:mu => 
        [2.0 2.0 2.0
         2.0 2.0 2.0
         2.0 2.0 2.0],
	:r => 1,
	:lambda => 1,
	:tau_div => [300, 300, 750],	# changes
	:sigma_div => 0.5,	
	:olap => 0.3,			# changes
	# :p => 0.009333,		# removed: no A cells in the B/C model
	:q => 0.004667,			# changes
	# :k => 5,			# removed: no differentiation feedback
	# changes
	:f_rep =>
        [3 3 3 ; 
         3 3 3 ; 
         3 3 5.1],
	:f_att => 
        [6   6   2.4; 
         6   6   1.6; 
         2.4 1.6 4.8],
	:t0 => 2,
	:r0 => 5,
	:f0 => 20
);
dt_do = 0.015;
save_each_do = round(Int64, 1 / dt_do);
