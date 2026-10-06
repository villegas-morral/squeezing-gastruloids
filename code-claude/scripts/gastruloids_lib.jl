# =============================================================================
# gastruloids_lib.jl -- model, parameters and analysis for the squeezing project
#
# EXPORTED FROM squeezing-gastruloids-claude.ipynb. Do not edit by hand; edit the
# notebook and re-export with scripts/export_lib.ps1 so the two cannot drift.
#
# Usage:
#   using Random, Statistics, Printf, Distributions, ProgressMeter
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
	com.cell_state = 1
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
"""
	push_fA!(com)

Refresh the model parameter `:fA` with the live fraction of A cells.

`agentRule` reads the global A fraction as `ni_a = fA`. It cannot compute it
itself: a `com.cell_state` inside the rule is resolved as a global lookup in
`Main`, not as the community being stepped. Call this once per step, before
`agentStepRule!`. Model parameters can be written while the community is loaded
on the platform, so this works inside the stepping loops.
"""
@inline function push_fA!(com)
	n = com.N[1]
	cs = com.cell_state
	nA = 0
	@inbounds for i in 1:n
		cs[i] == 1 && (nA += 1)
	end
	com.parameters[:fA][1] = nA / n
	return nothing
end;

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
		push_fA!(com)
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
function initialize_diff!(com; 
    g_on = true,
    t_reset = true,
    center = 0.8,
    b0 = 0.75,
    a0 = 0.15,
    cluster = true
    )

    if cluster
        # The A,B core is selected as a CYLINDER centred on the aggregate, not
        # as a sphere centred on the origin. Two reasons:
        #  * a confined aggregate sits at z ~ L/2, not at z = 0, so a sphere
        #    x^2+y^2+z^2 < (center*r_agg)^2 sits BELOW the tissue. At layers=4,
        #    N0=100 it seeded only 14% of the cells, all of them squashed
        #    against the bottom plate;
        #  * the experiment squeezes a sphere that already has a central core,
        #    so the core should span the FULL thickness of the pancake, which
        #    is what a cylinder gives.
        # r_agg is measured from the cells rather than passed in, and the 98th
        # percentile is used so a few stray cells cannot inflate it.
        xc, yc = mean(com.x), mean(com.y)
        rho = sqrt.((com.x .- xc).^2 .+ (com.y .- yc).^2)
        r_agg = quantile(rho, 0.98)

        com.cell_state .= 3
        r_center_2 = (center * r_agg)^2
        central_indices = findall(i -> (com.x[i]-xc)^2 + (com.y[i]-yc)^2 < r_center_2, eachindex(com.x))
        if !isempty(central_indices)
            shuffle!(central_indices)
            n_center = length(central_indices)
            n_b = round(Int, b0 * n_center)
            n_a = round(Int, a0 * n_center)
            com.cell_state[central_indices[1:n_b]] .= 2
            com.cell_state[central_indices[n_b+1:minimum([n_b+1+n_a, n_center])]] .= 1
        end
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

end;

# ---- notebook cell 9154a7b7 ----------------------------------------------------
function differentiate!(com, save_each, tf;
    preallocate = 2 * com.N * 2 ^ (tf / minimum(com.tau_div)),
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
			    println("$(round(100 * i / steps))% \t N=$(com.N)")
                plot_pancake(com, color_map_alpha, length(com), length(com); savefig=false, n_views=1, n=1)
            end
			progress += step_control
		end
        push_fA!(com)
        agentStepDE!(com)
        agentStepRule!(com)
        update!(com)
        computeNeighbors!(com)
        if i % save_each == 0
            saveRAM!(com)
        end
        if all(com.N .> preallocate)
            println(message1)
            break
        end
    end        

	if steps % save_each != 0
		saveRAM!(com)
    end
    
	bringFromPlatform!(com)
    
    if show_progress
        IJulia.clear_output(true)
        println("Final N: $(com.N)")
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

"Save the trajectory of `com` as plain arrays (session-independent)."
function save_frames(file, com)
	d  = getParameter(com, [:x, :y, :z, :r, :cell_state])
	nf = length(com)
	isempty(dirname(file)) || mkpath(dirname(file))
	jldsave(file;
		x = d[:x], y = d[:y], z = d[:z], r = d[:r], cell_state = d[:cell_state],
		t = [com[k].t for k in 1:nf], N = [com[k].N for k in 1:nf])
	return file
end;

"Load what `save_frames` wrote. Returns a NamedTuple of per-frame vectors."
load_frames(file) = jldopen(file, "r") do fh
	(x = fh["x"], y = fh["y"], z = fh["z"], r = fh["r"],
	 cell_state = fh["cell_state"], t = fh["t"], N = fh["N"])
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
	fig = Figure(size = (320*length(ks), 340))
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

    cut_top = true    drop the cells above the mid-plane in the TOP view, so the
                      islands are visible instead of the C roof. It is a
                      cutaway, not a slice -- the side view is left intact, so
                      the constant slab thickness is still checkable by eye.
"""
function plot_pancake_frames(d, color_map, ks;
		n_views = 2, showtime = true, shownumbers = false, cut_top = true,
		boxsize = nothing, filename = "pancake_frames", savefig = false,
		labels = nothing, size_ratio = 1.0)

	views = [
		(perspective = pi/2,        show_z = false, show_x = true,  name = "top"),
		(perspective = 2*pi+pi/32,  show_z = true,  show_x = false, name = "side"),
		(perspective = -pi/2,       show_z = false, show_x = true,  name = "bottom"),
	][1:n_views]

	bs = boxsize === nothing ?
		maximum(maximum(d.x[k]) - minimum(d.x[k]) for k in ks) / 1.9 : boxsize
	labelsize = 26*size_ratio
	fig = Figure(size = (round(Int, 520*size_ratio*length(ks)),
						 round(Int, 470*size_ratio*n_views)), figure_padding = 20)

	for (j, view) in enumerate(views), (i, k) in enumerate(ks)
		x, y, z, r, s = d.x[k], d.y[k], d.z[k], d.r[k], d.cell_state[k]
		zc = mean(z)
		keep = (view.name == "top" && cut_top) ? findall(v -> v <= zc + 0.2, z) :
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
		meshscatter!(ax, x[o] .- mean(x), y[o] .- mean(y), z[o] .- zc;
			markersize = r[o], color = [color_map[s[m]] for m in o])
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
# --- two-aggregate fusion assay (Sotiris' experiment) --------------------------
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
	c.fA = 0.0
	c.abm.neighbors.cellEdge = (2*maximum(c.r))*maximum(c.mu)
	loadToPlatform!(c, preallocateAgents = round(Int, 1.5n))
	while c.N[1] < n
		push_fA!(c)
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
		pairs_ = [(1,1),(2,2),(3,3),(2,3),(1,2),(1,3)], labels = ["A","B","C"],
		verbose = true)
	agg = Dict(st => make_aggregate(n, st, mu_, fatt_, frep_) for st in 1:3)
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
		c.confinement = false; c.fA = 0.0
		c.alpha_ecm = 0.9 * c.lambda
		setfield!(c, :t, 0.0)
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

# Reference (Sotiris) vs the set used below, at 24 h:
#   AA 1.614 / 1.702    BB 0.537 / 0.531    CC 0.562 / 0.536
#   BC 0.947 / 0.744    AB 0.470 / 0.138    AC 0.657 / 0.510


# ---- notebook cell 8a8800de ----------------------------------------------------
# --- batch driver: one confined run end to end --------------------------------
#
# Everything the Initialization + Differentiation cells do, as a function, so a
# phase scan does not need the notebook's global state. Writes a *_frames.jld2.
#
# Note `show_progress = false` throughout: the interactive `differentiate!`
# renders a plot_pancake every 10%, which is fine for one run in a notebook and
# very much not fine for a scan of twenty.

"The unified parameter set (see the Differentiation cell for the derivation)."
const MU_SET   = 2 * [1 0.91 1; 0.91 1 1; 1 1 1.07];
const FATT_SET = 4 * [0.9 1.0 0.6; 1.0 1.6 0.4; 0.6 0.4 0.75];
const FREP_SET = 3 * [3 1 1; 1 1 1; 1 1 1.7];

"""
	wall_settings(layers)

`(rep_wall, dt)` for a confined run.

A monolayer needs the stiff wall or the sheet buckles to ~1.4 cells thick, and
the stiff wall needs a smaller timestep. The naive linear-stability bound
`rep_wall*dt < 2` is NOT safe: (150, 0.0075) satisfies it and still sends the
positions to NaN at layers = 1 from N0 = 300 up, because the spring switches on
and off discontinuously at the contact point rather than acting as a linear
restoring force. Measured safe: `rep_wall*dt <= 0.75`, verified at layers = 1 up
to N0 = 750 (z-spread 0.139, i.e. still one cell thick).
Use scripts/probe_stability.jl if you change rep_wall.
"""
wall_settings(layers) = layers <= 1 ? (150.0, 0.005) : (60.0, 0.015);

"""
	run_confined(layers, N0; tf = 625.0, seed = 2345, fscale = 1.0,
	             mu_ = MU_SET, fatt_ = FATT_SET, frep_ = FREP_SET,
	             outdir = "coms", save_every = 40.0)

Grow `N0` cells inside the plates, seed the A,B core, differentiate for `tf`,
and save the trajectory with `save_frames`. Returns `(file, timeline, verdict)`.

`fscale` multiplies `fatt_` — the island-size dial. At `layers = 1`, `N0 = 500`:
x0.70 loose porous patches, x1.00 compact stable islands, x1.40 one blob.
"""
function run_confined(layers, N0; tf = 625.0, seed = 2345, fscale = 1.0,
		mu_ = MU_SET, fatt_ = FATT_SET, frep_ = FREP_SET,
		outdir = "coms", save_every = 40.0, verbose = true)

	rw, dtd = wall_settings(layers)
	Lgap = (1.05*2.0)*layers
	tag  = @sprintf("L%.1f_N%d_f%.2f_s%d", layers, N0, fscale, seed)

	Random.seed!(seed)
	c = Community(model, N = 1, dt = dt_tfm, simBox = [-200.0 200; -200 200; -200 200])
	for (par, val) in pairs(parameters_tfm); c[par] = val; end
	c.cell_state = 1
	c.x = 0.0; c.y = 0.0; c.z = 0.0
	c.vx = 0.0; c.vy = 0.0; c.vz = 0.0
	c.t_div = 1
	c.g_on = true; c.d_on = false
	c.confinement = true; c.L = Lgap; c.rep_wall = rw
	c.alpha_ecm = 0.9*c.lambda[1]
	c.fA = 0.0
	grow_size!(c, save_each_tfm, N0; show_progress = false)

	for (par, val) in pairs(parameters_do); c[par] = val; end
	setfield!(c, :dt, dtd)
	c.mu = mu_; c.f_att = fscale .* fatt_; c.f_rep = frep_
	c.confinement = true; c.L = Lgap; c.rep_wall = rw
	c.alpha_ecm = 0.9*c.lambda[1]

	se = round(Int, save_every/dtd)
	stabilize!(c, se, 100; show_progress = false)
	initialize_diff!(c; g_on = true, cluster = true)
	stabilize!(c, se, 25; show_progress = false)
	m1 = length(c)
	differentiate!(c, se, tf; show_progress = false)

	file = save_frames(joinpath(outdir, "com_$(tag).jld2"), c)

	d  = load_frames(file)
	tl = island_timeline(d, layers; from = m1, verbose = false)
	v  = island_verdict(tl)
	if verbose
		e = tl[end]
		@printf("%-28s -> %-9s  n_isl=%d fmax=%.2f clust=%.2f fAB=%.3f  ",
				tag, v, e.n_isl, e.frac_max, e.clustered, e.fAB)
		confinement_report(d.z[end], Lgap)
	end
	(file = file, timeline = tl, verdict = v)
end;

"""
	phase_scan(layers_v, N0_v; kw...)

`run_confined` over a grid, tabulated. This is the layers-vs-N0 phase diagram:
below a critical N0 the A,B mass stays one central blob, above it breaks into
islands, and the critical N0 grows with the number of layers.

Measured (seed 2345, tf 625, fscale 1.0): islands from N0 = 500 / 1700 / 2500 /
4000 at layers = 1 / 2 / 3 / 4; one blob at 750 / 1700 / 2500 respectively.
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

# ---- notebook cell b3813774 ----------------------------------------------------
# --- why the threshold moves with the number of layers -------------------------
#
# Break-up is NOT set by cell number, it is set by the ASPECT RATIO of the A,B
# domain. The domain is a flat disc of thickness h (the slab) and in-plane
# diameter D, and like a liquid sheet it only breaks up once it is wide enough
# relative to its thickness.
#
# Measured with the unified set (f_att x1.00), taking the largest D/h each run
# reached while the domain was still one piece:
#
#     layers 4, N0 2500   D/h 3.16   one blob
#     layers 3, N0 1700   D/h 4.08   one blob
#     layers 2, N0  750   D/h 4.75   one blob
#     layers 4, N0 4000   D/h 5.04   ISLANDS
#     layers 3, N0 2500   D/h 5.83   ISLANDS
#     layers 2, N0 1700   D/h 9.13   ISLANDS
#     layers 1, N0  500   D/h 16.9   ISLANDS
#
# so the threshold is D/h ~ 5. (At f_att x0.45 it was ~3.4: stronger adhesion
# holds the domain together to a wider aspect ratio. Recalibrate with
# `aspect_report` if you change the adhesion scale.)
#
# In cell counts. The domain holds n_AB cells in h/2r layers, so
#     D = 4*sqrt( n_AB * 2r / (h*pi) )   =>   D/h = 5  <=>  n_AB ~ 22 * layers^3
# which predicts 22 / 173 / 583 / 1382 cells for layers = 1..4. That separates
# every run above: the blobs peaked at 181, 422 and 596 cells, the islands ran
# at 906, 945 and 1474.
#
# THE h^3 IS THE EXPERIMENTAL RULE. "If we increase the number of layers, then we
# need to increase the number of cells in order to get patterns" -- and the
# coefficient is now predictable rather than fitted. The corresponding N0 at
# squeezing is larger than n_AB and depends on how much the aggregate grows
# before the readout; measured brackets are in the phase table.

"""
	ab_domain_aspect(x, y, z, s; cut = 2.6)

`(D, h, ratio, n_big)` for the largest A,B domain in one frame: equivalent-disc
in-plane diameter `D`, slab thickness `h`, and `D/h`.
"""
function ab_domain_aspect(x, y, z, s; cut = 2.6)
	sel = findall(i -> s[i] == 1 || s[i] == 2, eachindex(s))
	isempty(sel) && return (D = NaN, h = NaN, ratio = NaN, n_big = 0)
	h  = quantile(z, 0.99) - quantile(z, 0.01) + 2.0
	cs = cluster_sizes(x[sel], y[sel], z[sel], cut)
	nbig = cs[1]
	per_layer = nbig / max(h/2.1, 1.0)
	D = 2 * sqrt(per_layer/pi) * 2.0
	(D = D, h = h, ratio = D/h, n_big = nbig)
end;

"Predicted A,B domain size at break-up, from the D/h threshold."
nAB_threshold(layers; ratio = 5.0, r = 1.0) =
	(ratio/4)^2 * pi * (2.1*layers*r)^3 / (2r);

"""
	aspect_report(d, layers)

`D/h` frame by frame alongside the island classification, so you can read off
the aspect ratio at which THIS parameter set breaks up.
"""
function aspect_report(d, layers; from = 1)
	@printf("%7s %7s %7s %7s %7s  %s\n", "t", "n_big", "D", "h", "D/h", "class")
	for k in from:length(d.x)
		x, y, z, s = d.x[k], d.y[k], d.z[k], d.cell_state[k]
		a = ab_domain_aspect(x, y, z, s)
		cls = classify_islands(island_stats(x, y, z, s, layers))
		@printf("%7.0f %7d %7.1f %7.2f %7.2f  %s\n", d.t[k], a.n_big, a.D, a.h, a.ratio, cls)
	end
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
		# Differentiation
		:d_on => Bool,
		:b => Float64,
		:p => Float64,
		:q => Float64,
		:k => Float64,
		:fA => Float64,			# live fraction of A cells, pushed in by the caller
		# Reference values
		:t0 => Float64,
		:r0 => Float64,
		:f0 => Float64,
		# Confinement
		:confinement => Bool,
		:L => Float64,
		:rep_wall => Float64,
	),


	# Agent parameters
	agent = Dict(
		:t_div => Float64,
		:ni => Int64,
		:cell_state => Int64,
		:r => Float64,
		# Mechanics
		:vx => Float64,
		:vy => Float64,
		:vz => Float64,
		:fx => Float64,
		:fy => Float64,
		:fz => Float64,
		# Differentiation
		:ni_a => Float64,
		:r_ab => Float64,
		:r_bc => Float64,
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

		# Confinement: HOOKEAN plates.
		# The previous version applied a CONSTANT force (fz += rep_wall) while
		# the cell was inside the contact zone. A constant force has no stable
		# fixed point there: it drives the cell until z >= r and then switches
		# off. At layers = 1 the gap between "off at the bottom" and "on at the
		# top" is only L - 2r = 0.1, while one restoring step moves the cell
		# rep_wall*dt/(lambda*ni) ~ 0.075, so cells chattered across the gap and
		# 93% of them ended up OUTSIDE the plates -- the "monolayer" was really
		# ~1.95 cells thick, and that buckled sheet is why layers = 1 never
		# produced islands. A spring vanishes at contact, so the contact point
		# IS the fixed point, and it cannot be outrun by in-plane crowding.
		# Stability: the restoring step is rep_wall*(pen/r)*dt/(lambda*ni), and
		# for ni < 2 the code below drops the 1/ni. The linear bound that suggests
		# is rep_wall*dt < 2, and it is NOT safe -- the spring switches on and off
		# discontinuously at the contact point, so the linear analysis does not
		# apply. (150, 0.0075) gives rep_wall*dt = 1.125 and still sends positions
		# to NaN at layers = 1 from N0 = 300 up, surfacing as
		#     InexactError: Int64(NaN)   (from the neighbour grid)
		# MEASURED SAFE: rep_wall * dt <= 0.75
		#     layers  = 1  -> rep_wall 150, dt 0.005   (checked up to N0 = 750)
		#     layers >= 2  -> rep_wall  60, dt 0.015
		# See wall_settings() and scripts/probe_stability.jl.
		if confinement
			if z < r
				fz += rep_wall * (r - z) / r
			end

			diz_top = L - z
			if diz_top < r
				fz -= rep_wall * (r - diz_top) / r
			end
		end

		# Equations of motion
		ni_relv_min = 2
		if ni < ni_relv_min
			inv = 1 / lambda
			vx = inv * fx
			vy = inv * fy
			vz = inv * fz
		else
			inv = 1 / (lambda * ni)
			vx = inv * (fx + vsum_x*alpha_ecm)
			vy = inv * (fy + vsum_y*alpha_ecm)
			vz = inv * (fz + vsum_z*alpha_ecm)
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

		# Differentiation
		if d_on == true && cell_state != 3
			# ni_a is the GLOBAL fraction of A cells. It used to be written as
			#     ni_a = count(==(1), com.cell_state) / com.N
			# but `com` inside this quote is not the community passed to
			# agentStepRule! -- the generated kernel resolves it as a global
			# lookup in Main. Evolving a community under any other name (or a
			# deepcopy taken while `com` still pointed at the grown, all-type-1
			# aggregate) made ni_a come out 1.0, so 1 + k*ni_a = 6 and BOTH
			# differentiation rates ran 6x too slow. It was also O(N) work per
			# agent per step, i.e. O(N^2) per step.
			# Now it is a model parameter that the stepping loops refresh once
			# per step (see push_fA! in grow_size!/stabilize!/differentiate!).
			ni_a = fA
			ran = CBMDistributions.uniform(0, 1)

			if cell_state == 1
				r_ab = p / (1 + k * ni_a)
				if ran < r_ab * dt
					cell_state = 2
				end

			elseif cell_state == 2
				r_bc = q / (1 + k * ni_a)
				if ran < r_bc * dt
					cell_state = 3
				end
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
	:p => 0.25,
	:q => 0.125,
	:k => 4.76,
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
	:p => 0.009333,			# changes
	:q => 0.004667,			# changes
	:k => 5,			# changes
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

# ---- notebook cell f3b763eb ----------------------------------------------------
# UNCONFINED POLARISATION -- the sanity check that the island parameters are
# reasonable. Set `confined = false` in the Initialization cell and re-run from
# there, or call the self-contained version below.
#
# With the unified set (N0 = 300, tf = 700) the A,B centroid ends ~0.85 aggregate
# radii off centre with r_AB > r_C, i.e. A,B segregated to one pole and sitting
# at the surface, C on the other side. `offset` is still climbing at the end of
# the run, so a longer tf gives a more completely separated state.

function run_unconfined(N0 = 300; tf = 700.0, seed = 2345, fscale = 1.0,
		mu_ = MU_SET, fatt_ = FATT_SET, frep_ = FREP_SET,
		outdir = "coms", save_every = 40.0, show_progress = false)

	Random.seed!(seed)
	c = Community(model, N = 1, dt = dt_tfm, simBox = [-150.0 150; -150 150; -150 150])
	for (par, val) in pairs(parameters_tfm); c[par] = val; end
	c.cell_state = 1
	c.x = 0.0; c.y = 0.0; c.z = 0.0
	c.vx = 0.0; c.vy = 0.0; c.vz = 0.0
	c.t_div = 1
	c.g_on = true; c.d_on = false
	c.confinement = false; c.L = 4.2; c.rep_wall = 60
	c.alpha_ecm = 0.9*c.lambda[1]
	c.fA = 0.0
	grow_size!(c, save_each_tfm, N0; show_progress = show_progress)

	for (par, val) in pairs(parameters_do); c[par] = val; end
	setfield!(c, :dt, 0.015)
	c.mu = mu_; c.f_att = fscale .* fatt_; c.f_rep = frep_
	c.confinement = false
	c.alpha_ecm = 0.9*c.lambda[1]

	se = round(Int, save_every/0.015)
	stabilize!(c, se, 100; show_progress = show_progress)
	# NOTE the core is seeded as a SPHERE here, not the cylinder used when
	# confined: unconfined the aggregate is a ball, so a ball-shaped core is
	# what squeezing would have produced.
	xc, yc, zc = mean(c.x), mean(c.y), mean(c.z)
	rho = sqrt.((c.x .- xc).^2 .+ (c.y .- yc).^2 .+ (c.z .- zc).^2)
	r_agg = quantile(rho, 0.98)
	c.cell_state .= 3
	idx = shuffle!(findall(i -> rho[i] < 0.8*r_agg, eachindex(rho)))
	nb = round(Int, 0.75*length(idx)); na = round(Int, 0.15*length(idx))
	c.cell_state[idx[1:nb]] .= 2
	c.cell_state[idx[nb+1:min(nb+1+na, length(idx))]] .= 1
	c.g_on = true; c.d_on = true
	setfield!(c, :t, 0.0)
	for i in 1:c.N[1]
		sr = 2*c.tau_div[c.cell_state[i]]*c.sigma_div[] + c.dt[1]
		c.t_div[i] = CBMDistributions.uniform(c.dt[1], sr)
	end
	saveRAM!(c)
	m1 = length(c)
	# show_progress = false by default: differentiate!'s progress display calls
	# IJulia.clear_output and plot_pancake, neither of which exists in a plain
	# `julia script.jl` session.
	differentiate!(c, se, tf; show_progress = show_progress)

	file = save_frames(joinpath(outdir, "polar_N$(N0)_f$(fscale).jld2"), c)
	d = load_frames(file)
	@printf("%7s %6s %7s %7s %7s\n", "t", "fAB", "offset", "r_AB", "r_C")
	for k in m1:length(d.x)
		st = polar_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k])
		@printf("%7.0f %6.3f %7.3f %7.3f %7.3f\n", d.t[k], st.fAB, st.offset, st.r_AB, st.r_C)
	end
	(file = file, frames = d, m1 = m1)
end;

# res = run_unconfined(300);
# plot_frames(res.frames, color_map, round.(Int, range(res.m1, length(res.frames.x), length=6)))
