# =============================================================================
# claude_changes.jl -- every change Claude made to your notebook, as code.
#
# squeezing-gastruloids.ipynb          YOUR notebook, exactly as in commit f7fd270
# claude_changes.jl                    THIS FILE: the changes, to apply yourself
# squeezing-gastruloids_claude.ipynb   my full notebook from July (reference only)
#
# Read CATCH-UP.md (repo root) first: it explains each change and the evidence.
#
# HOW THIS FILE IS ORGANISED
#
#   PART 1  fixes to YOUR cells        a complete replacement for each cell, with
#                                      the lines that differ from yours tagged
#                                      `# [changed]` or `# [new]`
#   PART 2  the parameter set          constants + wall settings
#   PART 3  new tools                  new cells: saving, analysis, batch runs
#   PART 4  your run cells             the Initialization / Differentiation cells
#                                      as they would look -- inside a #= =#
#                                      block, so including this file does not
#                                      start a simulation
#
# Each part is ordered by importance. If you only take one thing, take 1.1.
#
# VERIFIED 2026-10-04: your notebook's package, function, model and parameter
# cells, followed by `include("claude_changes.jl")`, run the whole pipeline --
# confined runs at layers 1 and 2, unconfined polarisation, the fusion assay,
# and every plotting function.
#
# Extra packages needed beyond your package cell:   using Statistics
# (JLD2 you already load.)
# =============================================================================

using Statistics


# #############################################################################
# PART 1 -- FIXES TO YOUR CELLS
# #############################################################################

# -----------------------------------------------------------------------------
# 1.1  MODEL (cell 99c84e68) -- three changes
#
#   (a) `ni_a` no longer reads the global `com`                    [correctness]
#       Your line   ni_a = count(==(1), com.cell_state) / com.N
#       resolves `com` as whatever variable called `com` exists in Main, NOT the
#       community being stepped. Measured on 2026-10-04 with your exact model:
#         stepping `com` itself (your notebook flow)    ni_a = 0.14  (true 0.15) OK
#         stepping another community, global `com`
#         = an all-A aggregate                          ni_a = 0.89  -> rates ~4x too slow
#       So your notebook was getting the right answer, but any function that
#       steps a community under another name (every batch run, the fusion
#       assay) silently got the wrong rates. It was also O(N) per agent per
#       step, i.e. O(N^2) per step.
#       Fix: a model parameter `fA`, refreshed once per step by push_fA! (1.2).
#
#   (b) Hookean plates instead of a constant wall force           [correctness]
#       `fz += rep_wall` while z < r has no stable fixed point: it pushes until
#       z >= r, then switches off. At layers = 1 the band for cell centres is
#       L - 2r = 0.1 while one restoring step moves ~0.075, so cells chatter
#       across it: 93% of cells ended up OUTSIDE the plates and the "monolayer"
#       was a buckled sheet ~1.95 cells thick. (Layers >= 2 leaked 0.2-2%.)
#       A spring vanishes at contact, so contact IS the fixed point.
#       COST: a spring is softer than a constant force at small penetration --
#       see wall_settings (2.2) and CATCH-UP.md, "the wall at layers >= 2".
#
#   (c) nothing else. Agent parameters, forces, division, rates: all yours.
# -----------------------------------------------------------------------------
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
		:fA => Float64,			# [new] live fraction of A cells, set by push_fA!
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

		# [changed] Hookean plates: force proportional to the penetration.
		# Stability is NOT given by the linear bound rep_wall*dt < 2 (the
		# spring switches on/off discontinuously). See wall_settings.
		if confinement
			if z < r
				fz += rep_wall * (r - z) / r			# [changed] was: fz += rep_wall
			end

			diz_top = L - z
			if diz_top < r
				fz -= rep_wall * (r - diz_top) / r		# [changed] was: fz -= rep_wall
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
			ni_a = fA		# [changed] was: count(==(1), com.cell_state) / com.N
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


# -----------------------------------------------------------------------------
# 1.2  push_fA! [new] and grow_size! (cell f2f50d7d)
#
#   push_fA! goes with 1.1(a). Call it once per step, before agentStepRule!.
#   Model parameters can be written while the community is on the platform,
#   so it works inside the stepping loops.
#
#   grow_size! changes:
#     * push_fA!(com) before agentStepRule!                       [needed by 1.1a]
#     * the neighbour-grid size is set here                       [correctness]
#       You set cellEdge = 2r*max(mu) once, in initialize_growth, while mu is
#       still parameters_tfm's 2.0 -> cellEdge 4.0. The final mu goes up to
#       2.14, i.e. interactions out to 4.28 > 4.0, which the grid can miss.
#       Setting it at the top of every stepping function keeps it in sync with
#       whatever mu is current. (The line in initialize_growth is then
#       redundant but harmless -- keep it or drop it.)
#     * show_progress keyword: false for batch runs                [convenience]
# -----------------------------------------------------------------------------
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
	n_control=100, show_progress = true)							# [changed] show_progress

	com.abm.neighbors.cellEdge = (2 * maximum(com.r)) * maximum(com.mu)	# [new]

	loadToPlatform!(com, preallocateAgents = round(Int, n_cells * 1.1))
	i = 0
	n_total = 0

	while (com.N .< n_cells)
		if show_progress && com.N > n_total						# [changed]
			println("N > $(n_total)")
			n_total += n_control
		end
		i += 1
		push_fA!(com)												# [new]
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


# -----------------------------------------------------------------------------
# 1.3  initialize_diff! (cell f818a62d)                         [correctness]
#
#   The A,B core is selected as a CYLINDER about the aggregate's centroid,
#   not a sphere about the origin:
#     * a confined aggregate sits at z ~ L/2, not z = 0, so the sphere
#       x^2+y^2+z^2 < (center*r_agg)^2 sits partly BELOW the tissue. At
#       layers = 4, N0 = 100 it seeded only 14% of the cells, all squashed
#       against the bottom plate;
#     * the experiment squeezes a sphere that already has a core, so the core
#       should span the full thickness of the pancake -- a cylinder does that.
#   r_agg is measured from the cells (98th percentile of the in-plane radius,
#   so a few strays cannot inflate it) instead of being a keyword computed
#   from the x-extent.
#
#   Defaults changed to the values your Differentiation cell passes anyway
#   (center 0.45 -> 0.8, b0 0.7 -> 0.75, a0 0.2 -> 0.15).
#
#   NOTE for the UNCONFINED case a cylinder is the wrong shape -- run_unconfined
#   (3.7) seeds a sphere about the centroid instead.
# -----------------------------------------------------------------------------
function initialize_diff!(com;
    g_on = true,
    t_reset = true,
    center = 0.8,												# [changed] 0.45
    b0 = 0.75,													# [changed] 0.7
    a0 = 0.15,													# [changed] 0.2
    cluster = true
    )

	com.abm.neighbors.cellEdge = (2 * maximum(com.r)) * maximum(com.mu)

    if cluster
        xc, yc = mean(com.x), mean(com.y)							# [new]
        rho = sqrt.((com.x .- xc).^2 .+ (com.y .- yc).^2)			# [new]
        r_agg = quantile(rho, 0.98)									# [new]

        com.cell_state .= 3
        r_center_2 = (center * r_agg)^2
        central_indices = findall(i -> (com.x[i]-xc)^2 + (com.y[i]-yc)^2 < r_center_2, eachindex(com.x))	# [changed]
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


# -----------------------------------------------------------------------------
# 1.4  differentiate! (cell 9154a7b7)
#
#     * push_fA!(com) before agentStepRule!                       [needed by 1.1a]
#     * cellEdge set at the top (see 1.2)                         [correctness]
#     * show_progress keyword                                      [needed for scripts]
#       Your progress display calls IJulia.clear_output and plot_pancake, which
#       only exist inside Jupyter: in `julia script.jl` it crashes. Batch
#       functions pass show_progress = false.
#     * progress plot uses color_map_alpha (C transparent) so the  [cosmetic]
#       A,B pattern is visible while it runs.
# -----------------------------------------------------------------------------
function differentiate!(com, save_each, tf;
    preallocate = 2 * com.N * 2 ^ (tf / minimum(com.tau_div)),
    show_progress = true										# [new]
    )

	com.abm.neighbors.cellEdge = (2 * maximum(com.r)) * maximum(com.mu)	# [new]

    message1 = "Possible numerical instabilities.
            You might have to send argument preallocate = highnumber
            Default is preallocate = 2 * com.N * 2 ^ (tf / minimum(com.tau_div))"
    message2 = "Numerical instabilities.
            You might have to decrease the timestep"

	steps = round(Int64, tf / com.dt)
    step_control = steps / 10
    progress = step_control

    if com.g_on[]
        show_progress && println("Initial N: $(com.N)")			# [changed]
        loadToPlatform!(com, preallocateAgents = round(Int, preallocate))
    else
        loadToPlatform!(com)
    end

    for i in 1:steps
		if i > progress
            if show_progress										# [new]
                IJulia.clear_output(true)
			    println("$(round(100 * i / steps))% \t N=$(com.N)")
                plot_pancake(com, color_map_alpha, length(com), length(com); savefig=false, n_views=1, n=1)	# [changed] color_map_alpha
            end
			progress += step_control
		end
        push_fA!(com)												# [new]
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

    if show_progress												# [new]
        IJulia.clear_output(true)
        println("Final N: $(com.N)")
    end

end;


# -----------------------------------------------------------------------------
# 1.5  stabilize! (cell f06d6c5e)
#     * cellEdge set at the top (see 1.2)                         [correctness]
#     * show_progress keyword (the progress bar is noise in batch runs)
#     * the unused step_control / progress variables dropped
#   No push_fA! needed: stabilize! never calls agentStepRule!.
# -----------------------------------------------------------------------------
function stabilize!(com, save_each, tf; show_progress = true)		# [changed]

	com.abm.neighbors.cellEdge = (2 * maximum(com.r)) * maximum(com.mu)	# [new]

	steps = round(Int64, tf / com.dt)
	i = 0

	prog_bar = show_progress ? Progress(steps, desc="Stabilizing") : nothing	# [changed]

	loadToPlatform!(com)

	while (i<steps)
		i += 1
		agentStepDE!(com)
		update!(com)
		computeNeighbors!(com)
		if i % save_each == 0
			saveRAM!(com)
		end
		show_progress && next!(prog_bar)							# [changed]
	end
	if i % save_each != 0
		saveRAM!(com)
	end
	bringFromPlatform!(com)

end;


# #############################################################################
# PART 2 -- THE PARAMETER SET
# #############################################################################

# -----------------------------------------------------------------------------
# 2.1  Mechanics: one set for fusion, polarisation and patterns at layers 1-4
#
# This is Sotiris' FUSION set (the one commented out in your Differentiation
# cell) with TWO diagonal f_rep entries softened, and nothing else:
#     mu      FUSION, untouched
#     f_att   FUSION, untouched
#     f_rep   FUSION except   [2,2]  6 -> 3     (B self-repulsion)
#                             [3,3]  6 -> 5.1   (C self-repulsion)
#     p, q    parameters_do, untouched (p = 2q)
#
# Why the softening: with FUSION's stiff f_rep the A,B phase never coarsens --
# at layers = 1, N0 = 500 it stays salt-and-pepper (largest cluster 12 cells)
# and differentiates away. Softer BB and CC let it round up into islands.
#
# Fusion assay at 24 h (COM separation / (Rg1+Rg2); lower = more fused):
#             AA     BB     CC     BC     AB     AC
#   FUSION   1.614  0.537  0.562  0.947  0.470  0.657
#   this set 1.702  0.531  0.536  0.744  0.138  0.510
# BC is the least-fused of the pairs that fuse, by a clear margin: the "glue"
# phenotype, but more rounded than the reference (0.744 vs 0.947).
# -----------------------------------------------------------------------------
const MU_SET   = 2 * [1    0.91 1   ;
                      0.91 1    1   ;
                      1    1    1.07];
const FATT_SET = 4 * [0.9  1.0  0.6 ;
                      1.0  1.6  0.4 ;
                      0.6  0.4  0.75];
const FREP_SET = 3 * [3 1 1;
                      1 1 1;
                      1 1 1.7];

# 1 hour = 15 model time units  (p = 0.009333 <-> 7.1 h, q = 0.004667 <-> 14.1 h)
const HOUR = 15.0;

# -----------------------------------------------------------------------------
# 2.2  wall_settings(layers) -> (rep_wall, dt) for a confined run
#
# A monolayer needs a stiff wall or the sheet buckles, and a stiff wall needs a
# small timestep. Measured, not derived:
#   (150, 0.005)  layers = 1: z-spread ~0.13, i.e. one cell thick, stable to N0 = 750
#   (150, 0.0075) blew up (positions -> NaN, "InexactError: Int64(NaN)") in
#                 the July scan from N0 = 300; on 2026-10-04 one N0 = 500 run
#                 survived but with z-spread 0.45 -- so it is unreliable, and
#                 worse at its job even when it survives
#   (60, 0.015)   layers >= 2
# Rule of thumb from those runs: rep_wall * dt <= 0.75.
# Use scripts/probe_stability.jl if you change rep_wall.
# -----------------------------------------------------------------------------
wall_settings(layers) = layers <= 1 ? (150.0, 0.005) : (60.0, 0.015);


# #############################################################################
# PART 3 -- NEW TOOLS (add as new cells; none of them changes the model)
# #############################################################################

# -----------------------------------------------------------------------------
# 3.1  Saving and reloading runs
#
# `jldsave("com.jld2"; com)` does NOT round-trip. ABM() eval's the step
# functions as ANONYMOUS closures in Main, and the pickled Community refers to
# them by gensym name. Reloading needs a session whose gensym counter produced
# the same names -- any kernel restart or different cell order breaks it:
#     Warning: type Main.#19#20 does not exist in workspace; reconstructing
#     MethodError: Cannot `convert` an object of type
#                  JLD2.ReconstructedSingleton{Symbol("#19#20")} to ... Function
# Plain arrays have no such dependency, load without the model defined, and
# are ~7x smaller. The cost: you cannot continue evolving a reloaded run.
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# 3.2  Plotting saved runs
#
# plot_frames          flat 2D, mid-plane slice. Use the slice: a top
#                      projection through the whole slab lets the C roof hide
#                      the islands underneath.
# plot_pancake_frames  YOUR plot_pancake look (Axis3 + meshscatter), for the
#                      array format. Kept separate because plot_pancake needs a
#                      live Community. `cut_top` drops cells above the
#                      mid-plane in the top view so the islands are visible.
#                      This is what made the figures for the experimentalists.
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# 3.3  Island detection (one frame)
#
# The A,B cells are split into single-linkage clusters (two cells linked if
# closer than `cut` = 2.6). Then:
#   n_isl      clusters of at least nmin_island(layers) cells
#   frac_max   largest cluster / all A,B cells      (1 => one blob)
#   clustered  fraction of A,B cells in real islands
# `clustered` is the discriminator, not frac_max: a salt-and-pepper dispersion
# also has a tiny frac_max. The size floor grows with layers so a monolayer
# island is not penalised for holding fewer cells at the same footprint.
# Classes were validated by eye against ~60 runs.
# -----------------------------------------------------------------------------
"Sizes of the single-linkage clusters of the points (cutoff distance), descending. O(N) via a grid."
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

"One frame -> (n_isl, frac_max, clustered, nAB, top) for the A,B cells."
function island_stats(x, y, z, s, layers; cut = 2.6)
	sel = findall(i -> s[i] == 1 || s[i] == 2, eachindex(s))
	nAB = length(sel)
	nAB == 0 && return (n_isl = 0, frac_max = NaN, clustered = 0.0, nAB = 0, top = Int[])
	cs  = cluster_sizes(x[sel], y[sel], z[sel], cut)
	big = filter(>=(nmin_island(layers)), cs)
	(n_isl = length(big), frac_max = cs[1]/nAB,
	 clustered = sum(big; init = 0)/nAB, nAB = nAB, top = cs[1:min(end,8)])
end;

"`gone` / `dispersed` (salt and pepper) / `islands` / `single` (one central cluster)."
function classify_islands(st)
	st.nAB == 0        && return "gone"
	st.clustered < 0.6 && return "dispersed"
	(st.n_isl >= 2 && st.frac_max <= 0.7) && return "islands"
	return "single"
end;

# -----------------------------------------------------------------------------
# 3.4  Whole-trajectory readouts
#
# The A,B phase is a TRANSIENT: A -> B -> C is one-way, so fAB decays to 0 and
# the islands live in a window (t ~ 200-650 with the nominal p, q). Scoring a
# single late frame is how the monolayer was first mis-read as "no islands".
# Score every frame (island_timeline) and ask for islands at the END
# (island_verdict); a run that fragments then re-merges is "transient".
# -----------------------------------------------------------------------------
"Score every frame of `d` (from load_frames). Returns a vector of NamedTuples."
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

"One word for a timeline: ISLANDS (at the end), transient (formed, then gone), or the final class."
function island_verdict(tl)
	isempty(tl) && return "empty"
	nis = count(r -> r.class == "islands", tl)
	tl[end].class == "islands" ? "ISLANDS" :
		(nis >= 3 ? "transient" : uppercase(tl[end].class))
end;

"""
	polar_stats(x, y, z, s) -> (offset, r_AB, r_C, fAB, R)

Unconfined polarisation, lengths in units of the aggregate radius R.
`offset`: A,B centroid distance from the centre (0 = central core, ->1 = one
pole). `r_AB > r_C` means A,B has moved out to the surface.
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
	confinement_report(z, L; r = 1.0)

Are the plates confining? Prints the z-spread of cell centres against the
nominal band [r, L-r], and the effective number of layers (spread + 2r)/2.1r.
Look at the SPREAD, not the "outside" percentage: with Hookean plates every
cell in contact sits slightly inside the wall, so at layers = 1 "outside" is
~80% for a perfectly flat sheet (spread 0.13).
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

# -----------------------------------------------------------------------------
# 3.5  Adhesion energetics -- check a parameter set BEFORE simulating it
#
# Integrating your pair force over the attractive branch (2r < d < mu*2r):
#     W_ij = f_att[i,j] * 2r * A(mu_ij),   A(mu) = (1+mu)*ln(mu) - 2*(mu-1)
# f_rep does not appear (it only acts below contact). Then
#     sigma_i  = W_ii / 2                  surface tension of type i
#     gamma_ij = (W_ii + W_jj)/2 - W_ij    interfacial tension
# Steinberg: gamma > s_i + s_j -> separate; |s_i - s_j| > gamma -> engulfment;
# otherwise partial engulfment, a stable doublet ("glue"). BC must be "glue".
# Final set: AB 0.405, AC 0.232, BC 0.581 -- all three "glue".
# -----------------------------------------------------------------------------
A_mu(mu) = (1 + mu)*log(mu) - 2*(mu - 1);
adhesion_work(f_att, mu; r = 1.0) = f_att .* (2r) .* A_mu.(mu);

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

# -----------------------------------------------------------------------------
# 3.6  Two-aggregate fusion assay (Sotiris' doublet experiment)
#
# Two pure aggregates of n cells, placed just touching, no division and no
# differentiation. Metric: COM separation / (Rg1 + Rg2) -- ~1.6 at contact,
# ~0.5 once rounded into one sphere. LOWER = MORE FUSED. Look at the ORDERING:
# BC should be the least fused of the pairs that fuse (AA never fuses).
#
# make_aggregate relaxes each aggregate under ITS OWN type's mechanics: growth
# uses parameters_tfm (olap 0.75, heavily overlapping) and dropping that
# straight into a stiff f_rep makes the aggregates explode on contact.
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# 3.7  Batch drivers: a whole run as one function call
#
# Your Initialization + Differentiation cells, as functions, so a scan does not
# depend on the notebook's global state. They use LOCAL communities (`c`),
# which is exactly the situation where the old `ni_a` line (1.1a) broke.
#
# run_confined(layers, N0)   grow N0 cells between the plates, seed the core,
#                            differentiate for tf, save, score.
# run_unconfined(N0)         the polarisation check. Seeds a SPHERICAL core
#                            (unconfined the aggregate is a ball).
# phase_scan(layers_v, N0_v) run_confined over a grid.
#
# Measured, with seed 2345 and tf 625 (20 runs, results/phase_map.csv):
#     layers   no islands up to   islands from
#        1            200              300
#        2            750             1150
#        3           1700             2100
#        4           2500             3200
# Polarisation (N0 300, tf 700): offset 0.945, r_AB 0.973 vs r_C 0.637.
# -----------------------------------------------------------------------------
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
	# spherical core about the centroid
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

# -----------------------------------------------------------------------------
# 3.8  Phase diagram
#
# plot_phase(rows)    two categories, as you asked in July:
#                     "patterns" (islands at the end of the run) and
#                     "central cluster" (one A,B mass; `transient` runs are
#                     drawn here, since that is what they end as).
#                     A DISPERSED run would be neither: drawn grey, with a warning.
# critical_N0(rows)   the bracket per layer count.
# -----------------------------------------------------------------------------
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

# NOT carried over from my July notebook: the `ab_domain_aspect` /
# `nAB_threshold` cell (b3813774). It encoded a claim that break-up happens at a
# fixed aspect ratio D/h ~ 5 of the A,B domain, so n_AB ~ 22 layers^3. The full
# 20-run map refuted it (no single D/h separates blobs from islands). It is
# still in squeezing-gastruloids_claude.ipynb if you want the measurement.


# #############################################################################
# PART 4 -- YOUR RUN CELLS, AS THEY WOULD LOOK
#
# Not executed by include(). Differences from yours, in order:
#   * `confined` switch; layers = 1 by default
#   * wall and timestep from wall_settings(layers)    (yours: rep_wall 30, dt_do)
#   * com.fA = 0.0 initialised
#   * a confinement sanity check after growth
#   * the unified parameter set instead of the "CONFINED PATTERNS" matrices
#   * initialize_diff!(...; cluster = true)            (yours: cluster = false,
#     which seeds no A,B core at all -- every cell stays A)
#   * tf 625 instead of 500, and a save_frames at the end
# #############################################################################
#=

# ---- Initialization (cell 345b549e) -----------------------------------------
model = model_new;

dt = dt_tfm;
save_each = save_each_tfm;
parameters = parameters_tfm

confined = true      # true: squeezed pancake -> islands;  false: polarisation
layers = 1

# N0 = cells AT SQUEEZING. Thresholds (seed 2345): islands from ~300 / 1150 /
# 2100 / 3200 cells at layers 1 / 2 / 3 / 4.
n_cells = confined ? (layers == 1 ? 500 : 1700) : 300;

rep_wall_use, dt_do_use = confined ? wall_settings(layers) : (60.0, 0.015)

Random.seed!(2345)
com = initialize_growth(parameters; dt = dt);

com.confinement = confined
com.L = (1.05 * (2*com.r)) * layers
com.rep_wall = rep_wall_use
com.fA = 0.0            # refreshed every step by push_fA!

com.alpha_ecm = 0.9 * com.lambda

grow_size!(com, save_each, n_cells)
m0 = length(com);

growncom = deepcopy(com);

confined && confinement_report(com.z, com.L[1])
plot_pancake(com, color_map, 1, m0; savefig=false, filename="growth", n_views=2)


# ---- Differentiation (cell d69b054b) ----------------------------------------
com = deepcopy(growncom)
Random.seed!(2345);

for (par, val) in pairs(parameters_do)
		com[par] = val
end

setfield!(com,:dt, dt_do_use)
save_each = round(Int64, 1 / com.dt[1])

com.mu    = MU_SET
com.f_att = FATT_SET
com.f_rep = FREP_SET

com.confinement = confined
com.L = (1.05 * (2*com.r)) * layers
com.rep_wall = rep_wall_use
com.alpha_ecm = 0.9 * com.lambda

stabilize!(com, save_each, 100)
initialize_diff!(com; g_on=true, cluster=true, center=0.8, a0=0.15, b0=0.75)
stabilize!(com, save_each, 25)
m1 = length(com);
plot_pancake(com, color_map, m0, m1; savefig=false, filename="growth", n_views=1)

growncom = deepcopy(com)


# ---- Run (cell 9f98a77f) ----------------------------------------------------
m1 = length(com);
differentiate!(com, save_each, 625)
m2 = length(com);
evolvedcom = deepcopy(com);

save_frames("coms/com_L$(layers)_N$(n_cells).jld2", com)

plot_pancake(com, color_map, m1, m2; savefig=false)
plot_pancake(com, color_map_alpha, m1, m2; savefig=false, n_views=1)

# plot_pancake's top view looks through the whole slab, so at 2+ layers the C
# roof hides the islands. The mid-plane slice / cutaway shows them:
fr = load_frames("coms/com_L$(layers)_N$(n_cells).jld2")
plot_frames(fr, color_map, round.(Int, range(m1, length(fr.x), length = 6)))
plot_pancake_frames(fr, color_map, round.(Int, range(m1, length(fr.x), length = 4)))

=#
