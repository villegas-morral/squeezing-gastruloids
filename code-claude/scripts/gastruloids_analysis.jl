# =============================================================================
# gastruloids_analysis.jl -- scoring and plotting only, NO model.
#
# EXPORTED FROM squeezing-gastruloids-claude.ipynb by scripts/export_lib.ps1.
# Do not edit by hand.
#
# Loads without CellBasedModels, so it is fast and works on any machine that can
# read the saved trajectories:
#   using Statistics, Printf, JLD2, CairoMakie
#   include("gastruloids_analysis.jl")
#
# (`save_frames` is the one function here that needs a live Community; the rest
# operate on the array format from `load_frames`.)
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
