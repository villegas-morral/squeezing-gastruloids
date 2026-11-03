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
	fig = Figure(size = (round(Int, 520*size_ratio*length(ks)),
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
	fig = Figure(size = (round(Int, 520*size_ratio*ncol),
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
2026, your choice: 0.2, then "a little less transparent"): the whole ball is drawn, B opaque and C translucent, so the B core
shows through the C cells around it. As in your original plot_pancake: one
meshscatter with RGBA colours and no `transparency = true` (whose order-independent
blending turns many translucent layers dark brown), B drawn first so the C spheres
in front of it blend over it. Needs GLMakie; CairoMakie draws each sphere as many
small triangles whose edges add up. Confined aggregates are drawn opaque with the
mid-plane cutaway (`ghost = nothing`, the default of draw_cells!).
"""
GHOST = 0.3;

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

"""
	draw_cells!(pos, c; view = :top, ghost = nothing, cut_top = true, label = "", labelsize = 18, bs = nothing)

One 3D panel in the plot_pancake look (Axis3 + meshscatter) at the figure position
`pos`, for the cells `c` (anything with x, y, z, r, cell_state):

	:top      seen from above
	:side     seen from the side (shows the slab thickness)
	:front    the side view turned by 90 degrees (x runs left to right)
	:oblique  a three-quarter view from above

With ghost = nothing (the default) the cells are opaque, and with cut_top the top
view leaves out the cells above the mid-plane, so the B layout shows on the cut face
(confined aggregates). With ghost = GHOST every view shows the whole aggregate with C
translucent (unconfined aggregates, where B polarises; November 2026).

`bs` is the half-width of the box (default: fits the aggregate); give every panel
of a figure the same `bs` to draw them at the same scale.
"""
function draw_cells!(pos, c; view = :top, ghost = nothing, cut_top = true, label = "", labelsize = 18, labelcolor = :black, bs = nothing)
	x, y, z, r, s = c.x, c.y, c.z, c.r, c.cell_state
	xc, yc, zc = mean(x), mean(y), mean(z)
	bs = bs === nothing ? box_halfwidth(c) : bs
	gh = ghost
	keep = (view == :top && cut_top && gh === nothing) ? findall(w -> w <= zc + 0.2, z) : collect(eachindex(z))
	ax = Axis3(pos, aspect = :data,
		azimuth = view == :oblique ? pi/4 : view == :front ? pi/2 : 0,
		elevation = view == :top ? pi/2 : view in (:side, :front) ? 2*pi + pi/32 : pi/5,
		xlabel = "", ylabel = "", zlabel = "",
		xgridvisible = false, ygridvisible = false, zgridvisible = false,
		xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
		xticksvisible = false, yticksvisible = false, zticksvisible = false,
		xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false,
		titlevisible = !isempty(label), title = label, titlesize = labelsize, titlegap = 4, titlecolor = labelcolor)
	cells_mesh!(ax, x[keep] .- xc, y[keep] .- yc, z[keep] .- zc, r[keep], s[keep]; ghost = gh)
	xlims!(ax, -bs, bs); ylims!(ax, -bs, bs); zlims!(ax, -bs, bs)
	ax
end;

"""
	plot_panels(dir; size_ratio = 0.6, labelsize = 18)

A final 3D figure from its data folder (title.txt, panels.tsv, <panel>.csv). Each
panel sits at (row, col) of panels.tsv and gets one axis row per view
("top", "side", "front", "oblique", or several joined by "+"). An optional 7th
column `bs` fixes the box half-width of that panel (the same value in every panel
= the same scale); an empty title or label draws none. Unconfined aggregates (views
top+front) are drawn whole with C translucent (GHOST), confined ones opaque with the
mid-plane cutaway.
"""
function plot_panels(dir; size_ratio = 0.6, labelsize = 18)
	lines_ = filter(l -> !isempty(strip(l)), readlines(joinpath(dir, "panels.tsv")))
	pans = [split(l, '\t') for l in lines_[2:end]]       # panel, row, col, views, label, file[, bs]
	pbs(p) = length(p) >= 7 && !isempty(strip(p[7])) ? parse(Float64, p[7]) : nothing
	views(p) = Symbol.(split(p[4], "+"))
	nrow = maximum(parse(Int, p[2]) for p in pans)
	ncol = maximum(parse(Int, p[3]) for p in pans)
	# every figure row is as tall as its panel with the most views
	nv = [maximum([length(views(p)) for p in pans if parse(Int, p[2]) == i]; init = 1) for i in 1:nrow]
	start = cumsum(vcat(1, nv[1:end-1]))
	title = isfile(joinpath(dir, "title.txt")) ? strip(read(joinpath(dir, "title.txt"), String)) : ""
	toff = isempty(title) ? 0 : 1
	# a side view is flat: an axis row with only side views gets 40 % of the height
	h = round(Int, 470*size_ratio)
	flat = [all(p -> parse(Int, p[2]) != i || length(views(p)) < jv || views(p)[jv] == :side, pans)
			for i in 1:nrow for jv in 1:nv[i]]
	fig = Figure(size = (round(Int, 520*size_ratio*ncol), sum(f ? round(Int, 0.4h) : h for f in flat) + 70*toff),
				 figure_padding = 20)
	# the title shrinks to fit a narrow figure
	W = round(Int, 520*size_ratio*ncol)
	isempty(title) || Label(fig[1, 1:ncol], title; fontsize = min(labelsize + 6, round(Int, 1.8W / length(title))), font = :bold)
	for p in pans
		row, col = parse(Int, p[2]), parse(Int, p[3])
		c = load_panel(joinpath(dir, p[6]))
		gh = (:top in views(p) && :front in views(p)) ? GHOST : nothing
		for (jv, v) in enumerate(views(p))
			draw_cells!(fig[toff + start[row] + jv - 1, col], c; view = v, bs = pbs(p), ghost = gh,
						label = jv == 1 ? replace(p[5], "|" => "\n") : "", labelsize = labelsize)
		end
	end
	for (j, f) in enumerate(flat)
		f && rowsize!(fig.layout, toff + j, Fixed(round(Int, 0.4h)))
	end
	isempty(title) || rowgap!(fig.layout, 1, 50)
	fig
end;

"""
	write_panels(dir, title, panels; same_scale = false)

Write a final figure's data folder. `panels` is a vector of NamedTuples
(panel, row, col, views, label, cells, meta): cells has x, y, z, r, cell_state
(e.g. one frame of load_frames), meta a NamedTuple saved in the .jld2.
`same_scale` writes one box half-width (that of the largest panel) for every
panel, so plot_panels draws them at the same scale.
"""
function write_panels(dir, title, panels; same_scale = false)
	mkpath(dir)
	write(joinpath(dir, "title.txt"), title)
	bs = same_scale ? @sprintf("%.3f", maximum(box_halfwidth(p.cells) for p in panels)) : ""
	open(joinpath(dir, "panels.tsv"), "w") do io
		println(io, "panel\trow\tcol\tviews\tlabel\tfile\tbs")
		for p in panels
			c = p.cells
			f = save_panel(dir, p.panel, c.x, c.y, c.z, c.r, c.cell_state; p.meta...)
			println(io, join((p.panel, p.row, p.col, p.views, replace(p.label, "\n" => "|"), f, bs), '\t'))
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
	fig = Figure(size = (round(Int, 600*size_ratio*length(views)), round(Int, 600*size_ratio)), figure_padding = 20)
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
			elevation = v == :top ? pi/2 : 2*pi + pi/32,
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
free.csv, panels.tsv for the insets, readout.txt naming the class column to draw)
becomes four charts, figdir/phase_space*.png (y = plate gap, x = log10 N0),
*_N0.png (x = N0), *_diameter.png and *_diameter_N0.png (y = pancake diameter at
squeezing, when phase_space.csv has the column D).
"""
function replot_figures(figdir; save_png = true)
	figs = Dict{String,Any}()
	for name in sort(readdir(joinpath(figdir, "data")))
		dir = joinpath(figdir, "data", name)
		isdir(dir) || continue
		if startswith(name, "phase_space") && isfile(joinpath(dir, "phase_space.csv"))
			# which readout to draw: readout.txt names the class column (default class_02)
			ccol = isfile(joinpath(dir, "readout.txt")) ? String(strip(read(joinpath(dir, "readout.txt"), String))) : "class_02"
			raw, hdr = readdlm(joinpath(dir, "phase_space.csv"), ','; header = true)
			col = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
			# the same readout later (class_..._012), to break ties in the background
			lcol = replace(ccol, r"_02$" => "_012")
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
			for (suffix, x, y) in (("", :log, :layers), ("_N0", :linear, :layers),
								   ("_lin", :log, :layers_lin), ("_lin_N0", :linear, :layers_lin),
								   ("_diameter", :log, :diameter), ("_diameter_N0", :linear, :diameter))
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
