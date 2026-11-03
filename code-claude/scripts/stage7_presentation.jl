# Stage 7: the presentation figure, laid out like your sketch
# (notes/3d-simulations presentation idea.jpg), as an A4 PDF to finish in Inkscape.
#
#   julia --project=. code-claude/scripts/stage7_presentation.jl [figdir]
#
# figdir defaults to code-claude/figures. Nothing is simulated or picked here: the
# figure is drawn from the data that stage6_figures.jl wrote,
#
#   data/fig2_behaviours/   panel B (the same runs, frames at fB 0.6, 0.3, 0.12)
#   data/phase_space/       panel C (phase_space.csv, free.csv)
#
# and written to figdir/presentation.pdf (+ presentation.png, a preview).
#   A  model ingredients, drawn as vector shapes
#   B  no / medium / strong confinement through time (3D, each panel one embedded image)
#   C  the phase space, log10 N0 against the plate gap, in the colours of the sketch,
#      with the iso-area lines and each outcome as a 3D inset on its own region
# Text stays text and every panel is its own object, so all of it can be moved and
# edited in Inkscape.

using Statistics, Printf, DelimitedFiles, JLD2
using CairoMakie                 # the page (vector, editable in Inkscape)
import GLMakie                   # the 3D panels, rendered as images (translucent C, cells3d!)
# GLMakie renders in a small hidden window (a visible one fails without a desktop session, and
# GLMakie caps the window at the screen size); the page itself is CairoMakie
GLMakie.activate!(visible = false, scalefactor = 0.3, px_per_unit = 1)
CairoMakie.activate!()
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))

const FIG  = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "figures")
const DATA = joinpath(FIG, "data")

# the sketch's colours: T+ (B) green, T- (C) yellow (the C of every other figure);
# regions grey / cyan / pink (PHASE_COLORS, the same as the report's phase space)
const CELL = Dict(2 => Makie.to_color("#3fa535"), 3 => color_map[3])
const REGION = PHASE_COLORS
const REGION_TEXT = Dict("POLARISED" => "polarisation", "CENTRAL" => "no polarisation", "PATTERN" => "patterns")

# ---- small drawing helpers ------------------------------------------------------------

"An arrow from p to q (2D axis coordinates), with a filled head."
function arrow!(ax, p, q; color = :black, linewidth = 1.5, head = 0.18)
	d = (q[1] - p[1], q[2] - p[2]); n = hypot(d...); u = (d[1] / n, d[2] / n); w = (-u[2], u[1])
	b = (q[1] - head * u[1], q[2] - head * u[2])
	lines!(ax, [p[1], b[1]], [p[2], b[2]]; color = color, linewidth = linewidth)
	poly!(ax, Point2f[q, (b[1] + 0.45head * w[1], b[2] + 0.45head * w[2]), (b[1] - 0.45head * w[1], b[2] - 0.45head * w[2])];
		  color = color)
end

"A cell: a circle with a light hatch, as in the sketch (or filled with `fill`)."
function cell!(ax, c; r = 0.42, fill = nothing)
	th = range(0, 2pi, length = 60)
	poly!(ax, Point2f[(c[1] + r*cos(t), c[2] + r*sin(t)) for t in th];
		  color = fill === nothing ? RGBf(0.93, 0.93, 0.93) : fill, strokecolor = :black, strokewidth = 1.2)
	if fill === nothing
		for s in -0.6:0.3:0.6
			h = sqrt(max(0, 1 - s^2))
			lines!(ax, [c[1] + r*(s - 0.25h), c[1] + r*(s + 0.25h)], [c[2] - r*h*0.9, c[2] + r*h*0.9];
				   color = RGBf(0.45, 0.45, 0.45), linewidth = 0.6)
		end
	end
end

"A blank 2D axis for drawing in (data units, equal aspect)."
function sketch_axis(pos; title = "", lims = (-1.6, 1.6, -1.3, 1.3), kw...)
	ax = Axis(pos; aspect = DataAspect(), title = title, titlesize = 11, titlefont = :regular,
			  backgroundcolor = :transparent, kw...)
	hidedecorations!(ax); hidespines!(ax)
	limits!(ax, lims...)
	ax
end

"""
One 3D panel in the sketch colours, as one image in the PDF at `bbox`.
:top     seen from above
:oblique a three-quarter view
ghost = nothing: opaque, the cells above the mid-plane left out (the cutaway of the
confined aggregates); ghost = GHOST: the whole aggregate, C translucent (unconfined).
As in draw_cells!. The panel is rendered with GLMakie at `scale` pixels per point and
placed as an image, its white background made transparent; the rest of the page stays vector
(CairoMakie). `bs` is the half-width of the box: the same for every panel = the same scale.
"""
function cells3d!(fig, c; view = :oblique, bs = nothing, bbox, scale = 4, ghost = nothing)
	x, y, z, r, s = c.x, c.y, c.z, c.r, c.cell_state
	xc, yc, zc = mean(x), mean(y), mean(z)
	if ghost === nothing                        # the cutaway
		k = findall(w -> w <= zc + 0.2, z)
		x, y, z, r, s = x[k], y[k], z[k], r[k], s[k]
	end
	bs = bs === nothing ? max(maximum(x) - minimum(x), maximum(y) - minimum(y)) / 2 + 1.5 : bs
	w, h = round.(Int, scale .* Makie.widths(bbox))
	g = Figure(size = (w, h), figure_padding = 0, backgroundcolor = :white)
	ax = Axis3(g[1, 1]; aspect = :data, azimuth = view == :top ? 0.0 : -pi/2 + 0.35,
			   elevation = view == :top ? pi/2 : pi/6, protrusions = 0, viewmode = :fit)
	hidedecorations!(ax); hidespines!(ax)
	ax.xspinesvisible = false; ax.yspinesvisible = false; ax.zspinesvisible = false
	cells_mesh!(ax, x .- xc, y .- yc, z .- zc, r, s; ghost = ghost, cmap = CELL)
	# a flat slab gets a flat box, so it fills the panel
	zl = view == :oblique ? (minimum(z) - zc - 1.5, maximum(z) - zc + 1.5) : (-bs, bs)
	xlims!(ax, -bs, bs); ylims!(ax, -bs, bs); zlims!(ax, zl...)
	f = tempname() * ".png"
	save(f, g; backend = GLMakie)
	axi = Axis(fig; bbox = bbox, aspect = DataAspect(), backgroundcolor = :transparent)
	hidedecorations!(axi); hidespines!(axi)
	# GLMakie cannot draw on a transparent background: the white is made transparent here,
	# so a panel placed on the phase space shows the region's colour around the cells
	img = [min(p.r, p.g, p.b) > 0.99f0 ? RGBAf(1, 1, 1, 0) : RGBAf(p.r, p.g, p.b, 1) for p in Makie.FileIO.load(f)]
	image!(axi, rotr90(img))
	axi
end

# ---- data -----------------------------------------------------------------------------

b_dir = joinpath(DATA, "fig2_behaviours")
pans = [split(l, '\t') for l in readlines(joinpath(b_dir, "panels.tsv"))[2:end] if !isempty(strip(l))]
panel(row, col) = (p = only(filter(p -> p[2] == string(row) && p[3] == string(col), pans));
				   (cells = load_panel(joinpath(b_dir, p[6])), meta = jldopen(f -> f["meta"], joinpath(b_dir, replace(p[6], ".csv" => ".jld2")))))
COLS = (1, 3, 5)                     # fig2 columns: fB 0.6 (seeding), 0.3, 0.12 (readout)
COL_FB = (0.6, 0.3, 0.12)
ROWS = ((1, "no\nconfinement", :top), (2, "medium\nconfinement", :oblique), (3, "strong\nconfinement", :oblique))

p_dir = joinpath(DATA, "phase_space")
# the readout the phase-space figure draws (readout.txt: the class column, then the one that
# breaks ties between seeds; class_02 for the October figures)
rl = filter(!isempty, strip.(readlines(joinpath(p_dir, "readout.txt"))))
ccol = String(rl[1])
raw, hdr = readdlm(joinpath(p_dir, "phase_space.csv"), ','; header = true)
col = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
lcol = length(rl) >= 2 ? String(rl[2]) : replace(ccol, r"_02$" => "_012")
rows = [(layers = Float64(raw[i, col["layers"]]), N0 = Int(raw[i, col["N0"]]), class = String(raw[i, col[ccol]]),
		 late = haskey(col, lcol) ? String(raw[i, col[lcol]]) : String(raw[i, col[ccol]])) for i in 1:size(raw, 1)]
# the pancake area at squeezing, for the iso-area lines
cont = haskey(col, "A_um2") ? [(layers = r.layers, N0 = r.N0, A = Float64(raw[i, col["A_um2"]])) for (i, r) in enumerate(rows)] : nothing
raw, hdr = readdlm(joinpath(p_dir, "free.csv"), ','; header = true)
fc = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
free = [(N0 = Int(raw[i, fc["N0"]]), class = String(raw[i, fc[ccol]]),
		 late = haskey(fc, lcol) ? String(raw[i, fc[lcol]]) : String(raw[i, fc[ccol]])) for i in 1:size(raw, 1)]

# ---- the page: A4 in points, everything placed by hand (bbox = left, right, bottom, top) ----

W, H = 595, 842
fig = Figure(size = (W, H), fontsize = 10, figure_padding = 0, fonts = TEX_FONTS)
# a transparent page-sized axis for the free text and lines
page = Axis(fig; bbox = BBox(0, W, 0, H), backgroundcolor = :transparent)
hidedecorations!(page); hidespines!(page); limits!(page, 0, W, 0, H)
letter(x, y, t) = text!(page, x, y; text = t, fontsize = 20, font = :bold, align = (:left, :top))

# A: the model ingredients (y 705-815)
letter(28, 818, "A")
ax = sketch_axis(fig; title = "adhesion", bbox = BBox(60, 180, 705, 805))
cell!(ax, (-0.55, 0.15)); cell!(ax, (0.55, 0.15))
arrow!(ax, (-0.95, 0.15), (-0.2, 0.15); color = :orangered, head = 0.2)
arrow!(ax, (0.95, 0.15), (0.2, 0.15); color = :orangered, head = 0.2)
lines!(ax, [-0.55, 0.85], [-0.62, -0.62]; color = :black, linewidth = 1)
lines!(ax, [-0.55, -0.55], [-0.52, -0.72]; color = :black, linewidth = 1)
lines!(ax, [0.85, 0.85], [-0.52, -0.72]; color = :black, linewidth = 1)
text!(ax, 0.15, -0.75; text = "μ", align = (:center, :top), fontsize = 13)
ax = sketch_axis(fig; title = "proliferation", bbox = BBox(190, 310, 705, 805))
cell!(ax, (0, 0.7); r = 0.36)
arrow!(ax, (0, 0.28), (0, -0.12); head = 0.16)
cell!(ax, (-0.38, -0.55); r = 0.36); cell!(ax, (0.38, -0.55); r = 0.36)
text!(ax, 0.25, 0.1; text = "τ_div", align = (:left, :center), fontsize = 11)
ax = sketch_axis(fig; title = "differentiation", bbox = BBox(320, 440, 705, 805))
cell!(ax, (-0.95, 0.2); r = 0.32, fill = CELL[2]); cell!(ax, (0.95, 0.2); r = 0.32, fill = CELL[3])
arrow!(ax, (-0.5, 0.2), (0.5, 0.2); head = 0.18)
text!(ax, 0, 0.35; text = "q", align = (:center, :bottom), fontsize = 13)
text!(ax, -0.95, -0.25; text = "B (T+)", align = (:center, :top), fontsize = 11)
text!(ax, 0.95, -0.25; text = "C (T−)", align = (:center, :top), fontsize = 11)
ax = sketch_axis(fig; title = "friction", bbox = BBox(450, 570, 705, 805))
cell!(ax, (-0.45, -0.1); r = 0.32); cell!(ax, (0.45, 0.1); r = 0.32)
arrow!(ax, (-0.3, 0.12), (0.05, 0.75); head = 0.16); text!(ax, 0.12, 0.75; text = "v", fontsize = 11)
arrow!(ax, (-0.62, -0.38), (-0.98, -1.0); color = :orangered, head = 0.16)
arrow!(ax, (0.3, -0.12), (-0.05, -0.75); head = 0.16); text!(ax, 0.05, -0.85; text = "v", fontsize = 11, align = (:left, :top))
arrow!(ax, (0.62, 0.38), (0.98, 1.0); color = :orangered, head = 0.16)

# B: one aggregate under no / medium / strong confinement, through time (y 415-690)
letter(28, 690, "B")
hours = Float64[]
colx(j) = 150 + (j - 1) * 145                  # left edge of column j (width 140)
rowy(i) = 600 - (i - 1) * 85                   # bottom edge of row i (height 82)
for (i, (row, name, view)) in enumerate(ROWS)
	text!(page, 140, rowy(i) + 41; text = name, align = (:right, :center), justification = :right, fontsize = 10)
	for (j, c) in enumerate(COLS)
		p = panel(row, c)
		cells3d!(fig, p.cells; view = view, bbox = BBox(colx(j), colx(j) + 140, rowy(i), rowy(i) + 82),
				 ghost = row == 1 ? GHOST : nothing)     # C translucent in the unconfined row only
		i == 2 && push!(hours, p.meta["hours"])
	end
end
# the time axis under the columns, at the frames' times (medium row; B falls at the same
# rate in every run, so the other rows match within an hour or two)
yt = rowy(length(ROWS)) - 18
arrow!(page, (colx(1) - 5, yt), (W - 18, yt); head = 8, linewidth = 1.2)
for j in eachindex(COLS)
	x = colx(j) + 70
	lines!(page, [x, x], [yt - 4, yt + 4]; color = :black, linewidth = 1.2)
	fB = COL_FB[j]
	text!(page, x, yt - 7; text = @sprintf("%.0f h   (T+ %.0f %%)", hours[j], 100fB), align = (:center, :top), fontsize = 9)
end
text!(page, W - 18, yt + 5; text = "time", align = (:right, :bottom), fontsize = 11)

# C: the phase space (y 40-390), each example on its own region
letter(28, 382, "C")
# the same chart as the report's phase space (November): the regions' colours and the
# iso-area lines; stretched to the width of the page (no equal aspect)
axC = phase_space!(fig, rows; free = free, x = :log, y = :layers_lin, colors = REGION, background = 0.5,
				   fontsize = 12, contours = cont, contour_width = 1.4, contour_labelsize = 7,
				   axis = (bbox = BBox(80, 575, 50, 372), aspect = nothing))
# the readout frames of the rows of B (fig2 column 5), each titled with its outcome in the
# colour of its region and placed on that region, at (u, v) in the chart's coordinates
# (N0 across 0-10, the gap up 0-8, the unconfined row at 9.4; phase_space!)
Makie.update_state_before_display!(fig)
vp = axC.scene.viewport[]
lim = axC.finallimits[]
topage(u, v) = (vp.origin[1] + (u - lim.origin[1]) / lim.widths[1] * vp.widths[1],
				vp.origin[2] + (v - lim.origin[2]) / lim.widths[2] * vp.widths[2])
for (pnl, cls, (u, v)) in (((1, 5), "POLARISED", (1.8, 6.0)), ((2, 5), "CENTRAL", (8.7, 6.0)), ((3, 5), "PATTERN", (8.6, 1.0)))
	cx, cy = topage(u, v)
	cells3d!(fig, panel(pnl...).cells; view = :top, bbox = BBox(cx - 60, cx + 60, cy - 50, cy + 50),
			 ghost = cls == "POLARISED" ? GHOST : nothing)
	text!(page, cx, cy + 30; text = PHASE_NAMES[cls], align = (:center, :bottom), fontsize = 11, font = :bold,
		  color = darker(REGION[cls]))
end

save(joinpath(FIG, "presentation.pdf"), fig; pt_per_unit = 1, backend = CairoMakie)
save(joinpath(FIG, "presentation.png"), fig; px_per_unit = 3, backend = CairoMakie)
println("wrote ", joinpath(FIG, "presentation.pdf"))
