# Stage 7: the presentation figure, laid out like your sketch
# (notes/3d-simulations presentation idea.jpg), as an A4 PDF to finish in Inkscape.
#
#   julia --project=. code-claude/scripts/stage7_presentation.jl [figdir]
#
# figdir defaults to code-claude/figures. Nothing is simulated or picked here: the
# figure is drawn from the data that stage6_figures.jl wrote,
#
#   data/fig2_behaviours/   panel B (the same runs, frames at fB 0.6, 0.3, 0.2)
#   data/phase_space/       panel C (phase_space.csv, free.csv)
#
# and written to figdir/presentation.pdf (+ presentation.png, a preview).
#   A  model ingredients, drawn as vector shapes
#   B  no / medium / strong confinement through time (3D, each panel one embedded image)
#   C  the phase space, log10 N0 against the plate gap, in the colours of the sketch,
#      with the three outcomes as 3D insets
# Text stays text and every panel is its own object, so all of it can be moved and
# edited in Inkscape.

using Statistics, Printf, DelimitedFiles, JLD2, CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))

const FIG  = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "figures")
const DATA = joinpath(FIG, "data")

# the sketch's colours: T+ (B) green, T- (C) orange; regions grey / cyan / pink
const CELL = Dict(2 => Makie.to_color("#3fa535"), 3 => Makie.to_color("#f5b041"))
const REGION = Dict(
	"POLARISED" => Makie.to_color("#a9a9a9"),   # grey
	"CENTRAL"   => Makie.to_color("#1fb3c8"),   # cyan
	"PATTERN"   => Makie.to_color("#ee82b4"),   # pink
	"OTHER"     => Makie.to_color("#e6e6e6"))
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
One 3D panel in the sketch colours, rasterised (one image in the PDF).
:top     seen from above, cells above the mid-plane left out
:oblique a three-quarter view, cells above the mid-plane left out (the B pattern shows
         on the cut face, the slab keeps its thickness)
"""
function cells3d!(pos, c; view = :oblique, bs = nothing, kw...)
	x, y, z, r, s = c.x, c.y, c.z, c.r, c.cell_state
	xc, yc, zc = mean(x), mean(y), mean(z)
	bs = bs === nothing ? max(maximum(x) - minimum(x), maximum(y) - minimum(y)) / 2 + 1.5 : bs
	keep = findall(w -> w <= zc + 0.2, z)
	ax = Axis3(pos; aspect = :data, azimuth = view == :top ? 0.0 : -pi/2 + 0.35,
			   elevation = view == :top ? pi/2 : pi/6, protrusions = 0, viewmode = :fit, kw...)
	hidedecorations!(ax); hidespines!(ax)
	ax.xspinesvisible = false; ax.yspinesvisible = false; ax.zspinesvisible = false
	o = sort(keep; by = m -> s[m], rev = true)
	meshscatter!(ax, x[o] .- xc, y[o] .- yc, z[o] .- zc; markersize = r[o], color = [CELL[s[m]] for m in o],
				 rasterize = 4)
	# a flat slab gets a flat box, so it fills the panel
	zl = view == :oblique ? (minimum(z) - zc - 1.5, maximum(z) - zc + 1.5) : (-bs, bs)
	xlims!(ax, -bs, bs); ylims!(ax, -bs, bs); zlims!(ax, zl...)
	ax
end

"A darker shade of a colour, for text on the region colours."
darker(c) = RGBf(0.55c.r, 0.55c.g, 0.55c.b)

# ---- data -----------------------------------------------------------------------------

b_dir = joinpath(DATA, "fig2_behaviours")
pans = [split(l, '\t') for l in readlines(joinpath(b_dir, "panels.tsv"))[2:end] if !isempty(strip(l))]
panel(row, col) = (p = only(filter(p -> p[2] == string(row) && p[3] == string(col), pans));
				   (cells = load_panel(joinpath(b_dir, p[6])), meta = jldopen(f -> f["meta"], joinpath(b_dir, replace(p[6], ".csv" => ".jld2")))))
COLS = (1, 3, 4)                     # fig2 columns: fB 0.6 (seeding), 0.3, 0.2 (readout)
ROWS = ((1, "no\nconfinement", :top), (2, "medium\nconfinement", :oblique), (3, "strong\nconfinement", :oblique))

p_dir = joinpath(DATA, "phase_space")
raw, hdr = readdlm(joinpath(p_dir, "phase_space.csv"), ','; header = true)
col = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
rows = [(layers = Float64(raw[i, col["layers"]]), N0 = Int(raw[i, col["N0"]]), class = String(raw[i, col["class_02"]]))
		for i in 1:size(raw, 1)]
raw, hdr = readdlm(joinpath(p_dir, "free.csv"), ','; header = true)
fc = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
free = [(N0 = Int(raw[i, fc["N0"]]), class = String(raw[i, fc["class_02"]])) for i in 1:size(raw, 1)]

# ---- the page: A4 in points, everything placed by hand (bbox = left, right, bottom, top) ----

W, H = 595, 842
fig = Figure(size = (W, H), fontsize = 10, figure_padding = 0)
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
text!(ax, -0.95, -0.25; text = "T+", align = (:center, :top), fontsize = 11)
text!(ax, 0.95, -0.25; text = "T−", align = (:center, :top), fontsize = 11)
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
		cells3d!(fig, p.cells; view = view, bbox = BBox(colx(j), colx(j) + 140, rowy(i), rowy(i) + 82))
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
	fB = (0.6, 0.3, 0.2)[j]
	text!(page, x, yt - 7; text = @sprintf("%.0f h   (T+ %.0f %%)", hours[j], 100fB), align = (:center, :top), fontsize = 9)
end
text!(page, W - 18, yt + 5; text = "time", align = (:right, :bottom), fontsize = 11)

# C: the phase space (y 40-390), insets on the right
letter(28, 382, "C")
axC = phase_space!(fig, rows; free = free, x = :log, y = :layers, colors = REGION, background = 0.45,
				   fontsize = 12, pie_scale = 0.5, region_labels = REGION_TEXT,
				   axis = (bbox = BBox(80, 420, 50, 372),))
axC.xlabel = "log₁₀ N0 (number of cells at squeezing)"; axC.ylabel = "number of layers"
for (k, (pnl, name)) in enumerate((((1, 3), "polarisation (unconfined)"), ((2, 3), "no polarisation"), ((3, 3), "patterns")))
	y0 = 262 - (k - 1) * 105
	cells3d!(fig, panel(pnl...).cells; view = :top, bbox = BBox(445, 565, y0, y0 + 95))
	text!(page, 505, y0 + 97; text = name, align = (:center, :bottom), fontsize = 9, color = darker(REGION[("POLARISED", "CENTRAL", "PATTERN")[k]]))
end

save(joinpath(FIG, "presentation.pdf"), fig; pt_per_unit = 1)
save(joinpath(FIG, "presentation.png"), fig; px_per_unit = 3)
println("wrote ", joinpath(FIG, "presentation.pdf"))
