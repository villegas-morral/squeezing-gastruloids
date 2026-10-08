# Stage 6: the final figures (reports/FIGURES.md says how each one is made).
#
#   julia --project=. code-claude/scripts/stage6_figures.jl <cand> <core> <phase> [data|plot|movies|all]
#
#   cand, core  the chosen parameter set and initial core (e.g. july_frB5 sharp)
#   phase       the stage 5 results folder of its phase space (results/stage5/<phase>)
#
#   data    pick the runs and write code-claude/figures/data/<figure>/ (the cells of
#           every panel as .csv + .jld2, panels.tsv, title.txt; phase_space.csv)
#   plot    redraw every figure from figures/data alone (replot_figures) -> figures/*.png
#   movies  figures/movies/*.mp4, from the trajectories copied to figures/data/movies/
#   all     the three in turn (default)
#
# Nothing here simulates. The runs come from results/stage5/<phase>/ (confined) and
# results/stage5/tune/ (unconfined runs, fusion assay, the sharp/diluted comparison).

using Statistics, Printf, DelimitedFiles, JLD2, CairoMakie, LinearAlgebra
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))      # read_rows, with_readout!, frame_at

const CAND  = ARGS[1]
const CORE  = ARGS[2]
const PHASE = ARGS[3]
const MODE  = length(ARGS) >= 4 ? ARGS[4] : "all"
const S5    = joinpath(@__DIR__, "..", "results", "stage5")
const FIG   = joinpath(@__DIR__, "..", "figures")
const DATA  = joinpath(FIG, "data")
const HOURS = 15.0                       # time units per hour
const BASE  = String(split(CAND, "@")[1]) # "x@lw3" -> "x": unconfined runs and the fusion assay have no plates

fB_of(d, k) = count(==(2), d.cell_state[k]) / length(d.cell_state[k])
cells(d, k) = (x = d.x[k], y = d.y[k], z = d.z[k], r = d.r[k], cell_state = d.cell_state[k])

"""
Rotate a configuration so that the B centroid lies along +x from the centre: about z
only when confined (the slab stays flat), in 3D otherwise. Only for display, so a
polarised aggregate shows its B side in the same place in every panel.
"""
function align_B(c; confined = true, ref = c)
	xc, yc, zc = mean(ref.x), mean(ref.y), mean(ref.z)
	b = findall(==(2), ref.cell_state)
	isempty(b) && return c
	v = [mean(ref.x[b]) - xc, mean(ref.y[b]) - yc, confined ? 0.0 : mean(ref.z[b]) - zc]
	norm(v) < 1e-9 && return c
	v ./= norm(v)
	e = [1.0, 0.0, 0.0]
	ax = cross(v, e); s = norm(ax); co = dot(v, e)
	R = s < 1e-12 ? (co > 0 ? Matrix(1.0I, 3, 3) : Diagonal([-1.0, -1.0, 1.0])) :
		(K = [0 -ax[3] ax[2]; ax[3] 0 -ax[1]; -ax[2] ax[1] 0] ./ s;
		 Matrix(1.0I, 3, 3) + s .* K + (1 - co) .* (K * K))
	P = R * vcat((c.x .- xc)', (c.y .- yc)', (c.z .- zc)')
	(x = P[1, :] .+ xc, y = P[2, :] .+ yc, z = P[3, :] .+ zc, r = c.r, cell_state = c.cell_state)
end

pick(rows, kind, L, N0; seed = nothing, cand = kind == "free" ? BASE : CAND, core = CORE) =
	findfirst(r -> r["kind"] == kind && String(r["cand"]) == cand && String(r["core"]) == core &&
				   r["layers"] == L && r["N0"] == N0 && (seed === nothing || r["seed"] == seed), rows)

prow  = with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, PHASE))))
trow  = with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, "tune"))))
# unconfined runs: from tuning, and from the phase folder when it has some (a new core)
append!(trow, filter(r -> r["kind"] == "free", prow))
phase = with_diameter!(filter(r -> r["kind"] == "confined" && String(r["cand"]) == CAND && String(r["core"]) == CORE, prow))
println("$(length(phase)) phase-space runs of $CAND ($CORE)")
pcls(r) = phase_class(String(r["class_02"]))
major(L, N0) = (o = [pcls(r) for r in phase if r["layers"] == L && r["N0"] == N0];
				isempty(o) ? "" : argmax(c -> count(==(c), o), unique(o)))

"""
The example run of a class: at the first preferred point whose majority is that
class, its lowest seed with that class; if no preferred point qualifies, any run of
that class.
"""
function example(cls, preferred)
	for (L, N0) in preferred
		major(L, N0) == cls || continue
		rs = sort(filter(r -> r["layers"] == L && r["N0"] == N0 && pcls(r) == cls, phase); by = r -> r["seed"])
		isempty(rs) || return rs[1]
	end
	rs = filter(r -> pcls(r) == cls, phase)
	isempty(rs) ? nothing : rs[cld(length(rs), 2)]
end
load(r) = load_frames(joinpath(r["dir"], String(r["file"])))
lname(L) = @sprintf("%g layer%s", L, L == 1 ? "" : "s")
meta(r, d, k) = (run = String(r["file"]), cand = String(r["cand"]), core = String(r["core"]), layers = Float64(r["layers"]),
				 N0 = Int(r["N0"]), seed = Int(r["seed"]), frame = k, t = d.t[k], hours = d.t[k] / HOURS, fB = fB_of(d, k))

EX = Dict(
	"CENTRAL"   => example("CENTRAL",   [(4, 600), (3, 600), (2, 150), (2, 300), (3, 300), (4, 300)]),
	"PATTERN"   => example("PATTERN",   [(1, 600), (1, 1200), (2, 2400), (1, 2400)]),
	"POLARISED" => example("POLARISED", [(8, 300), (6, 150), (8, 150), (12, 600), (4, 150)]),
)
# a pattern in a slab of 3+ layers, if the set has one (no fallback: none if no preferred point qualifies)
let r = example("PATTERN", [(3, 2400), (4, 2400), (3, 4800), (4, 4800), (6, 2400), (3, 1200)])
	r !== nothing && r["layers"] >= 3 && (EX["PATTERN_THICK"] = r)
end
for (c, r) in EX
	println(@sprintf("%-10s example: ", c), r === nothing ? "none" : @sprintf("%s, N0 %d, seed %d (%s)",
			lname(r["layers"]), r["N0"], r["seed"], r["class_02"]))
end

"Pick the runs and write figures/data (every figure's panels, the phase space, the movie trajectories)."
function write_data()
	rm(DATA; force = true, recursive = true); mkpath(DATA)

	# ---- fig1: the initial condition, sharp and diluted -----------------------------
	pans = NamedTuple[]
	col = 0
	core_label(c) = c == "sharp" ? "sharp" : startswith(c, "fin") ? "diluted ($(round(Int, 100parse(Float64, c[4:end]))) % inside)" :
					 startswith(c, "cc") ? "C mixed in ($(round(Int, 100parse(Float64, c[3:end]))) % of the core)" : c
	for core in unique((CORE, "sharp", "fin0.75")), (kind, L, N0) in (("confined", 2, 2400), ("free", 0, 600), ("free", 0, 500))
		kind == "free" && N0 == 500 && any(p -> p.panel == "ic_$(core)_free", pans) && continue
		# the run of the chosen set (phase space, then tuning); for the diluted core,
		# which was only run in tuning, the run of the friction-free base set, or of the
		# Oct 7 set (the seeding does not depend on the B-B range; meta says which run)
		r = nothing
		for (rows, cand) in ((prow, kind == "free" ? BASE : CAND), (trow, kind == "free" ? BASE : CAND), (trow, BASE),
							 (trow, "july_frB4.5_frC6"))
			i = pick(rows, kind, L, N0; seed = 1, core = core, cand = cand)
			i === nothing || (r = rows[i]; break)
		end
		r === nothing && continue
		d = load(r); k = frame_at(d, 0.6); col += 1
		push!(pans, (panel = "ic_$(core)_$(kind)", row = 1, col = col, views = kind == "free" ? "top+front" : "top+side",
					 label = @sprintf("%s core|%s, %d cells|B = %.0f %%", core_label(core),
									  kind == "free" ? "unconfined" : lname(L), length(d.x[k]), 100fB_of(d, k)),
					 cells = cells(d, k), meta = meta(r, d, k)))
	end
	write_panels(joinpath(DATA, "fig1_initial_condition"),
				 "Initial condition: B (green) = 60 % of the cells, C (orange) outside", pans)

	# ---- fig2: one aggregate size under no / medium / strong confinement, through time ---------
	# (your presentation sketch, panel B): N0 600, unconfined (polarises), 4 layers (central
	# cluster that shrinks), 1 layer (islands)
	pans = NamedTuple[]
	fBs = (0.6, 0.45, 0.3, 0.2, 0.12)
	rows_ = Any[]
	i = pick(trow, "free", 0, 600; seed = 1)
	i === nothing || push!(rows_, (trow[i], "No confinement", false))
	EX["CENTRAL"] === nothing || push!(rows_, (EX["CENTRAL"], "Medium confinement", true))
	EX["PATTERN"] === nothing || push!(rows_, (EX["PATTERN"], "Strong confinement", true))
	for (row, (r, name, confined)) in enumerate(rows_)
		d = load(r); kend = frame_at(d, 0.12)
		for (j, f) in enumerate(fBs)
			k = frame_at(d, f)
			c = confined ? cells(d, k) : align_B(cells(d, k); confined = false, ref = cells(d, kend))
			where_ = confined ? "$(lname(r["layers"])), N0 $(r["N0"])" : "unconfined, N0 $(r["N0"])"
			lbl = (j == 1 ? "$name: $where_|" : "|") *
				  @sprintf("t = %.0f h, fB = %.2f%s", d.t[k] / HOURS, fB_of(d, k), f == 0.2 ? "  (readout)" : "")
			push!(pans, (panel = "row$(row)_$(j)", row = row, col = j, views = confined ? "top+side" : "top+front",
						 label = lbl, cells = c, meta = meta(r, d, k)))
		end
	end
	write_panels(joinpath(DATA, "fig2_behaviours"),
				 "No, medium and strong confinement (top: mid-plane cutaway; below: side)", pans)

	# ---- fig3: the patterns at the readout --------------------------------------------
	# every thickness that has patterns: its largest N0 first, then its smallest, up to 8
	pats = sort(filter(r -> pcls(r) == "PATTERN" && r["seed"] == 1 && major(r["layers"], r["N0"]) == "PATTERN", phase);
				by = r -> (r["layers"], r["N0"]))
	byL = [filter(r -> r["layers"] == L, pats) for L in sort(unique(r["layers"] for r in pats))]
	pats = unique(vcat([last(g) for g in byL], [first(g) for g in byL], pats))
	pats = sort(pats[1:min(end, 8)]; by = r -> (r["layers"], r["N0"]))
	pans = NamedTuple[]
	for (i, r) in enumerate(pats[1:min(end, 8)])
		d = load(r); k = frame_at(d, 0.2)
		st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k])
		push!(pans, (panel = @sprintf("L%g_N%d", r["layers"], r["N0"]), row = 1 + (i - 1) ÷ 4, col = 1 + (i - 1) % 4,
					 views = "top", label = @sprintf("%s, N0 %d|%s: %d islands, S = %.2f", lname(r["layers"]), r["N0"],
											lowercase(String(r["class_02"])), st.n_isl, st.S),
					 cells = cells(d, k), meta = meta(r, d, k)))
	end
	write_panels(joinpath(DATA, "fig3_patterns"),
				 "B patterns at the readout fB = 0.2 (seed 1, mid-plane cutaway)", pans)

	# ---- fig4: unconfined aggregates of every size, at the readout ---------------------------
	# (the "free" row of the phase space): B polarises up to N0 ~ 600, bigger ones are
	# still central when B runs out
	pans = NamedTuple[]
	for (j, N0) in enumerate((150, 300, 600, 1200, 2400))
		i = pick(trow, "free", 0, N0; seed = 1)
		i === nothing && continue
		r = trow[i]; d = load(r); k = frame_at(d, 0.2)
		st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = false)
		push!(pans, (panel = "free_N$(N0)", row = 1, col = j, views = "top+front",
					 label = @sprintf("N0 %d, fB = %.2f|%s, B offset %.2f", N0, fB_of(d, k), lowercase(classify_pattern(st)), st.offset),
					 cells = align_B(cells(d, k); confined = false), meta = meta(r, d, k)))
	end
	write_panels(joinpath(DATA, "fig4_unconfined"),
				 "Unconfined aggregates at the readout fB = 0.2 (top: mid-plane cut; below: whole)", pans)

	# ---- fig5: the fusion assay -------------------------------------------------------
	ff = joinpath(S5, "tune", "fusion_$(BASE).jld2")
	if isfile(ff)
		pans = NamedTuple[]
		jldopen(ff, "r") do fh
			m = fh["metric"]
			for (j, pr) in enumerate(("BB", "CC", "BC"))
				c = fh[pr]
				push!(pans, (panel = "doublet_$(pr)", row = 1, col = j, views = "top+front",
							 label = @sprintf("%s doublet, 24 h|separation / (Rg1 + Rg2) = %.2f", pr, m[pr][3]),
							 cells = c, meta = (assay = "fusion_assay(n = 150, hours = 24)", pair = pr,
												t0 = m[pr][1], h12 = m[pr][2], h24 = m[pr][3], cand = CAND)))
			end
		end
		write_panels(joinpath(DATA, "fig5_fusion"),
					 "Fusion: two 150-cell aggregates after 24 h (1.6 touching, 0.5 one sphere)", pans)
	end

	# ---- phase space: at the readout fB = 0.2, and the same runs later (fB = 0.12) ------------
	ns = [count(r -> (r["layers"], r["N0"]) == p, phase) for p in unique((r["layers"], r["N0"]) for r in phase)]
	nseed_txt = minimum(ns) == maximum(ns) ? "$(minimum(ns)) seeds" : "$(minimum(ns))–$(maximum(ns)) seeds per point"
	for (folder, ccol, fread, ttl) in (("phase_space", "class_02", 0.2, "Phase space: outcome at the readout fB = 0.2 ($nseed_txt)"),
									   ("phase_space_late", "class_012", 0.12, "The same runs later, at fB = 0.12"))
		dir = joinpath(DATA, folder); mkpath(dir)
		write(joinpath(dir, "readout.txt"), ccol)
		write(joinpath(dir, "title.txt"), ttl)
		keys_ = ["layers", "N0", "seed", "D", "class_04", "class_03", "class_02", "class_015", "class_012", "t_02",
				 "n_isl", "f_isl", "f_iso", "f_max", "S", "spread", "fill", "offset", "offset_max", "layers_eff", "N_end", "file"]
		open(joinpath(dir, "phase_space.csv"), "w") do io
			println(io, join(keys_, ","))
			for r in sort(phase; by = r -> (r["layers"], r["N0"], r["seed"]))
				println(io, join([string(r[k]) for k in keys_], ","))
			end
		end
		# the unconfined runs of the same mechanics (no plates): the "free" row of the chart
		open(joinpath(dir, "free.csv"), "w") do io
			println(io, "N0,seed,class_02,class_012,offset,file")
			grid_N0 = Set(r["N0"] for r in phase)      # only the N0 values of the grid
			for r in sort(filter(r -> r["kind"] == "free" && String(r["cand"]) == BASE && String(r["core"]) == CORE &&
								 r["N0"] in grid_N0, trow);
						  by = r -> (r["N0"], r["seed"]))
				println(io, join((r["N0"], r["seed"], r["class_02"], r["class_012"], @sprintf("%.4f", r["offset"]), r["file"]), ","))
			end
		end
		# insets: the examples of fig2 at this readout
		open(joinpath(dir, "panels.tsv"), "w") do io
			println(io, "panel\tlayers\tN0\tviews\tlabel\tfile")
			for (cls, name) in (("CENTRAL", "central cluster"), ("PATTERN", "pattern"), ("PATTERN_THICK", "pattern"), ("POLARISED", "polarised"))
				r = get(EX, cls, nothing); r === nothing && continue
				d = load(r); k = frame_at(d, fread)
				c = cls == "POLARISED" ? align_B(cells(d, k)) : cells(d, k)
				f = save_panel(dir, "inset_$(lowercase(cls))", c.x, c.y, c.z, c.r, c.cell_state; meta(r, d, k)...)
				println(io, join((lowercase(cls), r["layers"], r["N0"], "top+side",
								  "$name|$(lname(r["layers"])), N0 $(r["N0"]), fB $(fread)", f), '\t'))
			end
		end
	end

	# ---- trajectories for the movies ---------------------------------------------------
	# from seeding to the end only (the growth frames are most of the file), in the
	# save_frames format, so load_frames reads them
	mdir = joinpath(DATA, "movies"); mkpath(mdir)
	function keep_trajectory(src, dst)
		d = load_frames(src); ks = frame_at(d, 0.6):length(d.x)
		jldsave(dst; x = d.x[ks], y = d.y[ks], z = d.z[ks], r = d.r[ks], cell_state = d.cell_state[ks],
				t = d.t[ks], N = d.N[ks], source = basename(src))
	end
	for (cls, r) in EX
		r === nothing && continue
		keep_trajectory(joinpath(r["dir"], String(r["file"])), joinpath(mdir, "$(lowercase(cls)).jld2"))
	end
	i = pick(trow, "free", 0, 600; seed = 1)     # the unconfined row of fig2
	i === nothing || keep_trajectory(joinpath(trow[i]["dir"], String(trow[i]["file"])), joinpath(mdir, "unconfined.jld2"))
	println("wrote ", DATA)
end

"The movies, from the trajectories kept in figures/data/movies."
function make_movies()
	mdir = joinpath(DATA, "movies"); odir = joinpath(FIG, "movies"); mkpath(odir)
	for f in sort(readdir(mdir))
		endswith(f, ".jld2") || continue
		d = load_frames(joinpath(mdir, f))
		ks = frame_at(d, 0.6):frame_at(d, 0.1)
		record_pancake_frames(d, color_map, ks; fps = 6, filename = joinpath(odir, replace(f, ".jld2" => ".mp4")))
		println("wrote movies/", replace(f, ".jld2" => ".mp4"))
	end
end

MODE in ("data", "all")   && write_data()
MODE in ("plot", "all")   && replot_figures(FIG)
MODE in ("movies", "all") && make_movies()
