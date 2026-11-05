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

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using GLMakie                            # the translucent C cells of unconfined aggregates (GHOST); PNG output only
GLMakie.activate!(visible = false, scalefactor = 0.3, px_per_unit = 1)   # small hidden window: GLMakie caps the window at the screen size
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
# the 25 time units of relaxation after seeding have no division and no differentiation:
# times in the figures count from the start of differentiation (November 2026)
const T_RELAX = 25.0
hrs(d, k) = max(0.0, d.t[k] - T_RELAX) / HOURS
const BASE  = String(split(CAND, "@")[1]) # "x@lw3" -> "x": unconfined runs and the fusion assay have no plates
# the readout drawn in the figures (November 2026): the experiments' clustering and
# polarisation (classify_exp) with a minimum cluster of 10 cells instead of their 20,
# since our aggregates are 10-100 times smaller (with 20, monolayers of N0 <= 300 have
# no cluster at all), and a core cell needs half its neighbours B instead of their
# 40 % (CORE_FRAC, stage5_common.jl: rings held together by thin B threads count as
# clusters). Since Nov 3 it is read at fB = 0.12 (READ_FB; the experiments read at 96 h, with
# ~12 % T+), with three outcomes only: a run with no B cluster left there takes its last
# outcome (class_read, stage5_common.jl). The same readout at fB 0.2 (class_exp50_02) and
# the October classifier (class_02) are drawn too, for comparison.
const CCOL  = "class_read"
const CCOL_LATE = "class_read_late"
const READ_FB = 0.12

fB_of(d, k) = count(==(2), d.cell_state[k]) / length(d.cell_state[k])
cells(d, k) = (x = d.x[k], y = d.y[k], z = d.z[k], r = d.r[k], cell_state = d.cell_state[k])

"""
The frame the outcome is read from (class_read, stage5_common.jl): fB 0.12, else 0.15, as long as
at least half of the B cells are in clusters there (READ_MIN_CL); else fB 0.2.
"""
function read_frame(d)
	for f in (READ_FB, 0.15)
		k = frame_at(d, f)
		e = exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; nmin = 10, frac = CORE_FRAC)
		e.n_cl > 0 && e.f_cl >= READ_MIN_CL && return k
	end
	frame_at(d, 0.2)
end

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
pcls(r) = phase_class(String(r[CCOL]))
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
				 N0 = Int(r["N0"]), seed = Int(r["seed"]), frame = k, t = d.t[k], hours = hrs(d, k), fB = fB_of(d, k))

EX = Dict(
	"CENTRAL"   => example("CENTRAL",   [(4, 600), (3, 600), (2, 150), (2, 300), (3, 300), (4, 300)]),
	"PATTERN"   => example("PATTERN",   [(1, 600), (1, 1200), (2, 2400), (1, 2400)]),
	"POLARISED" => example("POLARISED", [(8, 300), (6, 150), (8, 150), (12, 600), (4, 150)]),
)
# a pattern in a slab of 3+ layers, if the set has one (no fallback: none if no preferred point qualifies)
let r = example("PATTERN", [(3, 2400), (4, 2400), (3, 4800), (4, 4800), (6, 2400), (3, 1200)])
	r !== nothing && r["layers"] >= 3 && (EX["PATTERN_THICK"] = r)
end
# the phase-space insets (November): a clearly central core, and clusters between a broken
# ring and fine islands
EX["CENTRAL_INSET"] = example("CENTRAL", [(3, 1200), (4, 1200), (3, 600), (4, 600)])
let r = example("PATTERN", [(2, 1200), (2, 2400), (1, 600)])
	r === nothing || (EX["PATTERN_INSET"] = r)
end
for (c, r) in EX
	println(@sprintf("%-10s example: ", c), r === nothing ? "none" : @sprintf("%s, N0 %d, seed %d (%s)",
			lname(r["layers"]), r["N0"], r["seed"], r[CCOL]))
end

"Pick the runs and write figures/data (every figure's panels, the phase space, the movie trajectories)."
function write_data()
	rm(DATA; force = true, recursive = true); mkpath(DATA)

	# ---- fig1: the initial condition (November: the core that is used, one N0, one scale, no titles)
	# 2 layers and unconfined, both N0 600 (as in fig2 and the overview), drawn at the same scale
	pans = NamedTuple[]
	for (col, (kind, L, N0)) in enumerate((("confined", 2, 600), ("free", 0, 600)))
		r = nothing
		for rows in (prow, trow)
			i = pick(rows, kind, L, N0; seed = 1)
			i === nothing || (r = rows[i]; break)
		end
		r === nothing && continue
		d = load(r); k = frame_at(d, 0.6)
		# one view each (Nov 3), the same oblique view for both (same view + same box = same scale);
		# the sphere whole with C translucent
		push!(pans, (panel = "ic_$(CORE)_$(kind)", row = 1, col = col, views = "oblique",
					 label = "", cells = cells(d, k), meta = meta(r, d, k), ghost = kind == "free"))
	end
	write_panels(joinpath(DATA, "fig1_initial_condition"), "", pans; same_scale = true)

	# ---- fig2: one aggregate size under no / medium / strong confinement, through time ---------
	# (your presentation sketch, panel B): N0 600, unconfined (polarises), 4 layers (central
	# cluster that shrinks), 1 layer (islands)
	pans = NamedTuple[]
	fBs = (0.6, 0.45, 0.3, 0.2, 0.12)
	rows_ = Any[]
	i = pick(trow, "free", 0, 600; seed = 1)
	i === nothing || push!(rows_, (trow[i], "No confinement", false))     # (titles: "No confinement, N_0 = 600")
	EX["CENTRAL"] === nothing || push!(rows_, (EX["CENTRAL"], "Medium confinement", true))
	EX["PATTERN"] === nothing || push!(rows_, (EX["PATTERN"], "Strong confinement", true))
	for (row, (r, name, confined)) in enumerate(rows_)
		d = load(r); kend = frame_at(d, 0.12)
		for (j, f) in enumerate(fBs)
			k = frame_at(d, f)
			c = confined ? cells(d, k) : align_B(cells(d, k); confined = false, ref = cells(d, kend))
			where_ = confined ? "$(lname(r["layers"])), \$N_0 = $(r["N0"])\$" : "\$N_0 = $(r["N0"])\$"
			lbl = (j == 1 ? "$name, $where_|" : "|") *
				  @sprintf("\$t = %.0f\$ h, \$f_B = %.2f\$%s", hrs(d, k), fB_of(d, k), f == READ_FB ? " (readout)" : "")
			# confined: top and side views; unconfined: one view, from the front (B to the side), whole
			push!(pans, (panel = "row$(row)_$(j)", row = row, col = j, views = confined ? "top+side" : "front",
						 label = lbl, cells = c, meta = meta(r, d, k), ghost = !confined))
		end
	end
	write_panels(joinpath(DATA, "fig2_behaviours"), "", pans)     # no title: the report's caption says it

	# ---- fig3: the patterns at the readout --------------------------------------------
	# every thickness that has patterns: its largest N0 first, then its smallest, up to 8
	pats = sort(filter(r -> pcls(r) == "PATTERN" && r["seed"] == 1 && major(r["layers"], r["N0"]) == "PATTERN", phase);
				by = r -> (r["layers"], r["N0"]))
	byL = [filter(r -> r["layers"] == L, pats) for L in sort(unique(r["layers"] for r in pats))]
	pats = unique(vcat([last(g) for g in byL], [first(g) for g in byL], pats))
	pats = sort(pats[1:min(end, 8)]; by = r -> (r["layers"], r["N0"]))
	pans = NamedTuple[]
	for (i, r) in enumerate(pats[1:min(end, 8)])
		d = load(r); k = read_frame(d)
		st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k])
		ex = exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; nmin = 10, frac = CORE_FRAC)
		push!(pans, (panel = @sprintf("L%g_N%d", r["layers"], r["N0"]), row = 1 + (i - 1) ÷ 4, col = 1 + (i - 1) % 4,
					 views = "top", label = @sprintf("%s, \$N_0 = %d\$", lname(r["layers"]), r["N0"]),
					 cells = cells(d, k), meta = meta(r, d, k)))
	end
	write_panels(joinpath(DATA, "fig3_patterns"), "", pans)

	# ---- fig4: unconfined aggregates of every size, at the readout ---------------------------
	# (the "free" row of the phase space): B polarises up to N0 ~ 600, bigger ones are
	# still central when B runs out
	pans = NamedTuple[]
	for (j, N0) in enumerate((150, 300, 600, 1200, 2400))
		i = pick(trow, "free", 0, N0; seed = 1)
		i === nothing && continue
		r = trow[i]; d = load(r); k = frame_at(d, READ_FB)
		ex = exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; nmin = 10, frac = CORE_FRAC)
		push!(pans, (panel = "free_N$(N0)", row = 1, col = j, views = "front",
					 label = @sprintf("\$N_0 = %d\$; \$P = %.2f\$|\$t = %.0f\$ h", N0, ex.P, hrs(d, k)),
					 cells = align_B(cells(d, k); confined = false), meta = meta(r, d, k), ghost = true))
	end
	write_panels(joinpath(DATA, "fig4_unconfined"), "", pans)

	# ---- fig5: the fusion assay -------------------------------------------------------
	# (Nov 3: results/stage5/fusion_nov3 has the starting doublets too; tune/ only the final ones)
	ff = first(filter(isfile, [joinpath(S5, "fusion_nov3", "fusion_$(BASE).jld2"), joinpath(S5, "tune", "fusion_$(BASE).jld2")]))
	if isfile(ff)
		pans = NamedTuple[]
		jldopen(ff, "r") do fh
			m = fh["metric"]
			for (j, pr) in enumerate(("BB", "CC", "BC"))
				me = (assay = "fusion_assay(n = 150, hours = 24)", pair = pr, t0 = m[pr][1], h12 = m[pr][2], h24 = m[pr][3], cand = CAND)
				# one column per pair: the two aggregates placed for fusion, an arrow, the doublet after 24 h
				haskey(fh, pr * "_0") && push!(pans, (panel = "start_$(pr)", row = 1, col = j, views = "front",
													 label = "\$$(pr)\$", cells = fh[pr * "_0"], meta = me))
				haskey(fh, pr * "_0") && push!(pans, (panel = "arrow_$(pr)", row = 2, col = j, views = "arrow",
													 label = "", cells = nothing, meta = me))
				push!(pans, (panel = "doublet_$(pr)", row = haskey(fh, pr * "_0") ? 3 : 1, col = j, views = "front",
							 label = haskey(fh, pr * "_0") ? "" : "\$$(pr)\$", cells = fh[pr], meta = me))
			end
		end
		write_panels(joinpath(DATA, "fig5_fusion"), "", pans)
	end

	# ---- fig6, fig7: a thick-slab ring through time, and B dissolving at mu BB 1.8 ------------
	# (Oct 8 runs; in the November report as ring-3layers.png and dissolves-muBB1.8.png)
	for (fname, file) in (("fig6_ring_3layers", joinpath(S5, PHASE, "com_L3_N2400_$(CAND)_rigid_w0_lw1.5_s1.jld2")),
						  ("fig7_dissolves_muBB1.8", joinpath(S5, "thick", "com_L2_N300_july_frB4.5_frC6_muB1.8@lw1.5_rigid_w0_lw1.5_s1.jld2")))
		isfile(file) || (println("missing ", file); continue)
		d = load_frames(file); pans = NamedTuple[]
		for (j, f) in enumerate((0.6, 0.4, 0.3, 0.2, 0.12))
			k = frame_at(d, f); c = cells(d, k)
			ex = exp_stats(c.x, c.y, c.z, c.cell_state; nmin = 10, frac = CORE_FRAC)
			push!(pans, (panel = "f$(j)", row = 1, col = j, views = "top+side",
						 label = @sprintf("\$f_B = %.2f\$, \$t = %.0f\$ h|%s, \$N_\\mathrm{eff} = %.1f\$", fB_of(d, k), hrs(d, k),
										  lowercase(classify_exp(ex)), ex.N_eff),
						 cells = c, meta = (run = basename(file), frame = k, t = d.t[k], hours = hrs(d, k), fB = fB_of(d, k))))
		end
		write_panels(joinpath(DATA, fname), "", pans)
	end

	# ---- phase space: at the readout fB = 0.2, and the same runs later (fB = 0.12) ------------
	ns = [count(r -> (r["layers"], r["N0"]) == p, phase) for p in unique((r["layers"], r["N0"]) for r in phase)]
	nseed_txt = minimum(ns) == maximum(ns) ? "$(minimum(ns)) seeds" : "$(minimum(ns))–$(maximum(ns)) seeds per point"
	# no titles on the charts (November): the readout is in readout.txt and the reports' captions
	# readout.txt: the class column, then the column that breaks ties between seeds (what wins later)
	for (folder, ccol, lcol, fread, ttl) in (("phase_space", CCOL, CCOL_LATE, nothing, ""),
											 ("phase_space_fB20", "class_exp50_02", "class_exp50_012", 0.2, ""),
											 ("phase_space_classic", "class_02", "class_012", 0.2, ""))
		dir = joinpath(DATA, folder); mkpath(dir)
		write(joinpath(dir, "readout.txt"), ccol * "\n" * lcol * "\n")
		write(joinpath(dir, "title.txt"), ttl)
		# A_um2: the pancake area at squeezing, pi (D r)^2 with r = 5 um (D in cell diameters)
		for r in phase; r["A_um2"] = round(pi * r["D"]^2 * 25; digits = 1); end
		keys_ = ["layers", "N0", "seed", "D", "A_um2", "class_04", "class_03", "class_02", "class_015", "class_012",
				 "class_exp_04", "class_exp_03", "class_exp_02", "class_exp_015", "class_exp_012", "class_exp10_02", "class_exp10_012",
				 "class_exp50_02", "class_exp50_015", "class_exp50_012", "class_exp50_010",
				 "class_read", "class_read_late", "f_cl50_02", "f_cl50_015", "f_cl50_012", "f_cl50_010", "t_02", "n_isl", "f_isl", "f_iso", "f_max", "S", "spread", "fill", "offset", "offset_max",
				 "n_cl", "N_eff", "f_cl", "P", "n_cl10", "N_eff10", "n_cl50", "N_eff50", "layers_eff", "N_end", "file"]
		open(joinpath(dir, "phase_space.csv"), "w") do io
			println(io, join(keys_, ","))
			for r in sort(phase; by = r -> (r["layers"], r["N0"], r["seed"]))
				println(io, join([string(r[k]) for k in keys_], ","))
			end
		end
		# the unconfined runs of the same mechanics (no plates): the "free" row of the chart
		fkeys = ["class_02", "class_012", "class_exp_02", "class_exp_012", "class_exp10_02", "class_exp10_012",
				  "class_exp50_02", "class_exp50_015", "class_exp50_012", "class_exp50_010", "class_read", "class_read_late"]
		open(joinpath(dir, "free.csv"), "w") do io
			println(io, "N0,seed,", join(fkeys, ","), ",offset,P,file")
			grid_N0 = Set(r["N0"] for r in phase)      # only the N0 values of the grid
			for r in sort(filter(r -> r["kind"] == "free" && String(r["cand"]) == BASE && String(r["core"]) == CORE &&
								 r["N0"] in grid_N0, trow);
						  by = r -> (r["N0"], r["seed"]))
				println(io, join(vcat(r["N0"], r["seed"], [r[k] for k in fkeys], @sprintf("%.4f", r["offset"]), @sprintf("%.4f", r["P"]), r["file"]), ","))
			end
		end
		# insets: one example per outcome at this readout, titled with the outcome's name
		# (the panel column is the class, which colours the title)
		open(joinpath(dir, "panels.tsv"), "w") do io
			println(io, "panel\tlayers\tN0\tviews\tlabel\tfile")
			for (key, cls) in (("POLARISED", "POLARISED"), ("CENTRAL_INSET", "CENTRAL"), ("PATTERN_INSET", "PATTERN"))
				r = get(EX, key, nothing); r === nothing && continue
				d = load(r); k = fread === nothing ? read_frame(d) : frame_at(d, fread)
				c = cls == "POLARISED" ? align_B(cells(d, k); confined = false) : cells(d, k)
				f = save_panel(dir, "inset_$(lowercase(cls))", c.x, c.y, c.z, c.r, c.cell_state; meta(r, d, k)...)
				println(io, join((lowercase(cls), r["layers"], r["N0"], "top+side", PHASE_NAMES[cls], f), '\t'))
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
		# unconfined runs whole with C translucent, confined ones with the mid-plane cutaway
		gh = f in ("unconfined.jld2", "polarised.jld2") ? GHOST : nothing
		record_pancake_frames(d, color_map, ks; fps = 6, t0 = T_RELAX, ghost = gh,
							  filename = joinpath(odir, replace(f, ".jld2" => ".mp4")))
		println("wrote movies/", replace(f, ".jld2" => ".mp4"))
	end
end

MODE in ("data", "all")   && write_data()
MODE in ("plot", "all")   && replot_figures(FIG)
MODE in ("movies", "all") && make_movies()
