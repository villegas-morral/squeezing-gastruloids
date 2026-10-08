# Stage 5b: the phase space (reports/REPORT.md).
#
#   julia --project=. code-claude/scripts/stage5_phase.jl jobs  <cand> <core> <out>   -> results/stage5/<out>.csv
#   julia --project=. code-claude/scripts/stage5_phase.jl reuse <cand> <core> <out>   -> copy the matching tuning runs in
#   julia --project=. code-claude/scripts/stage5_phase.jl table <cand> <core> <out>   -> summary + map + refine jobs
#
# jobs:  the rough grid, layers {1, 2, 3, 4, 6, 8, 12} x N0 {75, 150, 300, 600, 1200,
#        2400} plus N0 4800 at 1-4 layers, 3 seeds each; run them with stage5_run.jl.
# table: reads results/stage5/<out>/rows_w*.csv and writes there
#          phase_table.md   the majority class (and every seed) per point, at the readout fB = 0.2
#          phase_map.png    the chart (plot_phase_space), for checking -- the final one
#                           is drawn by stage6_figures.jl from figures/data
#          refine.csv       new jobs where the majority class changes between two
#                           neighbouring points: the geometric-mean N0 between them, and the
#                           in-between number of layers (5, 10) where there is one
#        Run refine.csv into the same <out> folder, then `table` again.

using Statistics, Printf, DelimitedFiles, JLD2, CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))      # read_rows, with_readout!

const MODE = ARGS[1]
const CAND = ARGS[2]
const CORE = ARGS[3]
const TAG  = ARGS[4]
const ROOT = joinpath(@__DIR__, "..", "results", "stage5")
const SEEDS = 1:3

if MODE == "jobs"
	open(joinpath(ROOT, "$(TAG).csv"), "w") do io
		println(io, "# stage 5b phase space of $CAND ($CORE core): kind,cand,core,layers,N0,seed")
		pts = [(L, N0) for L in (1, 2, 3, 4, 6, 8, 12) for N0 in (75, 150, 300, 600, 1200, 2400)]
		append!(pts, [(L, 4800) for L in (1, 2, 3, 4)])
		sort!(pts; by = p -> -p[1] * p[2])           # the longest runs first, so none is left for the end
		for s in SEEDS, (L, N0) in pts
			println(io, "confined,$CAND,$CORE,$L,$N0,$s")
		end
	end
	println("wrote ", joinpath(ROOT, "$(TAG).csv"))
end
if MODE in ("jobs", "reuse")
	# the tuning runs of the same set at grid points are the same simulations (same
	# seeds): copy them, so stage5_run.jl skips those jobs ("reuse" alone: copy only,
	# e.g. once more tuning runs have finished; the job list is left as it is)
	dst = joinpath(ROOT, TAG); mkpath(dst)
	grid = Set((Float64(L), N0) for L in (1, 2, 3, 4, 6, 8, 12) for N0 in (75, 150, 300, 600, 1200, 2400, 4800))
	tdir = joinpath(ROOT, "tune")
	reuse = filter(r -> r["ok"] == true && r["kind"] == "confined" && r["cand"] == CAND && r["core"] == CORE &&
						(Float64(r["layers"]), Int(r["N0"])) in grid, read_rows(tdir))
	hdr = readline(joinpath(tdir, first(filter(f -> startswith(f, "rows_w"), readdir(tdir)))))
	open(joinpath(dst, "rows_w0.csv"), "w") do io
		println(io, hdr)
		for r in reuse
			cp(joinpath(tdir, r["file"]), joinpath(dst, r["file"]); force = true)
			println(io, join([string(r[k]) for k in split(hdr, ',')], ","))
		end
	end
	println("reused $(length(reuse)) tuning runs -> ", joinpath(dst, "rows_w0.csv"))
	exit()
end

# ---- table ----------------------------------------------------------------------------
dir = joinpath(ROOT, TAG)
sel = filter(r -> r["ok"] == true && r["kind"] == "confined" && r["cand"] == CAND && r["core"] == CORE, read_rows(dir))
with_readout!(sel)                       # classes recomputed with the current classifier
rows = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), seed = Int(r["seed"]), class = r["class_02"],
		 class_03 = r["class_03"], class_012 = r["class_012"], S = Float64(r["S"]), n_isl = Int(r["n_isl"]),
		 offset = Float64(r["offset"]), spread = Float64(r["spread"]), fill = Float64(r["fill"]),
		 layers_eff = Float64(r["layers_eff"]), secs = Float64(r["secs"]), file = r["file"]) for r in sel]
println("$(length(rows)) runs of $CAND ($CORE)")

order = ["CENTRAL", "PATTERN", "POLARISED", "OTHER"]
Ls  = sort(unique(r.layers for r in rows))
N0s = sort(unique(r.N0 for r in rows))
at(L, N0) = filter(r -> r.layers == L && r.N0 == N0, rows)
function major(rs)
	isempty(rs) && return ""
	cnt = [count(r -> phase_class(r.class) == c, rs) for c in order]
	order[argmax(cnt)]
end
short(c) = Dict("CENTRAL" => "C", "POLARISED" => "P", "RING" => "R", "ISLANDS" => "I", "OTHER" => "o", "NONE" => "-")[c]

io = IOBuffer()
println(io, "# Phase space of $CAND ($CORE core), class at fB = 0.2\n")
println(io, "Majority over seeds (C central, PAT pattern = ring or islands, P polarised, o other); every seed in brackets",
		" (R ring, I islands). Rows: plate gap in cell layers; columns: initial cells N0.\n")
println(io, "| layers \\ N0 | ", join(N0s, " | "), " |")
println(io, "|---|", repeat("---|", length(N0s)))
for L in reverse(Ls)
	cells = String[]
	for N0 in N0s
		rs = sort(at(L, N0); by = r -> r.seed)
		m = major(rs)
		push!(cells, isempty(rs) ? "" : @sprintf("**%s** (%s)", Dict("CENTRAL" => "C", "PATTERN" => "PAT", "POLARISED" => "P",
							  "OTHER" => "o")[m], join(short.(getfield.(rs, :class)), " ")))
	end
	println(io, "| **", @sprintf("%g", L), "** | ", join(cells, " | "), " |")
end
# thickness: only slabs the aggregate fills (<= 3 layers); in a wider gap a small
# aggregate does not reach both plates, so its z-extent is below the gap anyway
thin = filter(r -> r.layers <= 3, rows)
println(io, @sprintf("\nThickness at the end, slabs of 1-3 layers: layers_eff / nominal = %.3f ± %.3f over %d runs. Run time %.0f s median, %.0f s max.",
		mean(r.layers_eff / r.layers for r in thin), std([r.layers_eff / r.layers for r in thin]), length(thin),
		median([r.secs for r in rows]), maximum(r.secs for r in rows)))
# how stable the class is around the readout
nchg = count(r -> phase_class(r.class_03) != phase_class(r.class) || phase_class(r.class_012) != phase_class(r.class), rows)
println(io, "Runs whose class at fB = 0.3 or 0.12 differs from the one at 0.2: $nchg of $(length(rows)).")
txt = String(take!(io)); print(txt)
write(joinpath(dir, "phase_table_$(CAND)_$(CORE).md"), txt)

# the unconfined runs of the same mechanics (no plates, so the base set) as the top row:
# from tuning, and from this folder when it has some
base = String(split(CAND, "@")[1])
free = [(N0 = Int(r["N0"]), class = r["class_02"]) for r in
		with_readout!(filter(r -> r["ok"] == true && r["kind"] == "free" && r["cand"] == base && r["core"] == CORE &&
							 Int(r["N0"]) in N0s,
							 vcat(read_rows(joinpath(ROOT, "tune")), read_rows(dir))))]
fig = plot_phase_space(rows; free = free, title = "$CAND ($CORE core) at fB = 0.2")
save(joinpath(dir, "phase_map_$(CAND)_$(CORE).png"), fig)

# ---- refinement jobs: between neighbours whose majority differs --------------------------
ref = Set{Tuple{Float64,Int}}()
for L in Ls
	Ns = filter(N0 -> !isempty(at(L, N0)), N0s)
	for i in 1:length(Ns)-1
		major(at(L, Ns[i])) != major(at(L, Ns[i+1])) || continue
		Nm = round(Int, sqrt(Ns[i] * Ns[i+1]))
		(Ns[i+1] / Ns[i] > 1.3) && push!(ref, (L, Nm))
	end
end
for N0 in N0s
	Lh = filter(L -> !isempty(at(L, N0)), Ls)
	for i in 1:length(Lh)-1
		major(at(Lh[i], N0)) != major(at(Lh[i+1], N0)) || continue
		Lm = round((Lh[i] + Lh[i+1]) / 2)
		(Lm != Lh[i] && Lm != Lh[i+1]) && push!(ref, (Lm, N0))
	end
end
ref = sort([p for p in ref if isempty(at(p...))]; by = p -> -p[1] * p[2])
open(joinpath(ROOT, "$(TAG)_refine.csv"), "w") do io
	println(io, "# refinement of $TAG: points between neighbours whose majority class differs")
	for s in SEEDS, (L, N0) in ref
		println(io, @sprintf("confined,%s,%s,%g,%d,%d", CAND, CORE, L, N0, s))
	end
end
println("$(length(ref)) refinement points -> ", joinpath(ROOT, "$(TAG)_refine.csv"))
