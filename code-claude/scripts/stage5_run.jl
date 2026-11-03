# Stage 5: run a list of simulations, split over several Julia processes
# (reports/REPORT.md).
#
#   julia --project=. code-claude/scripts/stage5_run.jl <jobs.csv> <worker> <n_workers> [out]
#
# jobs.csv has one job per line, `kind,cand,core,layers,N0,seed`:
#
#   kind    confined | free (unconfined, the polarisation check) | fusion (doublet assay)
#   cand    a name from stage4_candidates.jl, optionally with "@lwX": friction X with
#           the plates (run_confined lambda_wall; unconfined runs have no plates)
#           "@hookean" uses the July spring plates instead of the rigid ones
#           "@lamX" multiplies the friction lambda by X (Oct 8 test; lambda is meant to stay)
#           "@tdivX" mean division time X (time units) for every cell type (Nov test:
#           @tdiv195 = 13 h for B and C, instead of 300 and 750)
#   core    sharp | fin0.75 (75 % of the B cells inside the core, initialize_diff! f_in)
#           | fB0.65 (65 % B at seeding instead of 60 %)
#           | cc0.1 (10 % of the core left C, mixed in: the July seeding)
#           | several joined by + (fB0.65+cc0.1)
#
# Two ways to share the work between processes (each uses one core):
#
#   n_workers > 0   worker k of n takes jobs k, k+n, k+2n, ...
#   n_workers = 0   a queue: each worker claims the next job nobody has claimed yet
#                   (an empty folder claims/<job>; creating a folder is atomic), so
#                   workers can be added at any time
#
# Trajectories go to results/stage5/<out>/ (default <out> = the jobs file's name),
# one CSV row per run to rows_w<k>.csv there (one file per worker, so they never
# write to the same file). Resumable: a job with a row in any rows_*.csv is skipped;
# in queue mode a job that crashed keeps its claim -- delete claims/ to retry it.
#
# Every run goes on until the B fraction drops to 0.10. The row holds the class at
# fB = 0.4, 0.3, 0.2, 0.15 and 0.12, the full pattern stats at the readout fB = 0.2,
# and the effective number of layers at the end (the plates must hold it).

using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new
include(joinpath(@__DIR__, "stage4_candidates.jl"))      # variant, CANDIDATES

const JOBS   = ARGS[1]
const WORKER = parse(Int, ARGS[2])
const NWORK  = parse(Int, ARGS[3])
const OUT    = joinpath(@__DIR__, "..", "results", "stage5",
					  length(ARGS) >= 4 ? ARGS[4] : splitext(basename(JOBS))[1]); mkpath(joinpath(OUT, "claims"))
const CSV    = joinpath(OUT, "rows_w$(WORKER).csv")
const FBS    = (0.4, 0.3, 0.2, 0.15, 0.12)

const HEADER = "kind,cand,core,layers,N0,seed,ok,secs,N_end,layers_eff," *
			   "class_04,class_03,class_02,class_015,class_012," *
			   "t_02,n_isl,f_isl,f_iso,f_max,S,spread,offset,offset_max,offset_run,bb24,cc24,bc24,file"
isfile(CSV) || open(io -> println(io, HEADER), CSV, "w")

key(f) = join(f[1:6], ",")
done = Set{String}()
for f in readdir(OUT)
	(startswith(f, "rows_w") && endswith(f, ".csv")) || continue
	for line in Iterators.drop(eachline(joinpath(OUT, f)), 1)
		p = split(line, ','); length(p) >= 6 && push!(done, key(p))
	end
end

jobs = [split(strip(l), ',') for l in eachline(JOBS) if !isempty(strip(l)) && !startswith(strip(l), "#")]
# a job already done (e.g. a run copied in from another folder) counts as claimed, so
# stage5_pool.sh does not keep starting workers for it
for k in done; mkpath(joinpath(OUT, "claims", replace(k, "," => "_"))); end
mine = [j for (i, j) in enumerate(jobs) if (NWORK == 0 || (i - 1) % NWORK == WORKER - 1) && !(key(j) in done)]
println("worker $WORKER/$(NWORK == 0 ? "queue" : NWORK): $(length(mine)) jobs left"); flush(stdout)

"Queue mode: true if this worker got the job (nobody had claimed it)."
claim(j) = NWORK > 0 || (try; mkdir(joinpath(OUT, "claims", join(j[1:6], "_"))); true; catch; false; end)

"Frame where the B fraction first drops to f, after the seeding frame m1."
function frame_at(d, m1, f)
	k = findfirst(i -> i >= m1 && count(==(2), d.cell_state[i]) / length(d.cell_state[i]) <= f, eachindex(d.x))
	k === nothing ? length(d.x) : k
end

"Core spec -> (f_in, fB, c_core): sharp | fin0.75 | fB0.65 | cc0.1 | combined with + (fB0.65+cc0.1)."
function core_of(core)
	f_in, fB, cc = nothing, 0.6, 0.0
	for p in split(core, '+')
		p == "sharp" ? nothing :
		startswith(p, "fin") ? (f_in = parse(Float64, p[4:end])) :
		startswith(p, "fB")  ? (fB = parse(Float64, p[3:end])) :
		startswith(p, "cc")  ? (cc = parse(Float64, p[3:end])) : error("unknown core $p")
	end
	(f_in, fB, cc)
end

for j in mine
	claim(j) || continue
	kind, cand, core = String(j[1]), String(j[2]), String(j[3])
	L, N0, seed = parse(Float64, j[4]), parse(Int, j[5]), parse(Int, j[6])
	wall = occursin("@hookean", cand) ? :hookean : :rigid  # "x@hookean": the July spring plates
	m = match(r"@lam([0-9.]+)", cand)                      # "x@lam2": friction lambda x 2
	lam = m === nothing ? 1.0 : parse(Float64, m[1])
	m = match(r"@tdiv([0-9.]+)", cand)                     # "x@tdiv195": every cell divides every 195 (13 h)
	tdiv = m === nothing ? nothing : fill(parse(Float64, m[1]), 3)
	base, _, lw = parse_cand(replace(cand, "@hookean" => "", r"@lam[0-9.]+" => "", r"@tdiv[0-9.]+" => ""))   # "fusion@lw3" -> ("fusion", 0, 3)
	haskey(PSETS, cand) || (PSETS[cand] = CANDIDATES[base]())
	P = PSETS[cand]
	@printf("### %-8s %-18s %-8s L %-3g N0 %5d seed %d ... ", kind, cand, core, L, N0, seed); flush(stdout)
	t0 = time()
	try
		if kind == "fusion"
			pos = Dict()
			res = fusion_assay(P.mu, P.f_att, P.f_rep; n = N0, hours = 24.0,
							   pairs_ = [(2,2),(3,3),(2,3)], verbose = false, positions = pos)
			file = joinpath(OUT, "fusion_$(cand).jld2")
			jldsave(file; BB = pos["BB"], CC = pos["CC"], BC = pos["BC"], metric = res)
			b = (res["BB"][3], res["CC"][3], res["BC"][3])
			@printf("BB %.3f CC %.3f BC %.3f\n", b...)
			open(io -> @printf(io, "%s,%s,%s,%g,%d,%d,true,%.0f,%d,NaN,-,-,-,-,-,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,%.4f,%.4f,%.4f,%s\n",
							   kind, cand, core, L, N0, seed, time() - t0, 2N0, b..., basename(file)), CSV, "a")
			continue
		end
		confined = kind == "confined"
		f_in, fB0, cc = core_of(core)
		r = confined ?
			run_confined(L, N0; pset = cand, wall = wall, f_in = f_in, fB = fB0, c_core = cc, lambda_scale = lam, lambda_wall = lw,
						 tau_div = tdiv, tf = 900.0, fB_stop = 0.10, seed = seed, outdir = OUT, verbose = false) :
			run_unconfined(N0; pset = cand, f_in = f_in, fB = fB0, c_core = cc, lambda_scale = lam, tau_div = tdiv, tf = 900.0,
						   fB_stop = 0.10, seed = seed, outdir = OUT, verbose = false)
		d = load_frames(r.file)
		st_at(k) = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = confined)
		cls = [classify_pattern(st_at(frame_at(d, r.m1, f))) for f in FBS]
		k2 = frame_at(d, r.m1, 0.2); s = st_at(k2)
		off_run = maximum(st_at(k).offset for k in r.m1:length(d.x))
		Leff = confined ? thickness_timeline(d, 2.1L; from = length(d.z))[end].layers_eff : NaN
		@printf("%-9s %-9s %-9s %-9s %-9s (%.0fs)\n", cls..., time() - t0)
		open(io -> @printf(io, "%s,%s,%s,%g,%d,%d,true,%.0f,%d,%.3f,%s,%s,%s,%s,%s,%.1f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,NaN,NaN,NaN,%s\n",
						   kind, cand, core, L, N0, seed, time() - t0, length(d.x[end]), Leff, cls...,
						   d.t[k2], s.n_isl, s.f_isl, s.f_iso, s.f_max, s.S, s.spread, s.offset, s.offset_max,
						   off_run, basename(r.file)), CSV, "a")
	catch e
		println("FAILED: ", first(split(sprint(showerror, e), '\n')))
		open(io -> @printf(io, "%s,%s,%s,%g,%d,%d,false,%.0f,0,NaN,-,-,-,-,-,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,-\n",
						   kind, cand, core, L, N0, seed, time() - t0), CSV, "a")
	end
	flush(stdout)
end
println("worker $WORKER done")
