# Stage 4: screen parameter candidates (reports/REPORT.md).
#
#   julia -t 12 --project=. code-claude/scripts/stage4_screen.jl <wall> <cand1,cand2,...>
#
# Only the BB and CC entries of f_att, f_rep and mu change; the BC entries stay at
# the FUSION values (frozen) and q, lambda, alpha_ecm, tau_div are untouched.
# For every candidate:
#
#   islands_L1   1 layer,  N0 1500   -> want sparse ISLANDS
#   islands_L2   2 layers, N0 2500   -> want sparse ISLANDS
#   small_L1     1 layer,  N0 200    -> want a CENTRAL cluster
#   polar        unconfined, N0 500  -> want the B centroid to reach the surface
#   fusion       B and C aggregates, 24 h -> want BC to glue (two lobes)
#
# Trajectories go to results/stage4/, one CSV row per run to results/stage4/screen.csv
# (readout at fB = 0.2 plus the end of the run). Resumable: done rows are skipped.

using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new

const WALL  = Symbol(length(ARGS) >= 1 ? ARGS[1] : "rigid")
const CANDS = split(length(ARGS) >= 2 ? ARGS[2] : "july", ",")
const OUT   = joinpath(@__DIR__, "..", "results", "stage4"); mkpath(OUT)
const CSV   = joinpath(OUT, "screen.csv")

include(joinpath(@__DIR__, "stage4_candidates.jl"))      # variant, CANDIDATES, parse_cand

isfile(CSV) || open(io -> println(io,
	"cand,wall,screen,layers,N0,ok,secs,t_read,fB_read,class_read,n_isl_read,f_isl_read,f_iso_read,f_max_read,offmax_read,S_read," *
	"t_end,fB_end,class_end,n_isl_end,S_end,offset_max_run,bb24,cc24,bc24,file"), CSV, "w")
done = Set{Tuple{String,String,String}}()
for line in Iterators.drop(eachline(CSV), 1)
	f = split(line, ','); length(f) >= 3 && push!(done, (f[1], f[2], f[3]))
end

row(io, cand, scr, L, N0, ok, secs, rr, re, offrun, fus, file) = @printf(io,
	"%s,%s,%s,%g,%d,%s,%.0f,%.1f,%.4f,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.1f,%.4f,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%s\n",
	cand, WALL, scr, L, N0, ok, secs, rr.t, rr.fB, rr.class, rr.n_isl, rr.f_isl, rr.f_iso, rr.f_max,
	rr.offset_max, rr.S, re.t, re.fB, re.class, re.n_isl, re.S, offrun, fus..., file)

const NOROW = (t = NaN, fB = NaN, class = "-", n_isl = 0, f_isl = NaN, f_iso = NaN, f_max = NaN, offset_max = NaN, S = NaN)

for cand in CANDS
	# "name@w0.1" = candidate `name` with a diluted B core (initialize_diff! width 0.1)
	base, width, lw = parse_cand(cand)
	P = CANDIDATES[base]()
	PSETS[cand] = P
	println("\n######## candidate $cand  (wall $WALL, width $width)")
	println("  f_att BB $(P.f_att[2,2])  CC $(P.f_att[3,3])   f_rep BB $(P.f_rep[2,2])  CC $(P.f_rep[3,3])   mu BB $(P.mu[2,2])  CC $(P.mu[3,3])")
	flush(stdout)

	for (scr, L, N0) in [("small_L1", 1, 200), ("islands_L1", 1, 1500), ("islands_L2", 2, 2500)]
		(cand, string(WALL), scr) in done && continue
		ok = true; local r
		secs = @elapsed try
			r = run_confined(L, N0; pset = cand, wall = WALL, width = width, lambda_wall = lw, tf = 625.0, fB_stop = 0.12,
							 outdir = OUT, verbose = false)
		catch e
			ok = false; println("  $scr FAILED: ", first(split(sprint(showerror, e), '\n')))
		end
		if !ok
			open(io -> row(io, cand, scr, L, N0, false, secs, NOROW, NOROW, NaN, (NaN, NaN, NaN), "-"), CSV, "a")
			continue
		end
		d  = load_frames(r.file)
		tl = pattern_timeline(d; from = r.m1, confined = true)
		rr = readout_row(tl; fB_read = 0.2); re = tl[end]
		offrun = maximum(x -> isnan(x.offset_max) ? 0.0 : x.offset_max, tl)
		@printf("  %-11s t_read %4.0f  %-9s n_isl %2d  f_iso %.2f  f_max %.2f  S %.2f   | end %-9s (%.0fs)\n",
				scr, rr.t, rr.class, rr.n_isl, rr.f_iso, rr.f_max, rr.S, re.class, r.secs)
		flush(stdout)
		open(io -> row(io, cand, scr, L, N0, true, r.secs, rr, re, offrun, (NaN, NaN, NaN), basename(r.file)), CSV, "a")
	end

	if !((cand, string(WALL), "polar") in done)
		ok = true; local u
		secs = @elapsed try
			u = run_unconfined(500; pset = cand, width = width, tf = 625.0, fB_stop = 0.08, outdir = OUT)
		catch e
			ok = false; println("  polar FAILED: ", first(split(sprint(showerror, e), '\n')))
		end
		if ok
			tl = pattern_timeline(u.frames; from = u.m1, confined = false)
			rr = readout_row(tl; fB_read = 0.2); re = tl[end]
			offrun = maximum(x -> isnan(x.offset_max) ? 0.0 : x.offset_max, tl)
			@printf("  %-11s offset_max at fB 0.2: %.2f   max over run: %.2f   (%s at readout)\n",
					"polar", rr.offset_max, offrun, rr.class)
			open(io -> row(io, cand, "polar", 0, 500, true, u.secs, rr, re, offrun, (NaN, NaN, NaN), basename(u.file)), CSV, "a")
		else
			open(io -> row(io, cand, "polar", 0, 500, false, secs, NOROW, NOROW, NaN, (NaN, NaN, NaN), "-"), CSV, "a")
		end
		flush(stdout)
	end

	# the fusion assay does not depend on the core: run it once per base candidate
	if !((base, string(WALL), "fusion") in done)
		push!(done, (base, string(WALL), "fusion"))
		pos = Dict()
		secs = @elapsed res = fusion_assay(P.mu, P.f_att, P.f_rep; n = 150, hours = 24.0,
										   pairs_ = [(2,2),(3,3),(2,3)], verbose = false,
										   positions = pos)
		fus = (res["BB"][3], res["CC"][3], res["BC"][3])
		@printf("  %-11s 24 h:  BB %.3f   CC %.3f   BC %.3f   (lower = more fused)\n", "fusion", fus...)
		ffile = joinpath(OUT, "fusion_$(base).jld2")
		jldsave(ffile; BB = pos["BB"], CC = pos["CC"], BC = pos["BC"], metric = res)
		open(io -> row(io, base, "fusion", 0, 150, true, secs, NOROW, NOROW, NaN, fus, basename(ffile)), CSV, "a")
		flush(stdout)
	end
end
println("\ndone -> ", CSV)
