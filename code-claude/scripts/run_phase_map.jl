# Pattern phase map: (layers, N0) -> islands or not.
#
#   julia -t auto --project=. code-claude/scripts/run_phase_map.jl [outdir] [seed]
#
# Appends one row per run to <outdir>/phase_map.csv, and skips any (layers, N0,
# seed) already present, so the scan can be stopped and resumed.
#
# Runtime: ~5-12 min per point on 20 threads; the full default grid is ~4 h.

using Random, Statistics, Printf, Distributions, ProgressMeter
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new

const OUTDIR = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "results")
const SEED   = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2345
const COMS   = joinpath(OUTDIR, "coms")
mkpath(OUTDIR); mkpath(COMS)
const CSV = joinpath(OUTDIR, "phase_map.csv")

# The grid. N0 is the cell count at the moment of squeezing. Chosen to bracket
# the transition at each thickness; the threshold grows like layers^3 because
# break-up is set by the aspect ratio of the A,B domain (see the report).
const GRID = [
    (1.0, [100, 200, 300, 500, 750]),
    (2.0, [340, 500, 750, 1150, 1700, 2500]),
    (3.0, [750, 1150, 1700, 2100, 2500]),
    (4.0, [1700, 2500, 3200, 4000]),
]

if !isfile(CSV)
    open(CSV, "w") do io
        println(io, "layers,N0,seed,verdict,n_isl,frac_max,clustered,fAB_end,zspread,pct_out,secs,file")
    end
end

done = Set{Tuple{Float64,Int,Int}}()
for line in Iterators.drop(eachline(CSV), 1)
    f = split(line, ',')
    length(f) >= 3 && push!(done, (parse(Float64, f[1]), parse(Int, f[2]), parse(Int, f[3])))
end
isempty(done) || println("resuming; $(length(done)) points already done")

for (L, N0v) in GRID, N0 in N0v
    (L, N0, SEED) in done && continue
    print(@sprintf("### layers %.0f  N0 %5d ... ", L, N0)); flush(stdout)
    local r, secs
    # One point must not kill a multi-hour scan. The realistic failure is the
    # Hookean plates going unstable (positions -> NaN, surfacing as
    # `InexactError: Int64(NaN)` from the neighbour grid); record it and move on.
    failed = false
    secs = @elapsed try
        r = run_confined(L, N0; seed = SEED, outdir = COMS, verbose = false)
    catch e
        failed = true
        @printf("FAILED: %s\n", first(split(sprint(showerror, e), '\n'))[1:min(end, 70)])
        flush(stdout)
    end
    if failed
        open(CSV, "a") do io
            @printf(io, "%.1f,%d,%d,FAILED,0,NaN,NaN,NaN,NaN,NaN,%.0f,-\n", L, N0, SEED, secs)
        end
        continue
    end
    d  = load_frames(r.file)
    e  = r.timeline[end]
    zz = d.z[end]; Lgap = (1.05*2.0)*L
    zsp  = maximum(zz) - minimum(zz)
    pout = 100*count(v -> v < 1.0 - 1e-9 || v > Lgap - 1.0 + 1e-9, zz)/length(zz)
    @printf("%-9s n_isl=%d fmax=%.2f clust=%.2f fAB=%.3f zsp=%.2f (%.0fs)\n",
            r.verdict, e.n_isl, e.frac_max, e.clustered, e.fAB, zsp, secs)
    flush(stdout)
    open(CSV, "a") do io
        @printf(io, "%.1f,%d,%d,%s,%d,%.4f,%.4f,%.4f,%.3f,%.1f,%.0f,%s\n",
                L, N0, SEED, r.verdict, e.n_isl, e.frac_max, e.clustered,
                e.fAB, zsp, pout, secs, basename(r.file))
    end
end

println("\ndone -> ", CSV)
