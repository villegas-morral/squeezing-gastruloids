# Stage 2: compare the plate models (reports/REPORT.md).
#
#   julia -t 12 --project=. code-claude/scripts/stage2_walls.jl constant,rigid
#   julia -t 12 --project=. code-claude/scripts/stage2_walls.jl hookean,hookean_stiff
#   julia -t 12 --project=. code-claude/scripts/stage2_walls.jl rigid 2     # only the 2-layer case
#
# Same aggregate (July set, sharp B core, seed 2345) between plates of 1, 2 and 4
# layers, ~20 h of model time. For each wall model: how thick the slab stays
# (thickness_timeline), whether it is stable, what it costs, and whether the
# pattern outcome changes. Appends to results/stage2/walls.csv.

using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new

const WALLS = Symbol.(split(length(ARGS) >= 1 ? ARGS[1] : "constant,hookean,hookean_stiff,rigid", ","))
const CASES = filter(c -> length(ARGS) < 2 || string(c[1]) in split(ARGS[2], ","),
                     [(1, 600), (2, 1500), (4, 1500)])      # (layers, N0); optional 2nd arg: layers to run
const OUT   = joinpath(@__DIR__, "..", "results", "stage2"); mkpath(OUT)
const CSV   = joinpath(OUT, "walls.csv")
isfile(CSV) || open(io -> println(io, "wall,layers,N0,ok,secs,spread98_max,spread98_end,layers_eff_end,outside_end,fB_end,verdict,n_isl,frac_max,file"), CSV, "w")

for wall in WALLS, (layers, N0) in CASES
    print(@sprintf("### %-14s layers %d  N0 %5d ... ", wall, layers, N0)); flush(stdout)
    local r
    ok = true
    secs = @elapsed try
        r = run_confined(layers, N0; wall = wall, tf = 300.0, outdir = OUT, verbose = false)
    catch e
        ok = false
        println("FAILED: ", first(split(sprint(showerror, e), '\n')))
    end
    if !ok
        open(io -> @printf(io, "%s,%d,%d,false,%.0f,NaN,NaN,NaN,NaN,NaN,FAILED,0,NaN,-\n", wall, layers, N0, secs), CSV, "a")
        continue
    end
    d  = load_frames(r.file)
    L  = (1.05*2.0)*layers
    th = thickness_timeline(d, L; from = r.m1)
    e  = r.timeline[end]
    @printf("ok (%.0fs)  spread98 max %.2f end %.2f  layers_eff %.2f  outside %.0f%%  -> %s n_isl=%d\n",
            r.secs, maximum(x.spread98 for x in th), th[end].spread98, th[end].layers_eff,
            100*th[end].outside, r.verdict, e.n_isl)
    flush(stdout)
    open(CSV, "a") do io
        @printf(io, "%s,%d,%d,true,%.0f,%.4f,%.4f,%.4f,%.4f,%.4f,%s,%d,%.4f,%s\n",
                wall, layers, N0, r.secs, maximum(x.spread98 for x in th), th[end].spread98,
                th[end].layers_eff, th[end].outside, e.fAB, r.verdict, e.n_isl, e.frac_max,
                basename(r.file))
    end
end
println("done -> ", CSV)
