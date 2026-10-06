# Find a timestep that survives the largest monolayer in the grid.
#
#   julia -t auto --project=. code-claude/scripts/probe_stability.jl [N0] [dt,dt,...]
#
# The Hookean plates blow up (positions -> NaN, surfacing as
# `InexactError: Int64(NaN)` from the neighbour grid) well below the naive
# linear-stability bound rep_wall*dt < 2, because the spring switches on and off
# discontinuously at the contact point. This probes empirically.

using Random, Statistics, Printf, Distributions, ProgressMeter
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new

const N0  = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 750
const DTS = length(ARGS) >= 2 ? parse.(Float64, split(ARGS[2], ",")) : [0.0075, 0.005, 0.004]
const RW  = 150.0

@printf("monolayer N0 = %d, rep_wall = %.0f\n", N0, RW)
@printf("%8s %10s  %-10s %s\n", "dt", "rw*dt", "result", "z-spread")
for dtd in DTS
    @eval wall_settings(layers) = ($RW, $dtd)
    ok = true; msg = ""; zsp = NaN
    secs = @elapsed try
        r = run_confined(1.0, N0; outdir = joinpath(@__DIR__, "..", "results", "probe"),
                         verbose = false)
        d = load_frames(r.file)
        zsp = maximum(d.z[end]) - minimum(d.z[end])
        msg = r.verdict
    catch e
        ok = false
        msg = first(split(sprint(showerror, e), '\n'))[1:min(end, 40)]
    end
    @printf("%8.4f %10.3f  %-10s %s   (%.0fs)\n", dtd, RW*dtd,
            ok ? "OK" : "BLEW UP", isnan(zsp) ? "-" : @sprintf("%.3f", zsp), secs)
    flush(stdout)
    ok && break          # the first stable one is the cheapest stable one
end
