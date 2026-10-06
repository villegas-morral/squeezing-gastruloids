# Verification runs that are not part of the phase map.
#
#   julia -t auto --project=code code/scripts/run_checks.jl [what]
#
# what = "polar"  unconfined polarisation with the final parameter set
#        "wall"   does wall_settings(2) = (60, 0.015) reproduce (150, 0.0075)?
#        "fusion" the two-aggregate doublet assay
#        "all"    (default)

using Random, Statistics, Printf, Distributions, ProgressMeter
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new

const WHAT   = length(ARGS) >= 1 ? ARGS[1] : "all"
const OUTDIR = joinpath(@__DIR__, "..", "results")
mkpath(OUTDIR); mkpath(joinpath(OUTDIR, "coms"))

if WHAT in ("polar", "all")
    println("\n=== unconfined polarisation, final parameter set ===")
    println("f_rep = ", FREP_SET)
    res = run_unconfined(300; tf = 700.0, outdir = joinpath(OUTDIR, "coms"))
    d = res.frames
    st = polar_stats(d.x[end], d.y[end], d.z[end], d.cell_state[end])
    @printf("FINAL  offset=%.3f  r_AB=%.3f  r_C=%.3f  fAB=%.3f\n",
            st.offset, st.r_AB, st.r_C, st.fAB)
    flush(stdout)
    # written so make_report.jl can quote the numbers instead of hardcoding them
    open(joinpath(OUTDIR, "polar.csv"), "w") do io
        println(io, "t,N,fAB,offset,r_AB,r_C")
        for k in res.m1:length(d.x)
            s = polar_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k])
            @printf(io, "%.1f,%d,%.4f,%.4f,%.4f,%.4f\n",
                    d.t[k], length(d.x[k]), s.fAB, s.offset, s.r_AB, s.r_C)
        end
    end
    fig = plot_frames(d, color_map,
                      round.(Int, range(res.m1, length(d.x); length = 6)))
    mkpath(joinpath(@__DIR__, "..", "figures"))
    save(joinpath(@__DIR__, "..", "figures", "polarisation.png"), fig)
    println("wrote results/polar.csv and figures/polarisation.png")
end

if WHAT in ("wall", "all")
    println("\n=== wall-setting equivalence at layers = 2 ===")
    println("The phase map uses wall_settings(layers). At layers >= 2 that is")
    println("(rep_wall 60, dt 0.015); the exploratory scan used (150, 0.0075).")
    println("Both must give the same verdict or the map is not comparable.\n")
    for (rw, dtd) in [(60.0, 0.015), (150.0, 0.0075)]
        # temporary override of wall_settings for this comparison
        @eval wall_settings(layers) = ($rw, $dtd)
        secs = @elapsed r = run_confined(2.0, 1700; outdir = joinpath(OUTDIR, "coms"),
                                          verbose = false)
        e = r.timeline[end]
        @printf("rep_wall %5.0f  dt %.4f -> %-9s n_isl=%d fmax=%.2f clust=%.2f  (%.0fs)\n",
                rw, dtd, r.verdict, e.n_isl, e.frac_max, e.clustered, secs)
        flush(stdout)
    end
    @eval wall_settings(layers) = layers <= 1 ? (150.0, 0.005) : (60.0, 0.015)   # restore the real one
end

if WHAT in ("fusion", "all")
    println("\n=== two-aggregate fusion assay, final parameter set ===")
    fusion_assay(MU_SET, FATT_SET, FREP_SET; n = 150, hours = 24.0)
    println("\nreference (Sotiris' set) at 24 h:")
    println("  AA 1.614   BB 0.537   CC 0.562   BC 0.947   AB 0.470   AC 0.657")
end

println("\ndone")
