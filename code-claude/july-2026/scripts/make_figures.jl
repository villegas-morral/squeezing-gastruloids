# Figures for the report and for the experimentalists.
#
#   julia --project=code code/scripts/make_figures.jl [resultsdir]
#
# Reads results/phase_map.csv and results/coms/*.jld2 (written by
# run_phase_map.jl) and writes into figures/:
#
#   phase_map.png            layers vs N0, islands or not
#   pattern_L<n>.png         3D top-cutaway + side, ONE row per case,
#                            for each layer: below threshold and above
#   patterns_3D_all.png      the same cases side by side, one figure
#   timeline_L<n>_N<N0>.png  mid-plane slices through time
#
# No simulation is run here; everything comes from the saved trajectories.

using Statistics, Printf, DelimitedFiles
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_analysis.jl"))

const RES  = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "results")
const COMS = joinpath(RES, "coms")
const FIGS = joinpath(@__DIR__, "..", "figures")
mkpath(FIGS)

# --------------------------------------------------------------- read the map
csv = joinpath(RES, "phase_map.csv")
isfile(csv) || error("no $csv -- run scripts/run_phase_map.jl first")
raw, hdr = readdlm(csv, ','; header = true)
h = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
rows = [(layers = Float64(raw[r, h["layers"]]), N0 = Int(raw[r, h["N0"]]),
         verdict = strip(String(raw[r, h["verdict"]])),
         n_isl = Int(raw[r, h["n_isl"]]),
         frac_max = Float64(raw[r, h["frac_max"]]),
         clustered = Float64(raw[r, h["clustered"]]),
         file = strip(String(raw[r, h["file"]])))
        for r in 1:size(raw, 1)]

println("phase map: $(length(rows)) runs")
critical_N0(rows)

# ------------------------------------------------------------------ phase map
fig = plot_phase(rows)
save(joinpath(FIGS, "phase_map.png"), fig)
println("wrote figures/phase_map.png")

# ------------------------------------------- pick one no-island / one islands
"""
Per layer: the smallest N0 that gives islands, and the best no-island example.

For the no-island panel prefer a run classified SINGLE -- one central A,B
cluster that never breaks up, which is the case the experiment describes. A
`transient` run is excluded even though it also ends without islands, because it
*does* form them briefly and would misrepresent the sub-threshold behaviour.
"""
function pick_pair(rows, L)
    sub  = filter(r -> r.layers == L, rows)
    yes  = filter(r -> r.verdict == "ISLANDS", sub)
    ymin = isempty(yes) ? nothing : argmin(r -> r.N0, yes)
    below = ymin === nothing ? filter(r -> r.verdict != "ISLANDS", sub) :
                               filter(r -> r.N0 < ymin.N0, sub)
    singles = filter(r -> r.verdict == "SINGLE", below)
    pool = isempty(singles) ? filter(r -> r.verdict != "ISLANDS", below) : singles
    nmax = isempty(pool) ? nothing : argmax(r -> r.N0, pool)
    (without = nmax, with = ymin)
end

"""
Frame to show. The A,B phase is transient, so the LAST frame is well past the
clearest pattern -- pick the frame with the most islands instead, and for a
no-island run pick the frame at a comparable A,B fraction so the pair is a fair
visual comparison.
"""
function pick_frame(d, L; match_fAB = nothing, fAB_max = 0.40)
    fABs = [count(v -> v == 1 || v == 2, s)/length(s) for s in d.cell_state]

    # The saved trajectory starts with the GROWTH frames; initialize_diff! resets
    # t to 0, so the differentiation phase is everything from the last t == 0.
    # Scoring from frame 5 scored growth frames instead, and at layers 3-4 that
    # made argmax land on the freshly seeded core (t = 25, fAB = 0.6), which is
    # not a pattern.
    m1 = something(findlast(iszero, d.t), 1)

    # Restrict to the window where A,B is a clear minority but not yet gone.
    # Early frames are an unsorted A,B/C mixture whose "clusters" are seeding
    # noise, not islands.
    cand = [k for k in m1:length(d.x) if fABs[k] <= fAB_max]
    isempty(cand) && (cand = collect(m1:length(d.x)))

    if match_fAB !== nothing
        return cand[argmin(abs.(fABs[cand] .- match_fAB))], fABs
    end
    nisl = [island_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k], L).n_isl for k in cand]
    cand[argmax(nisl)], fABs
end

const LAYERS = sort(unique(r.layers for r in rows))
cases = NamedTuple[]

for L in LAYERS
    p = pick_pair(rows, L)
    # the islands case first, so the no-island case can be matched to its fAB
    fAB_ref = nothing
    for (tagname, r) in (("islands", p.with), ("no_islands", p.without))
        r === nothing && continue
        f = joinpath(COMS, r.file)
        isfile(f) || (@warn "missing $f"; continue)
        d = load_frames(f)
        k, fABs = tagname == "islands" ? pick_frame(d, L) :
                                         pick_frame(d, L; match_fAB = fAB_ref)
        tagname == "islands" && (fAB_ref = fABs[k])
        st = island_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k], L)
        push!(cases, (layers = L, N0 = r.N0, kind = tagname, n_isl = st.n_isl,
                      file = f, frame = k, t = d.t[k], fAB = fABs[k]))
        @printf("  layers %.0f  N0 %5d  %-10s -> frame %d (t=%.0f, fAB=%.2f, %d islands)\n",
                L, r.N0, tagname, k, d.t[k], fABs[k], st.n_isl)
    end
end
sort!(cases; by = c -> (c.layers, c.kind == "islands"))

# ---------------------------------------------- one 3D figure per layer, pair
for L in LAYERS
    sel = filter(c -> c.layers == L, cases)
    isempty(sel) && continue
    figs = []
    for c in sel
        d = load_frames(c.file)
        lbl = c.kind == "islands" ?
            @sprintf("layers %.0f,  N₀ = %d  →  %d islands   (t=%.0f, f_AB=%.2f)",
                     L, c.N0, c.n_isl, c.t, c.fAB) :
            @sprintf("layers %.0f,  N₀ = %d  →  no islands   (t=%.0f, f_AB=%.2f)",
                     L, c.N0, c.t, c.fAB)
        f = plot_pancake_frames(d, color_map, [c.frame]; n_views = 2, labels = [lbl],
                                cut_top = true, size_ratio = 1.3)
        push!(figs, (f, c))
    end
    # save them individually
    for (f, c) in figs
        nm = @sprintf("pattern_L%.0f_%s_N%d.png", c.layers, c.kind, c.N0)
        save(joinpath(FIGS, nm), f)
        println("wrote figures/", nm)
    end
end

# ------------------------------------- one combined figure, all layers at once
begin
    ncol = length(cases)
    fig = Figure(size = (620*ncol, 1120), figure_padding = 24)
    views = [(pi/2, false, true, "top"), (2*pi+pi/32, true, false, "side")]
    for (i, c) in enumerate(cases)
        d = load_frames(c.file)
        k = c.frame
        x, y, z, r, s = d.x[k], d.y[k], d.z[k], d.r[k], d.cell_state[k]
        zc = mean(z)
        bs = (maximum(x) - minimum(x))/1.9
        for (j, (persp, showz, showx, _)) in enumerate(views)
            keep = j == 1 ? findall(v -> v <= zc + 0.2, z) : collect(eachindex(z))
            ttl = j == 1 ? (c.kind == "islands" ?
                @sprintf("layers %.0f   N₀ = %d\n%d islands", c.layers, c.N0, c.n_isl) :
                @sprintf("layers %.0f   N₀ = %d\nno islands", c.layers, c.N0)) : ""
            ax = Axis3(fig[j, i], aspect = :data, azimuth = 0, elevation = persp,
                       xlabel = "", ylabel = "", zlabel = "",
                       xgridvisible = false, ygridvisible = false, zgridvisible = false,
                       xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
                       xticksvisible = false, yticksvisible = false, zticksvisible = false,
                       titlevisible = j == 1, titlesize = 26, titlegap = 8, title = ttl)
            o = sort(keep; by = m -> s[m], rev = true)
            meshscatter!(ax, x[o] .- mean(x), y[o] .- mean(y), z[o] .- zc;
                         markersize = r[o], color = [color_map[s[m]] for m in o])
            xlims!(ax, -bs, bs); ylims!(ax, -bs, bs); zlims!(ax, -bs, bs)
            ax.xticklabelsvisible = false
            ax.yticklabelsvisible = false
            ax.zticklabelsvisible = false
        end
    end
    Label(fig[0, 1:ncol],
          "Squeezed aggregates: A,B islands in a C sea - top cutaway (upper row), side view (lower row)";
          fontsize = 34, font = :bold, padding = (0, 0, 0, 30))
    rowgap!(fig.layout, 30)
    save(joinpath(FIGS, "patterns_3D_all.png"), fig)
    println("wrote figures/patterns_3D_all.png")
end

# ------------------------------------------------- time series for the islands
for c in filter(c -> c.kind == "islands", cases)
    d = load_frames(c.file)
    ks = round.(Int, range(5, length(d.x); length = 6))
    f = plot_frames(d, color_map, ks; slice = true, markersize = 4)
    nm = @sprintf("timeline_L%.0f_N%d.png", c.layers, c.N0)
    save(joinpath(FIGS, nm), f)
    println("wrote figures/", nm)
end

println("\ndone -> ", FIGS)
