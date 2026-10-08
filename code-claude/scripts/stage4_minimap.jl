# Stage 4: a coarse phase map per finalist, drawn as a 3D gallery laid out like the
# phase space (reports/REPORT.md).
#
#   julia -t 12 --project=. code-claude/scripts/stage4_minimap.jl <wall> <cand1,cand2,...> [run|plot]
#
# For every candidate: layers {1, 2, 4, 8} x N0 {150, 600, 2400}, one seed, run until
# the B fraction drops to 0.12. One CSV row per run (results/stage4/minimap.csv,
# resumable), with the class at fB = 0.3, 0.2 and 0.12. "plot" (or the end of
# "run") draws results/stage4/minimap_<cand>.png: rows = layers (thick on top),
# columns = N0, top cutaways at the readout fB = 0.2.

using Random, Statistics, Printf, Distributions, ProgressMeter, DelimitedFiles
using CellBasedModels
import CellBasedModels: update!
using JLD2, CairoMakie

include(joinpath(@__DIR__, "gastruloids_lib.jl"))
model = model_new
include(joinpath(@__DIR__, "stage4_candidates.jl"))

const WALL   = Symbol(length(ARGS) >= 1 ? ARGS[1] : "rigid")
const CANDS  = split(length(ARGS) >= 2 ? ARGS[2] : "fusion", ",")
const MODE   = length(ARGS) >= 3 ? ARGS[3] : "run"
const LAYERS = [1, 2, 4, 8]
const N0S    = [150, 600, 2400]
const OUT    = joinpath(@__DIR__, "..", "results", "stage4"); mkpath(OUT)
const CSV    = joinpath(OUT, "minimap.csv")

isfile(CSV) || open(io -> println(io, "cand,wall,layers,N0,ok,secs,class_03,class_02,class_012,n_isl_02,S_02,spread_02,offset_02,file"), CSV, "w")
done = Set{Tuple{String,Int,Int}}()
for line in Iterators.drop(eachline(CSV), 1)
	f = split(line, ','); length(f) >= 4 && push!(done, (f[1], parse(Int, f[3]), parse(Int, f[4])))
end

"Frame where the B fraction first drops to f, after the seeding frame m1."
function frame_at(d, m1, f)
	k = findfirst(i -> i >= m1 && count(==(2), d.cell_state[i]) / length(d.cell_state[i]) <= f, eachindex(d.x))
	k === nothing ? length(d.x) : k
end
st_at(d, k) = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = true)

if MODE == "run"
	for cand in CANDS
		base, width, lw = parse_cand(cand)
		PSETS[cand] = CANDIDATES[base]()
		for L in LAYERS, N0 in N0S
			(cand, L, N0) in done && continue
			print(@sprintf("### %-22s layers %d  N0 %5d ... ", cand, L, N0)); flush(stdout)
			ok = true; local r
			secs = @elapsed try
				r = run_confined(L, N0; pset = cand, wall = WALL, width = width, lambda_wall = lw, tf = 625.0,
								 fB_stop = 0.12, outdir = OUT, verbose = false)
			catch e
				ok = false; println("FAILED: ", first(split(sprint(showerror, e), '\n')))
			end
			if !ok
				open(io -> @printf(io, "%s,%s,%d,%d,false,%.0f,-,-,-,0,NaN,NaN,NaN,-\n", cand, WALL, L, N0, secs), CSV, "a")
				continue
			end
			d = load_frames(r.file)
			cls = [classify_pattern(st_at(d, frame_at(d, r.m1, f))) for f in (0.3, 0.2, 0.12)]
			s2 = st_at(d, frame_at(d, r.m1, 0.2))
			@printf("%-9s %-9s %-9s  (%.0fs)\n", cls..., r.secs); flush(stdout)
			open(io -> @printf(io, "%s,%s,%d,%d,true,%.0f,%s,%s,%s,%d,%.4f,%.4f,%.4f,%s\n", cand, WALL, L, N0,
							   r.secs, cls..., s2.n_isl, s2.S, s2.spread, s2.offset, basename(r.file)), CSV, "a")
		end
	end
end

# ---- the gallery, laid out like the phase space --------------------------------------
raw, hdr = readdlm(CSV, ','; header = true)
col = Dict(strip(String(hdr[i])) => i for i in eachindex(hdr))
for cand in CANDS
	panels = NamedTuple[]
	for L in reverse(LAYERS), N0 in N0S
		i = findfirst(r -> String(raw[r, col["cand"]]) == cand && raw[r, col["layers"]] == L &&
						   raw[r, col["N0"]] == N0 && raw[r, col["ok"]] == true, 1:size(raw, 1))
		i === nothing && continue
		d = load_frames(joinpath(OUT, String(raw[i, col["file"]])))
		m1 = something(findlast(iszero, d.t), 1)
		k = frame_at(d, m1, 0.2)
		st = st_at(d, k)
		push!(panels, (d = d, k = k, label = @sprintf("%d layer%s, N0 %d\n%s  S %.2f", L, L > 1 ? "s" : "", N0,
															classify_pattern(st), st.S)))
	end
	isempty(panels) && continue
	fig = plot_gallery(panels; ncol = length(N0S), views = (:top,), size_ratio = 0.55, labelsize = 16,
					   title = "$cand: coarse phase map at fB = 0.2 (rows: layers, columns: N0)")
	save(joinpath(OUT, "minimap_$(cand).png"), fig)
	println("wrote minimap_$(cand).png")
end
