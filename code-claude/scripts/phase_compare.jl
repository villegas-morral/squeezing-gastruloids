# Full phase spaces (the points and seeds of the current set's, phase_D) of the candidate sets of
# Nov 3-4 2026 next to the current set's, with the figures' readout (class_read; ties between seeds
# broken by class_read_late), on the same points (claude-answers/10 §2.6).
#
#   julia --project=. -t 8 code-claude/scripts/phase_compare.jl F        # now vs 6.1 + 13 h
#   julia --project=. -t 8 code-claude/scripts/phase_compare.jl F G H    # now vs all three
#
# Sets: D now (phase_D + tune), F f_att BB 6.1 + tau_div 13 h (phase_F), G f_att BB 7.5 + CC 3.6
# (phase_G + the seed-1 points of grid_nov3b/c), H f_att BB 7 + CC 3.3 (phase_H + grid_nov3b/c).
# Writes, with <tag> = the letters given (e.g. "F", "FGH"):
#   claude-answers/img/10-phase<tag>.png        the charts, with the pancake-area lines
#   claude-answers/img/10-phase<tag>-thick.png  readout frames of the thick-slab points, seed 1
#   results/stage5/phase_<last letter>/compare_<tag>.txt   every seed per point; thresholds; times

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using CairoMakie                 # everything here is opaque (confined panels, cut at the mid-plane)
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5  = joinpath(@__DIR__, "..", "results", "stage5")
const IMG = joinpath(@__DIR__, "..", "..", "claude-answers", "img")
const B   = "july_frB4.5_frC6_muB1.9"
# (name, title, folders, confined candidate, unconfined candidate)
const ALL = Dict(
	"D" => ("now", L"\mathrm{now}\ (f_\mathrm{att}^{BB} = 6.4)", ["phase_D", "tune"], B * "@lw1.5", B),
	"F" => ("6.1, 13 h", L"f_\mathrm{att}^{BB} = 6.1,\ \tau_\mathrm{div} = 13\,\mathrm{h}", ["phase_F"],
			B * "_faB6.1@lw1.5@tdiv195", B * "_faB6.1@tdiv195"),
	"G" => ("7.5/3.6", L"f_\mathrm{att}^{BB} = 7.5,\ f_\mathrm{att}^{CC} = 3.6", ["phase_G", "grid_nov3b", "grid_nov3c"],
			B * "_faB7.5_faC3.6@lw1.5", B * "_faB7.5_faC3.6"),
	"H" => ("7/3.3", L"f_\mathrm{att}^{BB} = 7,\ f_\mathrm{att}^{CC} = 3.3", ["phase_H", "grid_nov3b", "grid_nov3c"],
			B * "_faB7_faC3.3@lw1.5", B * "_faB7_faC3.3"))
const KEYS = vcat("D", ARGS)
const TAG  = join(ARGS)

pt(r) = (Float64(r["layers"]), Int(r["N0"]))
cls(r) = phase_class(String(r["class_read"]))
late(r) = phase_class(String(r["class_read_late"]))
abbr = Dict("POLARISED" => "o", "CENTRAL" => "c", "PATTERN" => "X", "OTHER" => "n")

cache = Dict{String,Vector}()
folder(d) = get!(() -> isdir(joinpath(S5, d)) ? with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, d)))) : Dict{String,Any}[], cache, d)
# one row per (point, seed): the first folder that has it wins
function dedup(rs, key)
	seen = Set(); out = similar(rs, 0)
	for r in rs
		k = (key(r), r["seed"]); k in seen && continue
		push!(seen, k); push!(out, r)
	end
	out
end
rows = Dict{String,Any}()
for k in KEYS
	_, _, dirs, cc, fc = ALL[k]
	rs = vcat([folder(d) for d in dirs]...)
	rows[k] = (conf = dedup(filter(r -> r["kind"] == "confined" && String(r["cand"]) == cc, rs), pt),
			   free = dedup(filter(r -> r["kind"] == "free" && String(r["cand"]) == fc, rs), r -> r["N0"]))
end
# the points of the current set's map, and the unconfined N0 values every set has
P0 = Set(pt.(rows["D"].conf))
Nf = intersect([Set(r["N0"] for r in rows[k].free) for k in KEYS]...)
extra = Dict(k => filter(r -> !(pt(r) in P0), rows[k].conf) for k in KEYS)       # e.g. N0 9600
for k in KEYS
	rows[k] = (conf = with_diameter!(filter(r -> pt(r) in P0, rows[k].conf)), free = filter(r -> r["N0"] in Nf, rows[k].free))
end

# the majority of a point; a tie goes to the tied class most frequent later (as phase_space! does)
function majority(rs)
	isempty(rs) && return ""
	order = ["POLARISED", "CENTRAL", "PATTERN", "OTHER"]
	n = [count(r -> cls(r) == c, rs) for c in order]
	tied = findall(==(maximum(n)), n)
	length(tied) == 1 && return order[tied[1]]
	l = [count(r -> late(r) == order[t], rs) for t in tied]
	order[tied[argmax(l)]]
end
at(k, p) = sort(filter(r -> pt(r) == p, rows[k].conf); by = r -> r["seed"])
atf(k, N) = sort(filter(r -> r["N0"] == N, rows[k].free); by = r -> r["seed"])

# ---- the charts --------------------------------------------------------------------------
nc = length(KEYS) <= 2 ? length(KEYS) : 2
nr = cld(length(KEYS), nc)
fig = Figure(size = (950 * nc, 860 * nr + 90), figure_padding = 20, fonts = TEX_FONTS)
for (j, k) in enumerate(KEYS)
	R = rows[k]
	cr = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), class = String(r["class_read"]), late = String(r["class_read_late"])) for r in R.conf]
	fr = [(N0 = Int(r["N0"]), class = String(r["class_read"]), late = String(r["class_read_late"])) for r in R.free]
	cont = [(layers = Float64(r["layers"]), N0 = Int(r["N0"]), A = pi * r["D"]^2 * 25) for r in R.conf]
	phase_space!(fig[1 + (j - 1) ÷ nc, 1 + (j - 1) % nc], cr; free = fr, title = ALL[k][2], contours = cont)
end
els = Any[PolyElement(color = PHASE_COLORS[c]) for c in ("POLARISED", "CENTRAL", "PATTERN")]
push!(els, LineElement(color = (:black, 0.7), linewidth = 3))
Legend(fig[nr + 1, 1:nc], els, vcat([PHASE_NAMES[c] for c in ("POLARISED", "CENTRAL", "PATTERN")], "pancake area at squeezing (µm²)");
	   orientation = :horizontal, framevisible = false, labelsize = 24)
save(joinpath(IMG, "10-phase$(TAG).png"), fig)

# ---- per point ---------------------------------------------------------------------------
names_ = join([ALL[k][1] for k in KEYS], " | ")
lines_ = String["o polarised, c B-centroid, X clusters, n other; one letter per seed; sets: $names_;",
				"* = the majority differs from now's", ""]
for L in sort(unique(first.(collect(P0))))
	parts = String[]
	for N in sort(unique(last.(filter(p -> p[1] == L, collect(P0)))))
		m0 = majority(at("D", (L, N)))
		ch = any(k -> majority(at(k, (L, N))) != m0, KEYS[2:end]) ? "*" : ""
		push!(parts, @sprintf("%d:%s%s", N, join([join(abbr[cls(r)] for r in at(k, (L, N))) for k in KEYS], "|"), ch))
	end
	push!(lines_, @sprintf("%-3g %s", L, join(parts, "  ")))
end
push!(lines_, "free " * join([@sprintf("%d:%s", N, join([join(abbr[cls(r)] for r in atf(k, N)) for k in KEYS], "|")) for N in sort(collect(Nf))], "  "))
for k in KEYS[2:end]
	n = count(p -> majority(at("D", p)) != majority(at(k, p)), collect(P0))
	push!(lines_, "$(ALL[k][1]): $n of $(length(P0)) points change their majority")
end
# where clusters start and how far polarisation reaches, per number of layers (majorities)
push!(lines_, "", "first N0 with a cluster majority (and above it all clusters?) / last N0 polarised before the first B-centroid:")
for L in sort(unique(first.(collect(P0))))
	Ns = sort(unique(last.(filter(p -> p[1] == L, collect(P0)))))
	s = String[]
	for k in KEYS
		m = [majority(at(k, (L, N))) for N in Ns]
		i = findfirst(==("PATTERN"), m); j = findfirst(==("CENTRAL"), m)
		pol = j === nothing ? "all" : j == 1 ? "-" : string(Ns[j - 1])
		push!(s, @sprintf("%s %s/%s", ALL[k][1], i === nothing ? ">" * string(Ns[end]) : string(Ns[i]) * (all(==("PATTERN"), m[i:end]) ? "" : "~"), pol))
	end
	push!(lines_, @sprintf("  %-3g %s", L, join(s, "   ")))
end
# readout time, thickness, run time
for k in KEYS
	R = rows[k]
	th = [r["layers_eff"] / r["layers"] for r in R.conf if r["layers"] <= 3]
	push!(lines_, @sprintf("%s: %d confined + %d free runs; readout fB 0.2 at %.1f h on average (%.1f-%.1f); layers_eff/layers (1-3 layers) %.2f-%.2f; run time median %.0f s, max %.0f s",
		ALL[k][1], length(R.conf), length(R.free), mean((r["t_02"] - 25) / 15 for r in R.conf), minimum((r["t_02"] - 25) / 15 for r in R.conf),
		maximum((r["t_02"] - 25) / 15 for r in R.conf), minimum(th), maximum(th),
		median([r["secs"] for r in R.conf]), maximum(r["secs"] for r in R.conf)))
end
for k in KEYS, r in sort(extra[k]; by = r -> (pt(r), r["seed"]))
	push!(lines_, @sprintf("extra point %s: %g layers, N0 %d, seed %d: %s (%.0f s)", ALL[k][1], r["layers"], r["N0"], r["seed"], cls(r), r["secs"]))
end
out = joinpath(S5, "phase_" * string(last(KEYS)), "compare_$(TAG).txt")
write(out, join(lines_, "\n") * "\n")
println(join(lines_, "\n"))

# ---- readout frames of the thick-slab points (seed 1, every set) ----------------------------
# the frame the readout used: fB 0.12 if half of B is in clusters there, else 0.15, else 0.2
read_fb(r) = r["f_cl50_012"] >= READ_MIN_CL ? 0.12 : r["f_cl50_015"] >= READ_MIN_CL ? 0.15 : 0.2
pts = [(3, 600), (3, 1200), (3, 2400), (4, 2400), (4, 4800), (6, 2400)]
g = Figure(size = (360 * length(pts), 430 * length(KEYS)), figure_padding = (15, 15, 15, 60))
for (j, p) in enumerate(pts), (i, k) in enumerate(KEYS)
	rs = at(k, p); isempty(rs) && continue
	r = rs[1]; d = load_frames(joinpath(r["dir"], String(r["file"]))); fb = read_fb(r); kf = frame_at(d, fb)
	c = (x = d.x[kf], y = d.y[kf], z = d.z[kf], r = d.r[kf], cell_state = d.cell_state[kf])
	draw_cells!(g[i, j], c; view = :top, labelsize = 15,
				label = @sprintf("%s\n%g layers, N0 %d, seed %d\n%s (read at fB %.2f)", ALL[k][1], p[1], p[2], r["seed"],
								 PHASE_NAMES[cls(r)], fb))
end
save(joinpath(IMG, "10-phase$(TAG)-thick.png"), g)
println("wrote ", IMG)
