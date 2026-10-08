# Stage 5a: summarise the tuning runs (reports/REPORT.md).
#
#   julia --project=. code-claude/scripts/stage5_tune_report.jl [tune] [plot] [cand1,cand2,...]
#
# Reads results/stage5/<tune>/rows_w*.csv and prints, per candidate and core, the
# class at the readout fB = 0.2 of every seed at each test point, against what the
# three-region picture asks for:
#
#   1 layer  N0 75     CENTRAL        2 layers N0 150    CENTRAL
#   1 layer  N0 1500   PATTERN        2 layers N0 2400   PATTERN
#   4 layers N0 150    POLARISED      unconfined N0 500  POLARISED
#   (2 layers N0 600, 4 layers N0 600: in between, shown for information)
#
# PATTERN = RING or ISLANDS (you: a ring counts). Writes the table to
# results/stage5/<tune>/summary.md; with "plot", also a 3D overview of seed 1
# (rows = candidates, columns = test points) as overview_<core>.png.

using Statistics, Printf, DelimitedFiles, JLD2, CairoMakie
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))      # read_rows, with_readout!, frame_at

const TAG  = length(ARGS) >= 1 ? ARGS[1] : "tune"
const PLOT = length(ARGS) >= 2 && ARGS[2] == "plot"
const ONLY = length(ARGS) >= 3 ? split(ARGS[3], ",") : nothing     # candidates to draw (default: all)
const DIR  = joinpath(@__DIR__, "..", "results", "stage5", TAG)

short(c) = Dict("CENTRAL" => "C", "POLARISED" => "P", "RING" => "R", "ISLANDS" => "I",
				"OTHER" => "o", "NONE" => "-")[String(c)]
ispattern(c) = String(c) in ("RING", "ISLANDS")

const POINTS = [(1, 75, "CENTRAL"), (2, 150, "CENTRAL"), (1, 1500, "PATTERN"), (2, 2400, "PATTERN"),
				(4, 150, "POLARISED"), (2, 600, "-"), (4, 600, "-")]
hit(c, want) = want == "PATTERN" ? ispattern(c) : String(c) == want

rows = with_readout!(read_rows(DIR))     # classes recomputed with the current classifier
ok = filter(r -> r["ok"] == true, rows)
cands = unique(String(r["cand"]) for r in ok if r["kind"] != "fusion")
cores = unique(String(r["core"]) for r in ok if r["kind"] != "fusion")

io = IOBuffer()
println(io, "# Stage 5a tuning: class at fB = 0.2 per seed (C central, R ring, I islands, P polarised, o other)\n")
println(io, "Score = runs that match the target (PATTERN = ring or islands) out of all runs at the six target points.")
println(io, "S = looseness of the B pattern at fB = 0.2 (mean fraction of B among a B cell's neighbours; lower = looser).\n")
print(io, "| candidate | core | score |")
for (L, N0, want) in POINTS; print(io, " L$L N$N0 ($(want == "-" ? "info" : lowercase(want[1:3]))) |"); end
println(io, " S (L1 N1500) | S (L2 N2400) | free: P / offset | BB / CC / BC 24 h |")
println(io, "|", repeat("---|", 3 + length(POINTS) + 4))
for cand in cands, core in cores
	sel = filter(r -> String(r["cand"]) == cand && String(r["core"]) == core, ok)
	isempty(sel) && continue
	nhit = 0; ntot = 0
	cells = String[]
	for (L, N0, want) in POINTS
		rs = sort(filter(r -> r["kind"] == "confined" && r["layers"] == L && r["N0"] == N0, sel); by = r -> r["seed"])
		cl = [String(r["class_02"]) for r in rs]
		if want != "-"
			nhit += count(c -> hit(c, want), cl); ntot += length(cl)
		end
		push!(cells, isempty(cl) ? "" : join(short.(cl), " "))
	end
	Smean(L, N0) = (v = [r["S"] for r in sel if r["kind"] == "confined" && r["layers"] == L && r["N0"] == N0];
					isempty(v) ? NaN : mean(v))
	# unconfined runs and the fusion assay have no plates: a friction variant ("x@lw3") uses those of "x"
	basec = String(split(cand, "@")[1])
	fr = filter(r -> r["kind"] == "free" && String(r["cand"]) == basec && String(r["core"]) == core, ok)
	nP = count(r -> String(r["class_02"]) == "POLARISED", fr)
	nhit += nP; ntot += length(fr)
	fus = filter(r -> r["kind"] == "fusion" && String(r["cand"]) == basec, ok)
	fstr = isempty(fus) ? "" : @sprintf("%.2f / %.2f / %.2f", fus[1]["bb24"], fus[1]["cc24"], fus[1]["bc24"])
	@printf(io, "| %s | %s | %d/%d | %s | %.2f | %.2f | %d/%d / %.2f | %s |\n", cand, core, nhit, ntot,
			join(cells, " | "), Smean(1, 1500), Smean(2, 2400), nP, length(fr),
			isempty(fr) ? NaN : mean(r["offset_run"] for r in fr), fstr)
end
nf = count(r -> r["ok"] != true, rows)
println(io, "\n$(length(ok)) runs done", nf > 0 ? ", $nf FAILED" : "", ". Thickness at the end (layers_eff / nominal): ",
		@sprintf("%.3f ± %.3f", mean(r["layers_eff"] / r["layers"] for r in ok if r["kind"] == "confined"),
				 std([r["layers_eff"] / r["layers"] for r in ok if r["kind"] == "confined"])))
txt = String(take!(io))
print(txt)
write(joinpath(DIR, "summary.md"), txt)

# ---- 3D overview of seed 1 -----------------------------------------------------------
if PLOT
	for core in cores
		panels = NamedTuple[]
		ncol = length(POINTS) + 1
		for cand in (ONLY === nothing ? cands : filter(in(ONLY), cands))
			sel = filter(r -> String(r["cand"]) == cand && String(r["core"]) == core && r["seed"] == 1, ok)
			length(sel) < 2 && continue
			for (L, N0, want) in vcat(POINTS, [(0, 500, "POLARISED")])
				kind = L == 0 ? "free" : "confined"
				i = findfirst(r -> r["kind"] == kind && r["layers"] == L && r["N0"] == N0, sel)
				if i === nothing
					push!(panels, (d = (x = [[0.0]], y = [[0.0]], z = [[0.0]], r = [[0.01]], cell_state = [[3]]), k = 1, label = "(missing)"))
					continue
				end
				d = load_frames(joinpath(DIR, String(sel[i]["file"])))
				lbl = @sprintf("%s\n%s, N0 %d: %s", cand, L == 0 ? "unconfined" : "$L layer$(L > 1 ? "s" : "")", N0,
							   short(sel[i]["class_02"]))
				push!(panels, (d = d, k = frame_at(d, 0.2), label = lbl))
			end
		end
		isempty(panels) && continue
		fig = plot_gallery(panels; ncol = ncol, views = (:top,), size_ratio = 0.45, labelsize = 13,
						   title = "stage 5a tuning, core $core, seed 1, at fB = 0.2 (top cutaways; unconfined: mid-plane cut)")
		save(joinpath(DIR, "overview_$(core).png"), fig)
		println("wrote overview_$(core).png")
	end
end
