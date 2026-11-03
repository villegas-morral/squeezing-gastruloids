# Division time 13 h for B and C (November 2026 test): compare with the final runs.
#
#   julia --project=. code-claude/scripts/tdiv_report.jl
#
# Reads results/stage5/tdiv (tau_div = 195 = 13 h for every cell, seed 1) and the same
# (layers, N0, seed 1) points of results/stage5/phase_D (tau_div B 300 = 20 h, C 750 =
# 50 h). Writes results/stage5/tdiv/compare.txt (classes at fB 0.2 and 0.12 with the
# figures' readout and the October one, readout time, growth) and a 3D gallery of the
# 13 h runs at fB 0.2, laid out like the phase space -> claude-answers/img/05-tdiv-*.png.

using Statistics, Printf, DelimitedFiles, JLD2, LinearAlgebra
using GLMakie; GLMakie.activate!(visible = false, scalefactor = 0.3, px_per_unit = 1)   # transparency (GHOST)
include(joinpath(@__DIR__, "gastruloids_analysis.jl"))
include(joinpath(@__DIR__, "stage5_common.jl"))

const S5  = joinpath(@__DIR__, "..", "results", "stage5")
const IMG = joinpath(@__DIR__, "..", "..", "claude-answers", "img")
const CAND = "july_frB4.5_frC6_muB1.9"

new = with_readout!(filter(r -> r["ok"] == true, read_rows(joinpath(S5, "tdiv"))))
old = with_readout!(filter(r -> r["ok"] == true && r["seed"] == 1 &&
						   ((r["kind"] == "confined" && r["cand"] == CAND * "@lw1.5") || (r["kind"] == "free" && r["cand"] == CAND)),
						   read_rows(joinpath(S5, "phase_D"))))
key(r) = (r["kind"] == "free" ? 0.0 : Float64(r["layers"]), Int(r["N0"]))
oldby = Dict(key(r) => r for r in old)
ab(c) = Dict("CENTRAL" => "central", "PATTERN" => "pattern", "POLARISED" => "polarised", "OTHER" => "other")[phase_class(String(c))]

lines_ = String[]
push!(lines_, @sprintf("%-11s %5s | %-10s %-10s | %-10s %-10s | %-10s %-10s | %6s %6s | %6s %6s",
					   "layers", "N0", "now .2", "13h .2", "now .12", "13h .12", "Oct now", "Oct 13h", "t02 h", "13h", "N/N0", "13h"))
for r in sort(new; by = key)
	k = key(r); o = get(oldby, k, nothing); o === nothing && continue
	push!(lines_, @sprintf("%-11s %5d | %-10s %-10s | %-10s %-10s | %-10s %-10s | %6.1f %6.1f | %6.2f %6.2f",
		k[1] == 0 ? "unconfined" : @sprintf("%g", k[1]), k[2],
		ab(o["class_exp10_02"]), ab(r["class_exp10_02"]), ab(o["class_exp10_012"]), ab(r["class_exp10_012"]),
		ab(o["class_02"]), ab(r["class_02"]), (o["t_02"] - 25) / 15, (r["t_02"] - 25) / 15,
		o["N_end"] / o["N0"], r["N_end"] / r["N0"]))
end
push!(lines_, "", "N/N0 is at the end of the run (fB = 0.10); t02 h = hours from the start of differentiation (25 time units after seeding) to fB = 0.2.")
write(joinpath(S5, "tdiv", "compare.txt"), join(lines_, "\n") * "\n")
println(join(lines_, "\n"))

# gallery of the 13 h runs at fB = 0.2, rows = layers (thick at the top), columns = N0
Ls = sort(unique(key(r)[1] for r in new); rev = true); Ns = sort(unique(key(r)[2] for r in new))
fig = Figure(size = (300 * length(Ns), 300 * length(Ls)), figure_padding = 10)
for (i, L) in enumerate(Ls), (j, N0) in enumerate(Ns)
	rs = filter(r -> key(r) == (L, N0), new); isempty(rs) && continue
	r = rs[1]; d = load_frames(joinpath(r["dir"], String(r["file"]))); k = frame_at(d, 0.2)
	c = (x = d.x[k], y = d.y[k], z = d.z[k], r = d.r[k], cell_state = d.cell_state[k])
	draw_cells!(fig[i, j], c; view = :top, labelsize = 14,
				label = @sprintf("%s, N0 %d|%s", L == 0 ? "unconfined" : @sprintf("%g layer%s", L, L == 1 ? "" : "s"), N0,
								 ab(r["class_exp10_02"])) |> s -> replace(s, "|" => "\n"))
end
save(joinpath(IMG, "05-tdiv-gallery.png"), fig)
println("wrote ", joinpath(IMG, "05-tdiv-gallery.png"))
