# Shared by the stage 5 and 6 analysis scripts (include after gastruloids_analysis.jl).
#
# read_rows(dir)        every run of a results/stage5/<dir> folder, once (rows_w*.csv)
# with_readout!(rows)   the readout of each run recomputed from its saved trajectory
#                       with the CURRENT classify_pattern
# with_diameter!(rows)  the pancake diameter at squeezing of each confined run
#
# Why recompute: stage5_run.jl stores the classes of the classifier its worker had
# loaded, and the classifier changed on Oct 7 (the `fill` test that recognises a
# loose central cluster). Recomputing makes every run comparable. The results are
# cached in <dir>/readout_cache.csv, so each trajectory is read only once per
# classifier version.

const FBS = (0.4, 0.3, 0.2, 0.15, 0.12)
const CLS_KEYS = ["class_04", "class_03", "class_02", "class_015", "class_012"]
const ST_KEYS = ["n_isl", "f_isl", "f_iso", "f_max", "S", "spread", "fill", "offset", "offset_max"]
# the experiment-like readout (classify_exp, November 2026) at the same fB, with the
# experiments' minimum cluster of 20 cells and with 10 (in between: small aggregates)
const EXP_KEYS = ["class_exp_04", "class_exp_03", "class_exp_02", "class_exp_015", "class_exp_012"]
const EXP10_KEYS = ["class_exp10_02", "class_exp10_012"]
# the readout of the figures since Nov 2 (second round): clusters of >= 10 cells and a
# core cell needs at least half B neighbours (the experiments: 40 %), which splits rings
# that hang together by thin B threads -- in between the October readout and the
# experiments' (claude-answers/07)
const CORE_FRAC = 0.5
const EXP50_KEYS = ["class_exp50_02", "class_exp50_015", "class_exp50_012", "class_exp50_010"]
const EXP50_FBS = (0.2, 0.15, 0.12, 0.10)
# the figures' readout since Nov 3: at fB = 0.12 (about 24 h after squeezing, as the experiments'
# 96 h with ~12 % T+), with three outcomes only. It needs at least half of the B cells in clusters
# (READ_MIN_CL); below that B is differentiating away -- the islands shrink below 10 cells and
# one survivor would read as a B-centroid -- and the run takes its outcome at fB 0.15, else 0.2.
# class_read_late: the same from fB 0.10 (the end of the run), to break ties between seeds.
const READ_KEYS = ["class_read", "class_read_late"]
const F50_KEYS = ["f_cl50_02", "f_cl50_015", "f_cl50_012", "f_cl50_010"]
const READ_MIN_CL = 0.5
const EXP_ST_KEYS = ["n_cl", "N_eff", "f_cl", "P", "n_cl10", "N_eff10", "n_cl50", "N_eff50", "f_cl50_02", "f_cl50_015", "f_cl50_012", "f_cl50_010"]
const STRING_KEYS = Set(vcat(CLS_KEYS, EXP_KEYS, EXP10_KEYS, EXP50_KEYS, READ_KEYS))
const READOUT_VERSION = "2026-11-03-read012b"

"""
The outcome read from `classes` (latest first) with `f_cl` the share of B in clusters at each:
the first one with f_cl >= READ_MIN_CL (and a cluster at all), else the last (earliest) one.
"""
function read_outcome(classes, f_cl)
	i = findfirst(j -> classes[j] != "NONE" && f_cl[j] >= READ_MIN_CL, eachindex(classes))
	classes[something(i, lastindex(classes))]
end

jobkey(r) = join((r["kind"], r["cand"], r["core"], r["layers"], r["N0"], r["seed"]), ",")

"""
	read_rows(dir)

All rows of the rows_w*.csv in `dir` as Dicts (keys from the CSV header, plus
"dir"). A job that ran twice is kept once (its first good row); a failure that a
rerun fixed is dropped.
"""
function read_rows(dir)
	rows = Dict{String,Any}[]; seen = Set{String}()
	for f in sort(readdir(dir))
		(startswith(f, "rows_w") && endswith(f, ".csv")) || continue
		raw, hdr = readdlm(joinpath(dir, f), ','; header = true)
		for i in 1:size(raw, 1)
			r = Dict{String,Any}(strip(String(hdr[j])) => raw[i, j] for j in eachindex(hdr))
			for k in ("kind", "cand", "core", "file"); r[k] = String(string(r[k])); end
			r["dir"] = dir
			if r["ok"] == true
				jobkey(r) in seen && continue
				push!(seen, jobkey(r))
			end
			push!(rows, r)
		end
	end
	filter!(r -> r["ok"] == true || !(jobkey(r) in seen), rows)
	rows
end

"Frame where the B fraction first drops to f, after seeding (t reset to 0 there)."
function frame_at(d, f)
	m1 = something(findlast(iszero, d.t), 1)
	f >= 0.6 && return m1
	k = findfirst(i -> i >= m1 && count(==(2), d.cell_state[i]) / length(d.cell_state[i]) <= f, eachindex(d.x))
	k === nothing ? length(d.x) : k
end

"The readout of one saved run: classes at FBS, stats at fB = 0.2, and the largest B offset."
function readout(file; confined = true)
	d = load_frames(file)
	out = Dict{String,Any}()
	for (key, xkey, f) in zip(CLS_KEYS, EXP_KEYS, FBS)
		k = frame_at(d, f)
		out[key] = classify_pattern(pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = confined))
		out[xkey] = classify_exp(exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]))
	end
	for (key, f) in zip(EXP10_KEYS, (0.2, 0.12))
		k = frame_at(d, f)
		out[key] = classify_exp(exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; nmin = 10))
	end
	for (key, fkey, f) in zip(EXP50_KEYS, F50_KEYS, EXP50_FBS)
		k = frame_at(d, f)
		e = exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; nmin = 10, frac = CORE_FRAC)
		out[key] = classify_exp(e); out[fkey] = e.f_cl
	end
	cl(fbs...) = [out["class_exp50_$(f)"] for f in fbs]; fc(fbs...) = [out["f_cl50_$(f)"] for f in fbs]
	out["class_read"] = read_outcome(cl("012", "015", "02"), fc("012", "015", "02"))
	out["class_read_late"] = read_outcome(cl("010", "012", "015", "02"), fc("010", "012", "015", "02"))
	k = frame_at(d, 0.2)
	st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = confined)
	for key in ST_KEYS; out[key] = getfield(st, Symbol(key)); end
	e20 = exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]); e10 = exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; nmin = 10)
	out["n_cl"] = e20.n_cl; out["N_eff"] = e20.N_eff; out["f_cl"] = e20.f_cl; out["P"] = e20.P
	out["n_cl10"] = e10.n_cl; out["N_eff10"] = e10.N_eff
	e50 = exp_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; nmin = 10, frac = CORE_FRAC)
	out["n_cl50"] = e50.n_cl; out["N_eff50"] = e50.N_eff
	out["t_02"] = d.t[k]
	m1 = frame_at(d, 0.6)
	off = 0.0
	for j in m1:length(d.x)
		x, y, z, s = d.x[j], d.y[j], d.z[j], d.cell_state[j]
		b = findall(==(2), s); isempty(b) && continue
		dz = confined ? 0.0 : 1.0
		R = quantile(sqrt.((x .- mean(x)).^2 .+ (y .- mean(y)).^2 .+ dz .* (z .- mean(z)).^2), 0.98)
		off = max(off, sqrt((mean(x[b]) - mean(x))^2 + (mean(y[b]) - mean(y))^2 + (mean(z[b]) - mean(z))^2) / R)
	end
	out["offset_run"] = off
	out
end

"""
	with_readout!(rows)

Overwrite the stored classes and stats of every confined / free run in `rows` with
readout() of its trajectory, using (and extending) each folder's readout_cache.csv.
"""
function with_readout!(rows)
	keys_ = vcat(CLS_KEYS, EXP_KEYS, EXP10_KEYS, EXP50_KEYS, READ_KEYS, ST_KEYS, EXP_ST_KEYS, ["t_02", "offset_run"])
	for dir in unique(r["dir"] for r in rows)
		cfile = joinpath(dir, "readout_cache.csv")
		cache = Dict{String,Dict{String,Any}}()
		if isfile(cfile)
			lines_ = readlines(cfile)
			if !isempty(lines_) && lines_[1] == "# $(READOUT_VERSION)"
				hdr = split(lines_[2], ',')
				for l in lines_[3:end]
					p = split(l, ',')
					cache[p[1]] = Dict{String,Any}(hdr[i] => (hdr[i] in STRING_KEYS ? String(p[i]) : parse(Float64, p[i]))
												   for i in 2:length(hdr))
				end
			end
		end
		todo = [r for r in rows if r["dir"] == dir && r["kind"] in ("confined", "free") && r["ok"] == true &&
							   !haskey(cache, r["file"])]
		isempty(todo) || println("reading $(length(todo)) new trajectories in $(basename(dir)) ...")
		res = Vector{Any}(undef, length(todo))      # in parallel when julia runs with threads (-t)
		Threads.@threads for i in eachindex(todo)
			res[i] = readout(joinpath(dir, todo[i]["file"]); confined = todo[i]["kind"] == "confined")
		end
		for (r, o) in zip(todo, res); cache[r["file"]] = o; end
		open(cfile, "w") do io
			println(io, "# $(READOUT_VERSION)")
			println(io, join(vcat("file", keys_), ","))
			for (f, c) in cache
				println(io, join(vcat(f, [c[k] isa AbstractString ? c[k] : @sprintf("%.5g", c[k]) for k in keys_]), ","))
			end
		end
		for r in rows
			(r["dir"] == dir && haskey(cache, r["file"])) || continue
			merge!(r, cache[r["file"]])
		end
	end
	rows
end

"""
	with_diameter!(rows)

Add "D" to every confined run in `rows`: the pancake diameter at squeezing (the
seeding frame), in cell diameters, 2 sqrt(2 <rho^2>) + 2r over 2r, where rho is the
in-plane distance of the cell centres from the aggregate centre (exact for a
uniform disc). Cached in <dir>/diameter_cache.csv.
"""
function with_diameter!(rows)
	for dir in unique(r["dir"] for r in rows)
		cfile = joinpath(dir, "diameter_cache.csv")
		cache = Dict{String,Float64}()
		if isfile(cfile)
			for l in readlines(cfile)[2:end]
				p = split(l, ','); cache[p[1]] = parse(Float64, p[2])
			end
		end
		todo = [r for r in rows if r["dir"] == dir && r["kind"] == "confined" && r["ok"] == true && !haskey(cache, r["file"])]
		isempty(todo) || println("measuring $(length(todo)) pancake diameters in $(basename(dir)) ...")
		for r in todo
			d = load_frames(joinpath(dir, r["file"])); k = frame_at(d, 0.6)
			x, y, rc = d.x[k], d.y[k], mean(d.r[k])
			cache[r["file"]] = (2sqrt(2mean((x .- mean(x)).^2 .+ (y .- mean(y)).^2)) + 2rc) / 2rc
		end
		open(cfile, "w") do io
			println(io, "file,D")
			for (f, v) in cache; println(io, f, ",", @sprintf("%.4f", v)); end
		end
		for r in rows
			(r["dir"] == dir && haskey(cache, r["file"])) && (r["D"] = cache[r["file"]])
		end
	end
	rows
end
