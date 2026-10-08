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
const READOUT_VERSION = "2026-10-07-fill-single"

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
	for (key, f) in zip(CLS_KEYS, FBS)
		k = frame_at(d, f)
		out[key] = classify_pattern(pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = confined))
	end
	k = frame_at(d, 0.2)
	st = pattern_stats(d.x[k], d.y[k], d.z[k], d.cell_state[k]; confined = confined)
	for key in ST_KEYS; out[key] = getfield(st, Symbol(key)); end
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
	keys_ = vcat(CLS_KEYS, ST_KEYS, ["t_02", "offset_run"])
	for dir in unique(r["dir"] for r in rows)
		cfile = joinpath(dir, "readout_cache.csv")
		cache = Dict{String,Dict{String,Any}}()
		if isfile(cfile)
			lines_ = readlines(cfile)
			if !isempty(lines_) && lines_[1] == "# $(READOUT_VERSION)"
				hdr = split(lines_[2], ',')
				for l in lines_[3:end]
					p = split(l, ',')
					cache[p[1]] = Dict{String,Any}(hdr[i] => (i - 1 <= length(CLS_KEYS) ? String(p[i]) : parse(Float64, p[i]))
												   for i in 2:length(hdr))
				end
			end
		end
		todo = [r for r in rows if r["dir"] == dir && r["kind"] in ("confined", "free") && r["ok"] == true &&
							   !haskey(cache, r["file"])]
		isempty(todo) || println("reading $(length(todo)) new trajectories in $(basename(dir)) ...")
		for r in todo
			cache[r["file"]] = readout(joinpath(dir, r["file"]); confined = r["kind"] == "confined")
		end
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
