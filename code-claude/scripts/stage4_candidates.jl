# Parameter candidates for stage 4 (reports/REPORT.md), shared by stage4_screen.jl and
# stage4_minimap.jl. Include after gastruloids_lib.jl (needs PSETS).
#
# A candidate is a base set from PSETS with some BB / CC entries replaced. The BC
# entries stay at the FUSION values (frozen) -- `variant` refuses anything else.
# "name@w0.1" means candidate `name` with a diluted B core (initialize_diff! width 0.1).

function variant(base; frep_BB = nothing, frep_CC = nothing, fatt_BB = nothing,
		fatt_CC = nothing, mu_BB = nothing, mu_CC = nothing)
	P = PSETS[base]
	mu = float.(copy(P.mu)); fa = float.(copy(P.f_att)); fr = float.(copy(P.f_rep))
	frep_BB === nothing || (fr[2,2] = frep_BB)
	frep_CC === nothing || (fr[3,3] = frep_CC)
	fatt_BB === nothing || (fa[2,2] = fatt_BB)
	fatt_CC === nothing || (fa[3,3] = fatt_CC)
	mu_BB   === nothing || (mu[2,2] = mu_BB)
	mu_CC   === nothing || (mu[3,3] = mu_CC)
	@assert fa[2,3] == 1.6 && fr[2,3] == 3.0 && mu[2,3] == 2.0 "BC entries must stay frozen"
	(mu = mu, f_att = fa, f_rep = fr)
end

const CANDIDATES = Dict(
	# wave 1
	"july"              => () -> variant("july"),
	"fusion"            => () -> variant("fusion"),
	"cp"                => () -> variant("confined_patterns"),
	"july_frB4.5"       => () -> variant("july"; frep_BB = 4.5),
	"july_frB4.5_frC6"  => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0),
	"july_fa0.7"        => () -> variant("july"; fatt_BB = 0.7*6.4, fatt_CC = 0.7*3.0),
	"july_fa0.5"        => () -> variant("july"; fatt_BB = 0.5*6.4, fatt_CC = 0.5*3.0),
	"cp_frB4.5"         => () -> variant("confined_patterns"; frep_BB = 4.5),
	# wave 2: around the FUSION set (stiff BB and CC), which gave sparse patterns,
	# polarisation and BC glue in wave 1
	"fusion_frC5.1"     => () -> variant("fusion"; frep_CC = 5.1),
	"fusion_frB9"       => () -> variant("fusion"; frep_BB = 9.0),
	"fusion_muB2.2"     => () -> variant("fusion"; mu_BB = 2.2),
	"fusion_muC2.0"     => () -> variant("fusion"; mu_CC = 2.0),
	# stage 5 tuning (Oct 7): you picked july_frB4.5, "a little bit more loose would
	# also be ok" -> BB between 4.5 and 6 (FUSION), CC 5.1 (July) or 6 (FUSION).
	# fusion_frC5.1 above is the same as "july_frB6".
	"july_frB5"         => () -> variant("july"; frep_BB = 5.0),
	"july_frB5.5"       => () -> variant("july"; frep_BB = 5.5),
	"july_frB5_frC6"    => () -> variant("july"; frep_BB = 5.0, frep_CC = 6.0),
	# Oct 8: patterns in slabs of 3+ layers? (the July A,B,C runs had them at 3 and 4)
	# single changes to the final set's BB / CC entries
	"july_frB3_frC6"            => () -> variant("july"; frep_BB = 3.0, frep_CC = 6.0),
	"july_frB4.5_frC6_faC3.6"   => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, fatt_CC = 3.6),
	"july_frB4.5_frC6_faB5.6"   => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, fatt_BB = 5.6),
	"july_frB4.5_frC6_muB1.8"   => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.8),
	# mu BB 1.8 gave a ring at 3 layers, N0 2400: a shorter B-B range. In between, and the
	# same range with the B-B adhesion energy of the final set kept (f_att BB 6.4 -> 11.1)
	"july_frB4.5_frC6_muB1.9"   => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.9),
	"july_frB4.5_frC6_muB1.85"  => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.85),
	"july_frB4.5_frC6_muB1.8_faB11" => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.8, fatt_BB = 11.1),
	"july_frB4.5_frC6_muB1.875" => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.875),
	# Nov 2: with the experiment-like readout the rings of 3+ layers count as central until
	# they break. Single changes to the final set that should split the B core sooner:
	# weaker B-B adhesion (energy W_BB 0.79 -> 0.72), stiffer B, more cohesive C (W_CC 0.65 -> 0.72)
	"july_frB4.5_frC6_muB1.9_faB5.8" => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.9, fatt_BB = 5.8),
	"july_frB4.5_frC6_muB1.9_faB6.1" => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.9, fatt_BB = 6.1),   # Nov 3: half-way to 5.8
	"july_frB5.5_frC6_muB1.9"        => () -> variant("july"; frep_BB = 5.5, frep_CC = 6.0, mu_BB = 1.9),
	"july_frB4.5_frC6_muB1.9_faC3.3" => () -> variant("july"; frep_BB = 4.5, frep_CC = 6.0, mu_BB = 1.9, fatt_CC = 3.3),
)

"""
	parse_cand("name@w0.1@lw10") -> ("name", 0.1, 10.0)

`@wX`  diluted B core of width X (initialize_diff!)
`@lwX` friction with the plates, lambda_wall = X (run_confined)
Both optional, in any order.
"""
function parse_cand(cand)
	parts = split(cand, "@")
	width = 0.0; lw = 0.0
	for p in parts[2:end]
		startswith(p, "lw") ? (lw = parse(Float64, p[3:end])) :
		startswith(p, "w")  ? (width = parse(Float64, p[2:end])) : error("unknown suffix @$p")
	end
	(String(parts[1]), width, lw)
end
