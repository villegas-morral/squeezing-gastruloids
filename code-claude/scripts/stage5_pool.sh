#!/usr/bin/env bash
# Keep the machine busy with stage5_run.jl queue workers (reports/REPORT.md, stage 5).
#
#   bash code-claude/scripts/stage5_pool.sh [n_workers] [first_id]   (from the repo root; default 24, 100)
#
# Reads code-claude/results/stage5/pool.txt every minute: one line per job list,
# "<jobs.csv> <out>", in priority order (lines can be added while it runs). While
# fewer than n_workers Julia processes are running, it starts one more queue worker
# on the first list that still has unclaimed jobs (one per minute, since a worker
# needs a few minutes to load). It stops when every list is fully claimed.

N=${1:-24}
S5=code-claude/results/stage5
POOL=$S5/pool.txt
ID=${2:-100}

unclaimed() {   # number of jobs in $1 not claimed yet in results/stage5/$2/claims
	local n=0
	while IFS= read -r line; do
		line=${line%$'\r'}
		[[ -z "$line" || "$line" == \#* ]] && continue
		[ -d "$S5/$2/claims/${line//,/_}" ] || n=$((n+1))
	done < "$1"
	echo $n
}

while true; do
	running=$(powershell -NoProfile -Command "(Get-Process julia -ErrorAction SilentlyContinue | Measure-Object).Count" | tr -d '\r')
	target=""
	left=0
	while read -r jobs out; do
		[[ -z "$jobs" || "$jobs" == \#* ]] && continue
		u=$(unclaimed "$jobs" "$out")
		left=$((left+u))
		# a worker claims its first job only after loading: count the ones still loading
		if [ -z "$target" ] && [ "$u" -gt 0 ]; then target="$jobs $out"; fi
	done < "$POOL"
	if [ "$left" -eq 0 ]; then echo "$(date +%H:%M) pool: every job claimed, stopping"; break; fi
	if [ "$running" -lt "$N" ] && [ -n "$target" ]; then
		set -- $target
		ID=$((ID+1))
		nohup julia --project=. code-claude/scripts/stage5_run.jl "$1" $ID 0 "$2" > "$S5/logs/pool_w$ID.log" 2>&1 &
		echo "$(date +%H:%M) pool: $running running, started worker $ID on $1 ($left jobs unclaimed)"
	fi
	sleep 60
done
