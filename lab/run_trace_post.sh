#!/usr/bin/env bash
# VCD -> LACPo features -> block power -> cache npz -> issue npz.
# Uses the Python 3.7 env shipped in the runtime pack (sklearn 0.20.4).
#
#   export SEGSEE_REPO=$PWD
#   export DECOUPLINGCAP=/path/to/DecouplingCap
#   export SEGSEE_TRACES=/opt/segsee/traces/lacpo_cy100/bench
#   bash lab/run_trace_post.sh
#
# Then, still with SEGSEE_TRACES set, from SEGSEE_REPO:
#   /opt/segsee/lacpo37/bin/python lacpo/workload.py
#   /opt/segsee/lacpo37/bin/python lacpo/workload5.py
set -euo pipefail
ROOT=${SEGSEE_ROOT:-/opt/segsee}
O=${SEGSEE_TRACES:-$ROOT/traces/lacpo_cy100/bench}
REPO=${SEGSEE_REPO:-$(cd "$(dirname "$0")/.." && pwd)}
DEC=${DECOUPLINGCAP:?set DECOUPLINGCAP to the DecouplingCap tree that holds tools/run_lacpo_flow.py and BOOM/power_modeling}
PY=$ROOT/lacpo37/bin/python
POST_JOBS=${POST_JOBS:-3}
WL="${WL:-qsort rsort mm dhrystone spmv multiply memcpy median vvadd towers}"
test -x "$PY"
test -f "$REPO/lacpo/fast_vcd2features.py"
test -f "$DEC/tools/run_lacpo_flow.py"
feat() {
  w=$1
  if [ ! -f "$O/$w.vcd" ]; then echo "SKIP $w no vcd"; return 0; fi
  echo "FEAT $w $(date +%T)"
  ( cd "$REPO/lacpo" && "$PY" fast_vcd2features.py --vcd "$O/$w.vcd" --out "$O/${w}_features.csv" --benchmark "$w" ) \
    > "$O/$w.feat.log" 2>&1 || { echo "feat $w FAILED"; return 0; }
  ( cd "$DEC" && "$PY" tools/run_lacpo_flow.py --features-csv "$O/${w}_features.csv" --out "$O/${w}_out" --name "$w" ) \
    > "$O/$w.lacpo.log" 2>&1 || { echo "lacpo $w FAILED"; return 0; }
  ( cd "$REPO/lacpo" && "$PY" cache_events.py --vcd "$O/$w.vcd" --out "$O/${w}_cache.npz" ) \
    > "$O/$w.cache.log" 2>&1 || { echo "cache $w FAILED"; return 0; }
  ( cd "$REPO/lacpo" && "$PY" issue_trace.py --vcd "$O/$w.vcd" --out "$O/${w}_issue.npz" ) \
    > "$O/$w.issue.log" 2>&1 || { echo "issue $w FAILED"; return 0; }
  echo "POST $w $(date +%T) $(grep 'composed P' "$O/$w.lacpo.log" || true)"
}
export -f feat
export O REPO DEC PY
printf '%s\n' $WL | xargs -P "$POST_JOBS" -I{} bash -c 'feat "$@"' _ {}
echo POST_DONE
