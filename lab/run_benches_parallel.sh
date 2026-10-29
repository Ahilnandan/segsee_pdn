#!/usr/bin/env bash
# One benchmark per core on the MediumBoom debug simulator.
# The binary is not built with verilator --threads. JOBS processes is the parallelism that works.
#
#   sudo MAX_CYCLES=5000000 JOBS=$(nproc) bash lab/run_benches_parallel.sh
#   MAX_CYCLES=200000 JOBS=$(nproc) bash lab/run_benches_parallel.sh
#
# Writes $SEGSEE_TRACES/<w>.vcd, <w>.sim.log, <w>.trace.gz
# "*** FAILED *** via trace_count" in the log is the cycle cap, not a crashed core.
set -euo pipefail
ROOT=${SEGSEE_ROOT:-/opt/segsee}
SIM=$ROOT/chipyard-cy100/sims/verilator/simulator-example-MediumBoomConfig-debug
B=$ROOT/chipyard-cy100/riscv-tools-install/riscv64-unknown-elf/share/riscv-tests/benchmarks
O=${SEGSEE_TRACES:-$ROOT/traces/lacpo_cy100/bench}
MAX_CYCLES=${MAX_CYCLES:-5000000}
JOBS=${JOBS:-$(nproc)}
WL="${WL:-qsort rsort mm dhrystone spmv multiply memcpy median vvadd towers}"
export LD_LIBRARY_PATH=$ROOT/chipyard-cy100/riscv-tools-install/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
mkdir -p "$O"
test -x "$SIM"
cd "$O"
echo "sim=$SIM cycles=$MAX_CYCLES jobs=$JOBS out=$O"
echo "wl=$WL"
run_one() {
  w=$1
  if [ ! -f "$B/$w.riscv" ]; then
    echo "SKIP $w (no $B/$w.riscv)"
    return 0
  fi
  echo "START $w $(date +%T)"
  "$SIM" -V -m "$MAX_CYCLES" -v "$O/$w.vcd" "$B/$w.riscv" > "$O/$w.sim.log" 2> >(gzip -1 > "$O/$w.trace.gz") || true
  if [ ! -s "$O/$w.vcd" ]; then echo "FAIL $w empty vcd"; return 0; fi
  echo "DONE $w $(date +%T) vcd=$(wc -c < "$O/$w.vcd") bytes"
}
export -f run_one
export SIM B O MAX_CYCLES
printf '%s\n' $WL | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
echo SIMS_DONE
ls -lh "$O"/*.vcd 2>/dev/null || true
