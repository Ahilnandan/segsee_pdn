#!/usr/bin/env bash
# Run riscv-tests on the lab machine's own chipyard MediumBoom simulator.
# One benchmark per core. This binary's verilator --threads setting is whatever that chipyard built.
#
#   CHIPYARD_ROOT=~/chipyard MAX_CYCLES=5000000 JOBS=$(nproc) bash lab/run_benches_parallel.sh
#
# Resolution order for the simulator:
#   $SIM if it is executable, else
#   $CHIPYARD_ROOT or $CHIPYARD_DIR, else
#   ~/chipyard, ~/chipyard-1.0.0, /opt/chipyard
# Then sims/verilator/simulator-*MediumBoom*debug* (newest).
# Benchmarks come from $RISCV, else from that chipyard's riscv toolchain.
# Does not look at /opt/segsee unless CHIPYARD_ROOT points there.
#
# Writes $SEGSEE_TRACES/<w>.vcd, <w>.sim.log, <w>.trace.gz
# "*** FAILED *** via trace_count" is the cycle cap, not a crashed core.
set -euo pipefail
MAX_CYCLES=${MAX_CYCLES:-5000000}
JOBS=${JOBS:-$(nproc)}
WL="${WL:-qsort rsort mm dhrystone spmv multiply memcpy median vvadd towers}"
O=${SEGSEE_TRACES:-$HOME/segsee_traces/lacpo_cy100/bench}

find_chipyard() {
  if [ -n "${CHIPYARD_ROOT:-${CHIPYARD_DIR:-}}" ]; then
    echo "${CHIPYARD_ROOT:-$CHIPYARD_DIR}"
    return
  fi
  local c
  for c in "$HOME/chipyard" "$HOME/chipyard-1.0.0" /opt/chipyard; do
    if [ -d "$c/sims/verilator" ]; then echo "$c"; return; fi
  done
  return 1
}

CY=$(find_chipyard) || {
  echo "No chipyard tree. Set CHIPYARD_ROOT to the lab machine's checkout." >&2
  exit 1
}

if [ -n "${SIM:-}" ] && [ -x "$SIM" ]; then
  :
else
  SIM=$(ls -1t "$CY"/sims/verilator/simulator-*MediumBoom*debug* 2>/dev/null | head -1 || true)
  if [ -z "$SIM" ]; then
    SIM=$(ls -1t "$CY"/sims/verilator/simulator-*MediumBoom* 2>/dev/null | head -1 || true)
  fi
fi
if [ -z "${SIM:-}" ] || [ ! -x "$SIM" ]; then
  echo "No MediumBoom simulator under $CY/sims/verilator." >&2
  echo "Build it with that chipyard: cd $CY/sims/verilator && make CONFIG=MediumBoomConfig debug -j\$(nproc)" >&2
  exit 1
fi

if [ -z "${RISCV:-}" ]; then
  if [ -d "$CY/riscv-tools-install" ]; then RISCV=$CY/riscv-tools-install
  elif [ -d "$CY/.conda-env/riscv-tools" ]; then RISCV=$CY/.conda-env/riscv-tools
  else
    echo "Set RISCV to this chipyard's riscv toolchain prefix." >&2
    exit 1
  fi
fi
B=${BENCH_DIR:-$RISCV/riscv64-unknown-elf/share/riscv-tests/benchmarks}
export LD_LIBRARY_PATH=$RISCV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export RISCV

mkdir -p "$O"
cd "$O"
echo "chipyard=$CY"
echo "sim=$SIM"
if [ -d "$CY/generators/boom/.git" ]; then
  echo "boom=$(git -C "$CY/generators/boom" rev-parse --short=12 HEAD) $(git -C "$CY/generators/boom" log -1 --format=%s)"
fi
echo "riscv=$RISCV"
echo "benchmarks=$B"
echo "cycles=$MAX_CYCLES jobs=$JOBS out=$O"
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
