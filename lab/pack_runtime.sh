#!/usr/bin/env bash
# Pack the already-built simulator and the Python envs for a native Ubuntu machine.
# Does not pack the chipyard source tree. Run inside the WSL distro on the laptop.
# Output: /mnt/d/WSL/segsee-runtime.tar.gz when D: is mounted, else ./segsee-runtime.tar.gz
set -euo pipefail
ROOT=/opt/segsee
SIM=chipyard-cy100/sims/verilator/simulator-example-MediumBoomConfig-debug
LIB=chipyard-cy100/riscv-tools-install/lib
BENCH=chipyard-cy100/riscv-tools-install/riscv64-unknown-elf/share/riscv-tests/benchmarks
OUT=${1:-}
if [ -z "$OUT" ]; then
  if [ -d /mnt/d/WSL ]; then OUT=/mnt/d/WSL/segsee-runtime.tar.gz
  else OUT=$PWD/segsee-runtime.tar.gz
  fi
fi
test -x "$ROOT/$SIM"
echo "packing -> $OUT"
tar -C "$ROOT" -czf "$OUT" \
  "$SIM" \
  "$LIB" \
  "$BENCH" \
  lacpo37 \
  venv
ls -lh "$OUT"
echo PACK_DONE
