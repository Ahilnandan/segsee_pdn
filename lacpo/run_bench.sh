#!/usr/bin/env bash
# Sequential: MediumBoomV3Config debug sim (VCD + BACT issue trace) -> scope alias -> fast features -> LACPo trees.
set -u
ROOT=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations
SIM=/home/ahilnandan/Documents/PMIC/chipyard/sims/verilator/simulator-chipyard.harness-MediumBoomV3Config-debug
BEN=/home/ahilnandan/Documents/PMIC/chipyard/.conda-env/riscv-tools/riscv64-unknown-elf/share/riscv-tests/benchmarks
OUT=$ROOT/traces/lacpo_local/bench
PY=$ROOT/.envs/lacpo37/bin/python
for wl in "$@"; do
  echo "=== $wl $(date -Is)"
  cd $OUT
  $SIM +permissive +verbose +vcdfile=$wl.vcd +permissive-off $BEN/$wl.riscv > $wl.sim.log 2>&1; echo "sim exit $?"
  ls -lh $wl.vcd
  python3 /home/ahilnandan/Documents/PMIC/pdn_work/vcd_alias.py $wl.vcd ${wl}_a.vcd && rm -f $wl.vcd
  $PY /home/ahilnandan/Documents/PMIC/pdn_work/lacpo/fast_vcd2features.py --vcd ${wl}_a.vcd --out ${wl}_features.csv --benchmark $wl 2>&1 | grep -v "Deprecation\|LARGE_SPARSE"
  mkdir -p $OUT/$wl && cd $ROOT && $PY tools/run_lacpo_flow.py --features-csv $OUT/${wl}_features.csv --name $wl --out $OUT/$wl 2>&1 | grep -v "Deprecation\|LARGE_SPARSE" | tail -4
  gzip -f $OUT/$wl.sim.log
  echo "=== $wl done $(date -Is)"
done
echo ALL_DONE
