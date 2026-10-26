#!/bin/bash
# LACPo training core (chipyard 717c4658 / BOOM 4e9d496d MediumBoomConfig): all benches, sims in parallel,
# features 3-at-a-time (RAM), then LACPo trees
R=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations
O=$R/traces/lacpo_cy100/bench; cd $O
export LD_LIBRARY_PATH=/home/ahilnandan/Documents/PMIC/cy100-deps/riscv-install/lib
S=/home/ahilnandan/Documents/PMIC/chipyard-1.0.0/sims/verilator/simulator-example-MediumBoomConfig-debug.boom4e9d496d
B=/home/ahilnandan/Documents/PMIC/chipyard/.conda-env/riscv-tools/riscv64-unknown-elf/share/riscv-tests/benchmarks
WL="${WL:-qsort rsort mm dhrystone spmv multiply memcpy median vvadd towers}"
for w in $WL; do
  ( $S -V -m 5000000 -v $w.vcd $B/$w.riscv 2> >(gzip -1 > $w.trace.gz) > $w.sim.log; echo "sim $w exit $?" ) &
done; wait
echo SIMS_DONE
feat() {
  w=$1; cd /home/ahilnandan/Documents/PMIC/pdn_work/lacpo
  $R/.envs/lacpo37/bin/python fast_vcd2features.py --vcd $O/$w.vcd --out $O/${w}_features.csv --benchmark $w > $O/$w.feat.log 2>&1 || { echo "feat $w FAILED"; return; }
  cd $R && .envs/lacpo37/bin/python tools/run_lacpo_flow.py --features-csv $O/${w}_features.csv --out $O/${w}_out --name $w > $O/$w.lacpo.log 2>&1 \
    && echo "lacpo $w $(grep 'composed P' $O/$w.lacpo.log)" || echo "lacpo $w FAILED"
}
export -f feat; export R O
echo $WL | tr ' ' '\n' | xargs -P3 -I{} bash -c 'feat {}'
echo ALL_DONE
