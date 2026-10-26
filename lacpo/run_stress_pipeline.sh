#!/bin/bash
# stressmark: VCD -> LACPo features/trees, cache events, issue trace -> tile currents + replay inputs (existing kernels)
R=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations; O=$R/traces/lacpo_cy100/bench; P=$R/.envs/lacpo37/bin/python; w=${1:-stress}
cd /home/ahilnandan/Documents/PMIC/pdn_work/lacpo
( $P fast_vcd2features.py --vcd $O/$w.vcd --out $O/${w}_features.csv --benchmark $w > $O/$w.feat.log 2>&1 && cd $R && $P tools/run_lacpo_flow.py --features-csv $O/${w}_features.csv --out $O/${w}_out --name $w > $O/$w.lacpo.log 2>&1; echo "lacpo exit $? $(grep 'composed P' $O/$w.lacpo.log)" ) &
( $P cache_events.py --vcd $O/$w.vcd --out $O/${w}_cache.npz > $O/$w.cache.log 2>&1; echo "cache exit $?" ) &
( $P issue_trace.py --vcd $O/$w.vcd --out $O/${w}_issue.npz > $O/$w.issue.log 2>&1; echo "issue exit $?" ) &
wait
OPENBLAS_NUM_THREADS=4 $P stress_build.py $w 2>&1 | grep -v Warn
echo STRESS_PIPELINE_DONE
