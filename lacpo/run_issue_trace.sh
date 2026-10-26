#!/bin/bash
P=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/.envs/lacpo37/bin/python
B=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench
cd /home/ahilnandan/Documents/PMIC/pdn_work/lacpo
for w in qsort rsort mm dhrystone spmv multiply memcpy median vvadd towers; do
  ( $P issue_trace.py --vcd $B/$w.vcd --out $B/${w}_issue.npz 2>&1 | grep -v Warn | sed "s/^/$w /" ) &
done; wait; echo ISSUE_DONE
