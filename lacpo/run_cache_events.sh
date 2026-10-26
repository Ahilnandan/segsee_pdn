#!/bin/bash
P=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/.envs/lacpo37/bin/python
B=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench
cd /home/ahilnandan/Documents/PMIC/pdn_work/lacpo
for w in qsort rsort mm dhrystone spmv multiply memcpy median vvadd towers; do
  ( $P cache_events.py --vcd $B/$w.vcd --out $B/${w}_cache.npz > $B/$w.cache.log 2>&1; echo "$w $(grep -h 'cycles\|MISSING\|Error' $B/$w.cache.log)" ) &
done; wait; echo CACHE_DONE
