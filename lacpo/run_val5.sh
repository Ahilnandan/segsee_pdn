#!/bin/bash
cd /home/ahilnandan/Documents/PMIC/pdn_work/lacpo
export OPENBLAS_NUM_THREADS=1
P=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/.envs/lacpo37/bin/python
until grep -q ALL_DONE ga5_runner.log 2>/dev/null; do sleep 60; done
for b in 0.35 0.2 0.5; do
  $P val5.py $b 6 > val5_$b.log 2>&1
  echo "val5 $b exit $?" >> val5_runner.log
done
echo VAL_ALL_DONE >> val5_runner.log
