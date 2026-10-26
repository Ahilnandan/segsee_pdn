#!/bin/bash
cd /home/ahilnandan/Documents/PMIC/pdn_work/lacpo
export OPENBLAS_NUM_THREADS=1
P=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/.envs/lacpo37/bin/python
for b in 0.35 0.2 0.5; do
  $P -c "import ga5; ga5.main($b, pop=36, gens=30, workers=5)" > ga5_$b.log 2>&1
  echo "budget $b exit $?" >> ga5_runner.log
done
echo ALL_DONE >> ga5_runner.log
