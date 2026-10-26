#!/bin/bash
P=/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/.envs/lacpo37/bin/python
S=/tmp/claude-1004/-home-ahilnandan-Documents-PMIC/4a0dca89-208b-401b-940b-844c3bbea50c/scratchpad
cd /home/ahilnandan/Documents/PMIC/pdn_work/lacpo; rm -f $S/rep_*.log
for w in qsort rsort mm dhrystone spmv multiply memcpy median vvadd towers; do
  for v in A,R,H4,H6,H10,D3 D5,D8,RH6,HD6,RHD6,RHD4; do
    ( $P replay3.py $w $v 2>&1 | grep -v Warn > $S/rep_${w}_${v//,/_}.log ) &
  done
done; wait
cat $S/rep_*.log; echo REPLAY_DONE
