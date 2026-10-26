import json, sys
from multiprocessing import Pool
import smt_exp
smt_exp.MIXES.update({
 'smt8-b': [('qsort', 100000), ('dhrystone', 80000), ('mm', 80000), ('mm', 230000),
            ('qsort', 300000), ('dhrystone', 120000), ('mm', 150000), ('qsort', 200000)],
 'smt8-c': [('qsort', 10000), ('qsort', 120000), ('dhrystone', 0), ('dhrystone', 100000),
            ('mm', 10000), ('mm', 90000), ('mm', 170000), ('qsort', 220000)]})
if __name__ == '__main__':
    jobs = [(m, p, 8) for m in ('smt8-b', 'smt8-c') for p in ('A', 'T')]
    with Pool(4) as p, open('smt_exp2.jsonl', 'w') as f:
        for r in p.imap_unordered(smt_exp.job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
