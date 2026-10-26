"""Scheduling knobs alone: every replay policy on all 10 full traces with a FIXED passive network (0.6 mm^2 decap)."""
import json, numpy as np, pdn3, ga3
from multiprocessing import Pool
WL = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
def job(a):
    w, v = a; z = np.load(ga3.B + '%s_rep_%s.npz' % (w, v)); s = ga3.STARTS[w]
    m = pdn3.simulate(z['I'][s:].astype(np.float32), dict(A_x=[0.6 * x for x in pdn3.SHARE]))
    return w, v, m['p99'], m['p999'], m['max']
if __name__ == '__main__':
    with Pool(4) as P: res = P.map(job, [(w, v) for w in WL for v in ga3.VARS])
    R = {(w, v): r for w, v, *r in res}
    print('droop p99.9 (mV) with 0.6 mm2 passive decap; change vs A in brackets; last col = worst-bench max droop')
    print('%-6s' % 'policy' + ''.join('%11s' % w[:9] for w in WL) + '%9s %9s %8s %8s' % ('mean d', 'worst', 'slow%', 'dummy%'))
    S = ga3.SUMM
    for v in ga3.VARS:
        row = [R[(w, v)][1] for w in WL]; base = [R[(w, 'A')][1] for w in WL]
        print('%-6s' % v + ''.join('%6.1f(%+4.1f)' % (x, x - b) for x, b in zip(row, base)) + '%9.2f %9.1f %8.2f %8.1f' % (
            np.mean(np.array(row) - np.array(base)), max(R[(w, v)][2] for w in WL),
            np.mean([S[w]['slow'][v] for w in WL]), np.mean([S[w]['dummy'][v] for w in WL])))
    json.dump([list(r) for r in res], open('policy_scan.json', 'w'))
