"""Per-tile droop map, passive designs, new flip-chip package (pdn5 defaults), full execution traces."""
import os, sys, numpy as np
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pdn5
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
STARTS = dict(qsort=44786, rsort=60805, mm=33082, dhrystone=39788, spmv=87078, multiply=15819, memcpy=43571, median=19595, vvadd=20319, towers=15683)
def job(a):
    w, Ax = a; I = np.load(B + w + '_tileI.npy')[STARTS[w]:].astype(np.float32)
    s = pdn5.Sim(dict(A_x=Ax * pdn5.CELL_SHARE), Im_A=I[:3000].mean(0) * 1e-3); vmin, vavg, _ = s.run(I)
    d = (pdn5.VNOM - vmin) * 1e3
    return w, Ax, np.percentile(d, 99.9, axis=0), d.max(0), np.percentile(d.max(1), 99.9), d.max()
if __name__ == '__main__':
    import multiprocessing as mp
    WL = list(STARTS); jobs = [(w, a) for a in (0.0, 0.5, 0.9) for w in WL]
    with mp.get_context('spawn').Pool(10) as P: res = P.map(job, jobs)
    for Ax in (0.0, 0.5, 0.9):
        R = [r for r in res if r[1] == Ax]
        p999 = np.max([r[2] for r in R], 0).reshape(8, 8); mx = np.max([r[3] for r in R], 0).reshape(8, 8)
        print('\npassive decap %.1f mm2 (by cell share): worst-tile p99.9 per bench %s' % (Ax, ' '.join('%s:%.0f/%.0f' % (r[0][:5], r[4], r[5]) for r in R)))
        print('  per-tile p99.9 droop, worst over benches (mV):')
        for row in p999: print('   ' + ' '.join('%5.1f' % x for x in row))
        print('  per-tile max droop (mV): min %.1f max %.1f' % (mx.min(), mx.max()))
    np.save('tile_map_passive.npy', np.array([[r[2], r[3]] for r in res]))
