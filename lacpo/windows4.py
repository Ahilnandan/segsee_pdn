"""Worst-droop windows per (bench, scheduling variant): full execution trace through a reference passive network
(0.4 mm^2 decap, pdn4), top-K droop events separated by >= 5000 cycles, window +-W/2 around each. Saves windows4.json."""
import os, sys, json, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn4, ga3
WL = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
K, W = 3, 3000
def job(a):
    w, v = a; z = np.load(ga3.B + '%s_rep_%s.npz' % (w, v)); s = ga3.STARTS[w]; I = z['I'][s:].astype(np.float32)
    m = pdn4.simulate(I, dict(A_x=[0.4 * x for x in pdn4.SHARE]), want_v=True)
    d = (pdn4.VNOM - m['vmin'].min(1)) * 1e3; order = np.argsort(d)[::-1]; picks = []
    for k in order:
        if all(abs(k - q) >= 5000 for q in picks): picks.append(int(k))
        if len(picks) == K: break
    return w, v, [[max(0, s + q - W // 2), s + q + W // 2] for q in picks], [float(d[q]) for q in picks]
if __name__ == '__main__':
    import multiprocessing as mp
    with mp.get_context('spawn').Pool(10) as P: res = P.map(job, [(w, v) for w in WL for v in ga3.VARS])
    out = {}
    for w, v, wins, dd in res: out.setdefault(w, {})[v] = dict(win=wins, droop=dd)
    json.dump(out, open('windows4.json', 'w'), indent=1)
    for w in WL: print(w, 'A', ['%.0f' % x for x in out[w]['A']['droop']], 'H4', ['%.0f' % x for x in out[w]['H4']['droop']])
