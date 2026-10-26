"""Multi-core open-loop test (pdn7, 2x2 cores): how burst alignment across cores sets droop on a shared grid."""
import os, sys, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn6, pdn7, eval6
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'; W = 60000
def win(w, off=0):
    I = np.load(B + w + "_tileI.npy", mmap_mode="r"); s = 70000 + off; return np.asarray(I[s:s + W]).astype(float)
def run(mod, d, I):
    sim = mod.Sim(d, Im_A=I[:3000].mean(0) * 1e-3); vmin, _, _ = sim.run(I.astype(np.float32), vref=np.full(len(I), mod.VNOM))
    return (mod.VNOM - vmin) * 1e3
D, _ = eval6.designs(); rng = np.random.RandomState(3)
for dn in ('decap', 'banks'):
    d1 = D[dn]; d4 = pdn7.place_design(d1); print('== design', dn, '(per core)')
    x = run(pdn6, d1, win('qsort'))[2000:].max(1); print('  1 core qsort                    p99.9 %5.1f  max %5.1f' % (np.percentile(x, 99.9), x.max()))
    cases = [('4x qsort, aligned (lockstep)', [win('qsort')] * 4), ('4x qsort, random offsets', [win('qsort', int(o)) for o in rng.randint(0, 20000, 4)]),
             ('4 different (qsort rsort mm dhry)', [win(w) for w in ('qsort', 'rsort', 'mm', 'dhrystone')])]
    for nm, arrs in cases:
        dr = run(pdn7, d4, pdn7.place(arrs))[2000:]; wt = dr.max(1); pc = [dr[:, pdn7.tiles_of(k)].max(1) for k in range(4)]
        tot = sum(a.sum(1) for a in arrs); dI = np.abs(np.diff(tot))
        print('  %-31s p99.9 %5.1f  max %5.1f | per-core p99.9 %s | total-current step p99.9 %.0f mA' % (nm, np.percentile(wt, 99.9), wt.max(), ' '.join('%.1f' % np.percentile(p, 99.9) for p in pc), np.percentile(dI, 99.9)))
