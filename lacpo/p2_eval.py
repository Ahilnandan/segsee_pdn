"""Phase 2 first evaluation: baseline vs token protocol (kernel / oracle predictor) +- charge-pump LDO banks +- load line,
on 20k-cycle windows around each bench's worst droop event (closed loop, tiled PDN)."""
import os, sys, json, numpy as np
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
EV = json.load(open(os.path.dirname(os.path.abspath(__file__)) + '/windows4.json'))
BENCH = ['rsort', 'memcpy', 'qsort', 'dhrystone', 'mm']
BANK_TILES = [r * 8 + c for c in (3, 4, 5, 6) for r in range(8) if (r % 2 == 0)]       # 16 tiles over core + D$ columns
def designs():
    import pdn5
    base = dict(A_x=0.3 * pdn5.CELL_SHARE, R_d=np.full(64, 0.3))
    bank = {t: dict(C_res=0.6e-9 / len(BANK_TILES), R_pass=0.5 * len(BANK_TILES) / 4, gm=12.0 / len(BANK_TILES), tau=1.2e-9, off=0.02,
                    I_ch=0.02 / len(BANK_TILES)) for t in BANK_TILES}
    return {'decap0.3': base, 'decap0.3+bank0.6nF': dict(base, bank=bank)}
POL = {'A': ('A', {}), 'TOK': ('TOK', {}), 'TOK-reo-only': ('TOK', dict(hold=0)), 'TOK-oracle': ('TOK', dict(oracle=1)),
       'TOK-oracle-reo': ('TOK', dict(oracle=1, hold=0)), 'TOK+LL': ('TOK', dict(R_LL=0.15, vg=0.0))}
def job(a):
    w, dn, pn = a
    import closedloop as CLM
    cl = CLM.CL(w); D = designs()[dn]; pol, prm = POL[pn]
    k = EV[w]['A']['win'][0]; c = (k[0] + k[1]) // 2; lo = max(0, c - 10000); hi = min(cl.n, c + 10000)
    r = cl.run(D, pol, prm, a=lo, b=hi)
    return w, dn, pn, r
if __name__ == '__main__':
    import multiprocessing as mp
    jobs = [(w, d, p) for w in BENCH for d in ('decap0.3', 'decap0.3+bank0.6nF') for p in POL]
    with mp.get_context('spawn').Pool(10) as P: res = P.map(job, jobs)
    print('%-10s %-20s %-15s %7s %7s %8s %8s %7s %6s %6s' % ('bench', 'design', 'policy', 'p99.9', 'max', 'slow%', 'energy%', 'holds', 'reo', 'arm%'))
    for w, dn, pn, r in res:
        print('%-10s %-20s %-15s %7.1f %7.1f %8.2f %8.2f %7d %6d %6.1f' % (w, dn, pn, r['p999'], r['max'], 100 * r['slowdown'], 100 * r['energy_over'], r['nhold'], r['nreo'], 100 * r['arm_frac']))
    json.dump([(w, d, p, {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}) for w, d, p, r in res], open('p2_eval.json', 'w'), indent=1, default=float)
