"""Phase 4.1: full-trace closed-loop validation of GA5 Pareto designs + all-decap baseline, all 10 benches.
Usage: python val5.py <budget> [workers]. Reads ga5_<budget>.json, writes val5_<budget>.json.
Each job is one (design, bench) full trace; a worker caches only the last bench to bound memory."""
import os, sys, json, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
BENCHES = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
NREP = 4
def front(budget):
    d = json.load(open(HERE + 'ga5_%s.json' % budget)); X, F = d['X'], np.array(d['F'], float)
    ok = np.where((F < 1e2).all(1))[0]; Fo = F[ok]
    nd = [ok[i] for i in range(len(ok)) if not any(all(Fo[j] <= Fo[i]) and any(Fo[j] < Fo[i]) for j in range(len(ok)))]
    nd = np.array(nd); Fn = F[nd]; sc = (Fn - Fn.min(0)) / (np.ptp(Fn, 0) + 1e-12); knee = nd[np.argmin(sc.sum(1))]
    picks = [nd[np.argmin(Fn[:, 0])], nd[np.argmin(Fn[:, 1])], nd[np.argmin(Fn[:, 2])], knee]
    names = ['min_p999', 'min_max', 'min_energy', 'knee']; out = []
    seen = set()
    for nm, i in zip(names, picks):
        if int(i) in seen: continue
        seen.add(int(i)); out.append((nm, [float(v) for v in X[i]], F[i].tolist()))
    return out[:NREP], len(nd)
def baseline(): return ('all_decap', [0.0] * 16, None)
_cache = {}
def job(a):
    budget, name, g, bench = a
    import closedloop as CLM, ga5
    if bench not in _cache: _cache.clear(); _cache[bench] = CLM.CL(bench)
    cl = _cache[bench]; D, pol, prm, info = ga5.decode(g, float(budget))
    r = cl.run(D, pol, prm, a=0, b=cl.n)
    return dict(budget=budget, design=name, bench=bench, pol=pol, n=cl.n, p999=r['p999'], max=r['max'], p99=r['p99'],
                slowdown=r['slowdown'], energy_over=r['energy_over'], nhold=r['nhold'], nreo=r['nreo'], info=info)
if __name__ == '__main__':
    import multiprocessing as mp
    budget = sys.argv[1]; workers = min(6, int(sys.argv[2]) if len(sys.argv) > 2 else 6)
    fr, nnd = front(budget); designs = [baseline()] + fr
    print('budget %s: %d nondominated, validating %d designs x %d benches' % (budget, nnd, len(designs), len(BENCHES)), flush=True)
    jobs = [(budget, nm, g, b) for b in BENCHES for (nm, g, _) in designs]
    with mp.get_context('spawn').Pool(workers) as P: res = P.map(job, jobs, chunksize=1)
    json.dump(dict(budget=budget, designs=[(nm, g, ga) for nm, g, ga in designs], results=res), open(HERE + 'val5_%s.json' % budget, 'w'), indent=1, default=float)
    print('%-12s %-10s %7s %7s %7s %8s %9s' % ('design', 'bench', 'p99.9', 'max', 'p99', 'slow%', 'energy%'))
    for r in res: print('%-12s %-10s %7.1f %7.1f %7.1f %8.2f %9.2f' % (r['design'], r['bench'], r['p999'], r['max'], r['p99'], 100 * r['slowdown'], 100 * r['energy_over']))
    print('VAL_DONE', budget, flush=True)
