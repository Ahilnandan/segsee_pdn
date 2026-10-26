"""Phase 3: NSGA-II at a FIXED total capacitor budget (decap area + bank reservoir area), closed loop, tiled PDN.
Objectives (min): worst p99.9 droop, worst max droop (over 5 benches' worst-event windows), energy overhead (V^2 + device
losses vs nominal), slowdown. Genes: decap/bank split, bank count/placement, bank params, damping, token-protocol params."""
import os, sys, json, time, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
EV = json.load(open(HERE + 'windows4.json'))
BENCH = ['rsort', 'memcpy', 'qsort', 'dhrystone', 'mm']
NG = 16; WIN = 20000
_cl = {}
def cl_for(w):
    if w not in _cl:
        import closedloop as CLM; _cl[w] = CLM.CL(w)
    return _cl[w]
def decode(g, budget):
    import pdn5
    g = np.clip(np.asarray(g, float), 0, 1)
    bank_frac = 0.6 * g[0]                                       # share of the capacitor budget given to bank reservoirs
    nbank = [0, 2, 4, 8, 16][min(4, int(g[1] * 5))]
    if nbank == 0: bank_frac = 0.0
    A_dec = budget * (1 - bank_frac); C_res_tot = budget * bank_frac * pdn5.CP_DENS
    cols = [3, 4, 5, 6] if g[2] < 0.5 else list(range(8))          # banks over core/D$ columns or whole die
    cand = [r * 8 + c for c in cols for r in range(8)]; step = max(1, len(cand) // max(1, nbank)); tiles = cand[::step][:nbank]
    D = dict(A_x=A_dec * pdn5.CELL_SHARE, R_d=np.full(64, 0.6 * g[3]))
    if nbank:
        D['bank'] = {t: dict(C_res=C_res_tot / nbank, R_pass=0.2 + 1.8 * g[4], gm=(2 + 30 * g[5]) / nbank, tau=0.5e-9 * 10 ** g[6],
                             off=0.008 + 0.03 * g[7], I_ch=(0.004 + 0.04 * g[8]) / nbank) for t in tiles}
    prm = dict(th=4 + 20 * g[9], kv=0.6 * g[10], kb=8 * g[11], maxd=1 + int(5 * g[12]), hold=int(g[13] > 0.5), reo=int(g[14] > 0.3),
               th_arm=4 + 20 * g[15])
    pol = 'TOK' if (prm['hold'] or prm['reo'] or nbank) else 'A'
    return D, pol, prm, dict(nbank=nbank, bank_frac=bank_frac, cols=cols)
def evaluate(args):
    g, budget = args
    try:
        D, pol, prm, info = decode(g, budget); p999 = mx = en = sl = 0.0
        for w in BENCH:
            cl = cl_for(w); k = EV[w]['A']['win'][0]; c = (k[0] + k[1]) // 2
            r = cl.run(D, pol, prm, a=max(0, c - WIN // 2), b=min(cl.n, c + WIN // 2))
            p999 = max(p999, r['p999']); mx = max(mx, r['max']); en += r['energy_over'] / len(BENCH); sl += max(0, r['slowdown']) / len(BENCH)
        return [p999, mx, en * 100, sl * 100]
    except Exception as e:
        return [1e3, 1e3, 1e3, 1e3]
def dom(a, b): return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))
def nsga(F):
    n = len(F); S = [[] for _ in range(n)]; cnt = [0] * n; fronts = [[]]; rank = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j:
                if dom(F[i], F[j]): S[i].append(j)
                elif dom(F[j], F[i]): cnt[i] += 1
        if cnt[i] == 0: fronts[0].append(i)
    k = 0
    while fronts[k]:
        nx = []
        for i in fronts[k]:
            for j in S[i]:
                cnt[j] -= 1
                if cnt[j] == 0: rank[j] = k + 1; nx.append(j)
        k += 1; fronts.append(nx)
    crowd = [0.0] * n
    for fr in fronts:
        if not fr: continue
        for m in range(len(F[0])):
            f2 = sorted(fr, key=lambda i: F[i][m]); lo, hi = F[f2[0]][m], F[f2[-1]][m]; crowd[f2[0]] = crowd[f2[-1]] = 1e9
            for a in range(1, len(f2) - 1): crowd[f2[a]] += (F[f2[a + 1]][m] - F[f2[a - 1]][m]) / (hi - lo + 1e-12)
    return rank, crowd
def main(budget, pop=36, gens=30, workers=4, seed=0):
    import multiprocessing as mp
    rng = np.random.RandomState(seed); X = rng.rand(pop, NG); X[0] = 0; X[0, 1] = 0.0          # seed: all decap, baseline policy
    pool = mp.get_context('spawn').Pool(workers); F = pool.map(evaluate, [(x, budget) for x in X]); out = HERE + 'ga5_%s.json' % budget
    for gen in range(gens):
        t0 = time.time(); rank, crowd = nsga(F)
        def tour():
            i, j = rng.randint(pop, size=2); return i if (rank[i], -crowd[i]) < (rank[j], -crowd[j]) else j
        kids = []
        while len(kids) < pop:
            a, b = X[tour()], X[tour()]; u = rng.rand(NG); beta = np.where(u < 0.5, (2 * u) ** (1 / 16.), (1 / (2 * (1 - u))) ** (1 / 16.))
            for c in (0.5 * ((1 + beta) * a + (1 - beta) * b), 0.5 * ((1 - beta) * a + (1 + beta) * b)):
                mu = rng.rand(NG) < 2.0 / NG; c[mu] += rng.normal(0, 0.15, mu.sum()); kids.append(np.clip(c, 0, 1))
        kids = np.array(kids[:pop]); Fk = pool.map(evaluate, [(x, budget) for x in kids])
        XA = np.vstack([X, kids]); FA = F + Fk; rank, crowd = nsga(FA); order = sorted(range(len(FA)), key=lambda i: (rank[i], -crowd[i]))[:pop]
        X = XA[order]; F = [FA[i] for i in order]
        best = sorted(F)[:4]
        print('[%.2f] gen %2d front0 %d  best p99.9 %s  (%.0fs)' % (budget, gen, sum(1 for i in range(pop) if rank[order[0]] == 0 and dom is not None),
              ' '.join('(%.1f,%.1f,%.2f%%,%.2f%%)' % tuple(f) for f in best), time.time() - t0), flush=True)
        json.dump(dict(budget=budget, X=X.tolist(), F=F), open(out, 'w'))
    pool.close()
if __name__ == '__main__':
    main(float(sys.argv[1]), pop=int(sys.argv[2]) if len(sys.argv) > 2 else 36, gens=int(sys.argv[3]) if len(sys.argv) > 3 else 30)
