"""NSGA-II co-design of the on-die PDN + issue policy for the MediumBoom tile (LACPo + cache-model currents).

Objectives (minimise): f1 = on-die PDN area (mm^2: explicit decap + active-decap reservoir + on-die inductors/switches)
                       f2 = overhead = slowdown + (dummy charge + injector loss + active-decap power) / load energy
Constraints (every fit bench): droop p99.9 <= 55 mV and max <= 80 mV (constraint-domination).
Fit set: droop windows of 5 benches; validation of the final front on all 10 full traces (ga3_validate.py)."""
from __future__ import print_function
import json, os, sys, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn3
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
FIT = ['rsort', 'memcpy', 'dhrystone', 'qsort', 'mm']; WIN = 60000
VARS = ['A', 'R', 'H4', 'H6', 'H10', 'D3', 'D5', 'D8', 'RH6', 'HD6', 'RHD6', 'RHD4']
P999, PMAX = 55.0, 80.0
NG = 25
SUMM = json.load(open(HERE + 'replay_summary.json'))
STARTS = {r[0]: r[2] for r in json.load(open(HERE + 'passive_scan.json'))}
_cache = {}

def window(w, v):
    key = (w, v)
    if key not in _cache:
        z = np.load(B + '%s_rep_%s.npz' % (w, v)); s = STARTS[w]
        _cache[key] = (z['I'][s:s + WIN].astype(np.float32), z['Ihat'][s:s + WIN].astype(np.float32))
    return _cache[key]

def decode(g):
    g = np.clip(np.asarray(g, float), 0, 1)
    d = dict(A_x=[0.5 * x for x in g[0:4]], R_d=[1.0 * x for x in g[4:8]],
             Ca=[max(0.0, (x - 0.3) / 0.7) * 2e-9 for x in g[8:12]], G=[2, 4, 8][min(2, int(g[12] * 3))],
             amp=0.5 + 1.5 * g[13], inj=g[14] > 0.5, ind_pkg=g[15] > 0.5, L=[5e-9, 10e-9, 20e-9][min(2, int(g[16] * 3))],
             Vin=[1.3, 1.5][min(1, int(g[17] * 2))], Ron=[0.5, 1.0, 2.0][min(2, int(g[18] * 3))], Ipk=0.05 + 0.25 * g[19],
             th_b=0.010 + 0.040 * g[20], pred=g[21] > 0.5, th_p=2 + 13 * g[22], lead=1 + int(g[23] * 5.999),
             var=VARS[min(len(VARS) - 1, int(g[24] * len(VARS)))])
    return d

def spec_of(d):
    spec = []
    for r in range(4):
        s = {}
        if d['Ca'][r] > 0: s['ad'] = dict(Ca=d['Ca'][r], G=d['G'], tau_a=1e-9 / d['amp'], Pq=1.8e-3 * d['G'] / 4 * d['amp'])
        if d['inj']: s['buck'] = dict(L=d['L'], Vin=d['Vin'], Ron=d['Ron'], Ipk=d['Ipk'] * pdn3.SHARE[r], th=d['th_b'], pred=int(d['pred']), dwell=2)
        spec.append(s or None)
    return spec

def area_of(d):
    a = sum(d['A_x']) + sum(d['Ca']) / pdn3.HV_DENS
    if d['inj']: a += 4 * (0.005 if d['ind_pkg'] else pdn3.IND_MM2_PER_NH * d['L'] / 1e-9)
    return a

def evaluate(g, benches=FIT, full=False):
    d = decode(g); spec = spec_of(d); p = dict(A_x=d['A_x'], R_d=d['R_d'])
    cv = 0.0; ov = []; per = {}
    for w in benches:
        if full:
            z = np.load(B + '%s_rep_%s.npz' % (w, d['var'])); s = STARTS[w]; I = z['I'][s:].astype(np.float32); Ih = z['Ihat'][s:]
        else:
            I, Ih = window(w, d['var'])
        pred = pdn3.predict_triggers(Ih, d['th_p'], d['lead']) if (d['inj'] and d['pred']) else None
        m = pdn3.simulate(I, p, spec, pred=pred)
        sl = SUMM[w]['slow'][d['var']] / 100.0; dm = SUMM[w]['dummy'][d['var']] / 100.0
        o = max(0.0, sl) + dm + m['loss_buck'] + m['loss_ad']
        ov.append(o); cv += max(0, m['p999'] - P999) / P999 + max(0, m['max'] - PMAX) / PMAX
        per[w] = dict(p999=m['p999'], max=m['max'], p99=m['p99'], over=m['over_max'], slow=sl, dummy=dm, buck=m['loss_buck'], ad=m['loss_ad'])
    return area_of(d), float(np.mean(ov)), cv, per

def _ev(g):
    try:
        a, o, cv, per = evaluate(g); return a, o, cv
    except Exception as e:
        return 1e9, 1e9, 1e9

def dominates(a, b):
    if a[2] < b[2] - 1e-12: return True
    if b[2] < a[2] - 1e-12: return False
    return (a[0] <= b[0] and a[1] <= b[1]) and (a[0] < b[0] or a[1] < b[1])

def nsga_sort(F):
    n = len(F); S = [[] for _ in range(n)]; cnt = [0] * n; fronts = [[]]; rank = [0] * n
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if dominates(F[i], F[j]): S[i].append(j)
            elif dominates(F[j], F[i]): cnt[i] += 1
        if cnt[i] == 0: fronts[0].append(i)
    k = 0
    while fronts[k]:
        nxt = []
        for i in fronts[k]:
            for j in S[i]:
                cnt[j] -= 1
                if cnt[j] == 0: rank[j] = k + 1; nxt.append(j)
        k += 1; fronts.append(nxt)
    crowd = [0.0] * n
    for fr in fronts:
        if not fr: continue
        for m in (0, 1):
            fr2 = sorted(fr, key=lambda i: F[i][m]); lo, hi = F[fr2[0]][m], F[fr2[-1]][m]
            crowd[fr2[0]] = crowd[fr2[-1]] = 1e9
            for a in range(1, len(fr2) - 1):
                crowd[fr2[a]] += (F[fr2[a + 1]][m] - F[fr2[a - 1]][m]) / (hi - lo + 1e-12)
    return rank, crowd

def main(pop=40, gens=40, seed=1, out='ga3_run.json'):
    rng = np.random.RandomState(seed)
    X = rng.rand(pop, NG)
    # seeds: passive-only designs (no active decap, no injector, policy A) at several decap sizes
    for i, a in enumerate((0.4, 0.6, 0.8, 1.0)):
        X[i] = 0.0; X[i, 0:4] = a / 2.0; X[i, 4:8] = 0.0; X[i, 8:12] = 0.0; X[i, 14] = 0.0; X[i, 24] = 0.0
    pool = Pool(20); F = pool.map(_ev, list(X)); hist = []
    for gen in range(gens):
        t0 = time.time(); rank, crowd = nsga_sort(F)
        def tour():
            i, j = rng.randint(pop, size=2)
            return i if (rank[i], -crowd[i]) < (rank[j], -crowd[j]) else j
        kids = []
        while len(kids) < pop:
            a, b = X[tour()], X[tour()]
            eta = 15.0; u = rng.rand(NG); beta = np.where(u < 0.5, (2 * u) ** (1 / (eta + 1)), (1 / (2 * (1 - u))) ** (1 / (eta + 1)))
            c1 = 0.5 * ((1 + beta) * a + (1 - beta) * b); c2 = 0.5 * ((1 - beta) * a + (1 + beta) * b)
            for c in (c1, c2):
                mut = rng.rand(NG) < 1.0 / NG * 2
                c[mut] += rng.normal(0, 0.15, mut.sum()); kids.append(np.clip(c, 0, 1))
        kids = np.array(kids[:pop]); Fk = pool.map(_ev, list(kids))
        XA = np.vstack([X, kids]); FA = F + Fk; rank, crowd = nsga_sort(FA)
        order = sorted(range(len(FA)), key=lambda i: (rank[i], -crowd[i]))[:pop]
        X = XA[order]; F = [FA[i] for i in order]
        feas = [f for f in F if f[2] == 0]
        best_a = min([f[0] for f in feas]) if feas else None
        hist.append(dict(gen=gen, feasible=len(feas), min_area=best_a, sec=time.time() - t0))
        print('gen %2d  feasible %2d  min feasible area %s  front %s  (%.0fs)' % (gen, len(feas), ('%.3f' % best_a) if feas else '-',
              ' '.join('(%.3f,%.3f)' % (f[0], f[1]) for f in sorted(feas)[:8]), time.time() - t0), flush=True)
        json.dump(dict(X=X.tolist(), F=F, hist=hist), open(HERE + out, 'w'))
    pool.close()

if __name__ == '__main__':
    main(pop=int(sys.argv[1]) if len(sys.argv) > 1 else 40, gens=int(sys.argv[2]) if len(sys.argv) > 2 else 40)
