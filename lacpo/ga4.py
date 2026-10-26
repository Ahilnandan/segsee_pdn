"""NSGA-II per architecture on pdn4 (role-switchable regulators, LDOs, predictive boost, scheduling policies).

arch A: off-chip sets V + on-die transient injector (package inductor) [+ predictor]
arch B: off-chip current helper + on-die multiphase regulator (fast loop, predictive boost)
arch D: series LDO per region (rail = V + headroom), predictive boost
arch E: arch-A network + parallel LDO droop clamp on a separate rail (bump share taken from VDD), predictive boost
Common genes: decap + damping per region, active decap, scheduling variant, predictor threshold/lead.
Objectives: on-die PDN area (mm^2), overhead (conversion/LDO/amp losses + boost energy + dummy charge + slowdown).
Constraints: p99.9 <= 55 mV on 30k-cycle windows of 5 benches; max <= 80 mV on the 3 worst-event windows of all 10."""
from __future__ import print_function
import json, os, sys, time, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn4, ga3
B = ga3.B; HERE = ga3.HERE; VARS = ga3.VARS; SUMM = ga3.SUMM; STARTS = ga3.STARTS
FIT = ['rsort', 'memcpy', 'dhrystone', 'qsort', 'mm']; WIN = 30000
WL = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
EV = json.load(open(HERE + 'windows4.json'))
P999, PMAX = 55.0, 80.0
NC = 18                                    # common genes
NA = {'A': 7, 'B': 6, 'D': 6, 'E': 8, 'F': 15, 'G': 6}
_cache = {}

def trace(w, v):
    if (w, v) not in _cache:
        z = np.load(B + '%s_rep_%s.npz' % (w, v)); _cache[(w, v)] = (z['I'].astype(np.float32), z['Ihat'].astype(np.float32))
    return _cache[(w, v)]

def pick(lst, x): return lst[min(len(lst) - 1, int(x * len(lst)))]

def decode(arch, g):
    g = np.clip(np.asarray(g, float), 0, 1); c = g[:NC]; a = g[NC:]
    d = dict(arch=arch, A_x=[0.5 * x for x in c[0:4]], R_d=[1.0 * x for x in c[4:8]],
             Ca_tot=max(0.0, (c[8] - 0.4) / 0.6) * 2e-9, G=pick([2, 4, 8], c[9]), amp=0.5 + 1.5 * c[10],
             var=pick(VARS, c[11]), th_p=2 + 13 * c[12], lead=1 + int(c[13] * 5.999), pred_on=c[14] > 0.5,
             boost=0.025 * c[15] if c[14] > 0.5 else 0.0)
    if arch == 'A':
        d.update(inj=a[0] > 0.3, L=pick([5e-9, 10e-9, 20e-9, 50e-9], a[1]), Vin=pick([1.3, 1.5, 1.8], a[2]), Ron=0.3 + 1.7 * a[3],
                 Ipk=0.05 + 0.25 * a[4], th_b=0.02 + 0.04 * a[5], inj_pred=a[6] > 0.5)
    if arch == 'B':
        d.update(N=pick([2, 4, 8], a[0]), L=pick([20e-9, 50e-9, 100e-9, 200e-9], a[1]), fsw=20e6 * 10 ** a[2], Ron=0.5 + 3.5 * a[3],
                 fc=2e6 * 15 ** a[4], f_share=5e3 + 45e3 * a[5])
    if arch == 'D':
        d.update(rail_off=0.04 + 0.11 * a[0], R_pass=0.1 + 0.4 * a[1], gm=2 + 18 * a[2], tau=0.5e-9 * 20 ** a[3], A_rail_tot=0.3 * a[4])
    if arch == 'F':                      # full hybrid: parallel LDO (first droop) + on-die buck injector (tail)
        e = a[:8]; j = a[8:]
        d.update(V_rail=pick([1.2, 1.25, 1.3], e[0]), bump_f=0.1 + 0.3 * e[1], R_pass=0.1 + 0.9 * e[2], gm=1 + 19 * e[3],
                 tau=0.5e-9 * 20 ** e[4], off=0.005 + 0.035 * e[5], A_rail_tot=0.2 * e[6], ldo_on=e[7] > 0.2,
                 inj=j[0] > 0.3, L=pick([5e-9, 10e-9, 20e-9, 50e-9], j[1]), Vin=pick([1.3, 1.5, 1.8], j[2]), Ron=0.3 + 1.7 * j[3],
                 Ipk=0.05 + 0.25 * j[4], th_b=0.02 + 0.04 * j[5], inj_pred=j[6] > 0.5)
    if arch == 'G':                      # on-die charge pump (2:1 from Vdd) -> reservoir -> regulated pass device
        d.update(R_pass=0.1 + 0.9 * a[0], gm=1 + 19 * a[1], tau=0.5e-9 * 20 ** a[2], off=0.005 + 0.035 * a[3],
                 A_cp_tot=0.4 * a[4], I_ch_tot=0.001 + 0.039 * a[5])
    if arch == 'E':
        d.update(V_rail=pick([1.2, 1.25, 1.3], a[0]), bump_f=0.1 + 0.3 * a[1], R_pass=0.1 + 0.9 * a[2], gm=1 + 19 * a[3],
                 tau=0.5e-9 * 20 ** a[4], off=0.005 + 0.035 * a[5], A_rail_tot=0.2 * a[6], ldo_on=a[7] > 0.2)
    return d

def build_design(d):
    S = pdn4.SHARE; p = dict(A_x=d['A_x'], R_d=d['R_d']); spec = [dict() for _ in range(4)]
    if d['Ca_tot'] > 0:
        for r in range(4): spec[r]['ad'] = dict(Ca=d['Ca_tot'] * S[r], G=d['G'], amp=d['amp'])
    role_off = None
    if d['arch'] in ('A', 'F') and d['inj']:
        for r in range(4): spec[r]['reg'] = dict(role=1, L=d['L'], on_die=False, Vin=d['Vin'], Ron=d['Ron'], Ipk=d['Ipk'] * S[r], th=d['th_b'], dwell=2)
    if d['arch'] == 'B':
        for r in range(4): spec[r]['reg'] = dict(role=2, L=d['L'], on_die=False, Vin=1.8, Ron=d['Ron'], phases=d['N'], fsw=d['fsw'], fc=d['fc'], boost=d['boost'])
        p['f_share'] = d['f_share']; role_off = 'helper'
    if d['arch'] == 'D':
        p.update(series_ldo=True, rail_off=d['rail_off'], A_rail=[d['A_rail_tot'] * s for s in S])
        for r in range(4): spec[r]['ldo'] = dict(series=1, R_pass=d['R_pass'], gm=d['gm'] * S[r] * 4, tau=d['tau'], off=0.0, boost=d['boost'], Iq=0.5e-3 * (1e-9 / d['tau']))
    if d['arch'] in ('E', 'F') and d['ldo_on']:
        f = d['bump_f']; p.update(ldo_rail=True, L_p=0.5e-9 / (1 - f), R_p=50e-3 / (1 - f), L_p2=0.5e-9 / f, R_p2=50e-3 / f, A_rail=[d['A_rail_tot'] * s for s in S])
        for r in range(4): spec[r]['ldo'] = dict(series=0, V=d['V_rail'], R_pass=d['R_pass'], gm=d['gm'] * S[r] * 4, tau=d['tau'], off=d['off'], boost=d['boost'], Iq=0.5e-3 * (1e-9 / d['tau']))
    if d['arch'] == 'G':
        p.update(cp_rail=True, A_cp=[d['A_cp_tot'] * s for s in S])
        for r in range(4): spec[r]['ldo'] = dict(series=0, R_pass=d['R_pass'], gm=d['gm'] * S[r] * 4, tau=d['tau'], off=d['off'], boost=d['boost'],
                                                 Iq=0.5e-3 * (1e-9 / d['tau']), cp=dict(I_ch=d['I_ch_tot'] * S[r], R_cp=5.0, eta=0.9, Vmul=2.0))
    spec = [s if s else None for s in spec]
    return p, spec, role_off

def area_of(d, p, spec):
    return pdn4.area(dict(pdn4.DEFAULT, **p), spec)['total']

def sim(d, p, spec, role_off, I, Ih):
    pred = pdn4.predict_triggers(Ih, d['th_p'], d['lead']) if d['pred_on'] else None
    if d['arch'] == 'A' and not d.get('inj_pred', False) and d['boost'] == 0 and pred is not None and not d.get('inj'): pred = None
    ro = np.ones(len(I), np.uint8) if role_off == 'helper' else None
    return pdn4.simulate(I, p, spec, pred=pred, role_off=ro)

def evaluate(arch, g, full=False):
    d = decode(arch, g); p, spec, ro = build_design(d); v = d['var']
    cv = 0.0; ov = []; per = {}
    benches = WL if full else FIT
    for w in benches:
        I, Ih = trace(w, v); s = STARTS[w]
        seg = slice(s, None) if full else slice(s, s + WIN)
        m = sim(d, p, spec, ro, I[seg], Ih[seg])
        o = max(0.0, SUMM[w]['slow'][v] / 100) + SUMM[w]['dummy'][v] / 100 + m['loss_total'] + max(0.0, 2 * (m['mean_v'] - pdn4.VNOM) / pdn4.VNOM)
        ov.append(o); cv += max(0, m['p999'] - P999) / P999
        per[w] = dict(p999=m['p999'], max=m['max'], loss=m['loss'], o=o)
    wmax = 0.0
    if not full:
        for w in WL:
            I, Ih = trace(w, v)
            for a, b in EV[w][v]['win']:
                m = sim(d, p, spec, ro, I[a:b], Ih[a:b]); wmax = max(wmax, m['max'])
        cv += max(0, wmax - PMAX) / PMAX
    else:
        wmax = max(x['max'] for x in per.values()); cv += max(0, wmax - PMAX) / PMAX
    return area_of(d, p, spec), float(np.mean(ov)), cv, dict(per=per, wmax=wmax)

def _ev(args):
    arch, g = args
    try:
        a, o, cv, _ = evaluate(arch, g); return a, o, cv
    except Exception as e:
        return 1e9, 1e9, 1e9

def main(arch, pop=48, gens=60, seed=1, workers=3):
    import multiprocessing as mp
    NG = NC + NA[arch]; rng = np.random.RandomState(seed); X = rng.rand(pop, NG)
    for i, a in enumerate((0.3, 0.5, 0.7)):            # seeds: passive-ish starting points
        X[i, :NC] = 0; X[i, 0:4] = a * pdn4.SHARE / 0.5 * 4 / 4; X[i, 11] = 0.0
    pool = mp.get_context('spawn').Pool(workers)
    F = pool.map(_ev, [(arch, x) for x in X]); hist = []; out = 'ga4_%s.json' % arch
    for gen in range(gens):
        t0 = time.time(); rank, crowd = ga3.nsga_sort(F)
        def tour():
            i, j = rng.randint(pop, size=2)
            return i if (rank[i], -crowd[i]) < (rank[j], -crowd[j]) else j
        kids = []
        while len(kids) < pop:
            a, b = X[tour()], X[tour()]; u = rng.rand(NG); eta = 15.0
            beta = np.where(u < 0.5, (2 * u) ** (1 / (eta + 1)), (1 / (2 * (1 - u))) ** (1 / (eta + 1)))
            for c in (0.5 * ((1 + beta) * a + (1 - beta) * b), 0.5 * ((1 - beta) * a + (1 + beta) * b)):
                mut = rng.rand(NG) < 2.0 / NG; c[mut] += rng.normal(0, 0.15, mut.sum()); kids.append(np.clip(c, 0, 1))
        kids = np.array(kids[:pop]); Fk = pool.map(_ev, [(arch, x) for x in kids])
        XA = np.vstack([X, kids]); FA = F + Fk; rank, crowd = ga3.nsga_sort(FA)
        order = sorted(range(len(FA)), key=lambda i: (rank[i], -crowd[i]))[:pop]; X = XA[order]; F = [FA[i] for i in order]
        feas = sorted(f[:2] for f in F if f[2] == 0)
        hist.append(dict(gen=gen, feasible=len(feas), min_area=feas[0][0] if feas else None, sec=time.time() - t0))
        print('[%s] gen %2d feasible %2d front %s (%.0fs)' % (arch, gen, len(feas), ' '.join('(%.3f,%.3f)' % tuple(f) for f in feas[:6]), time.time() - t0), flush=True)
        json.dump(dict(arch=arch, X=X.tolist(), F=F, hist=hist), open(HERE + out, 'w'))
    pool.close()

if __name__ == '__main__':
    main(sys.argv[1], pop=int(sys.argv[2]) if len(sys.argv) > 2 else 48, gens=int(sys.argv[3]) if len(sys.argv) > 3 else 60)
