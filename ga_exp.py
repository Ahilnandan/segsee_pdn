"""Genetic-algorithm tuning of droop controllers, for a FAIR head-to-head:
  OURS : predictive (MPC) issue holds + optional reorder underneath, bank-SOC-aware
  GRO  : Grochowski/Ayers/Tiwari HPCA'02 reactive convolution controller (lower clamp: stop issue;
         optional upper clamp: force units on = dummy ops, energy penalised)
Each family is tuned by the GA to its best droop at <=1% slowdown on TRAINING slices, then the champions
are validated on HELD-OUT slices/mixes. Both use the same realistic hardware view: 4-bit-amp x 4-bit-dur
power tokens and a PDN state sensed LAT cycles late.

fitness (lower = better) = p99.9/A + 0.5*worst/A + 10*max(0, slow% - 1) + 0.05*dummy-energy%  (per mix, averaged)
"""
import json, sys, random, numpy as np
from multiprocessing import Pool
import smt_exp, smt_exp2
from smt_exp import MIXES
from soc_exp import MINB

LAT = 2; FRAC = 0.5
TRAIN = [('smt8-mix', 0, 20000), ('smt4-mix', 0, 20000), ('smt8-mix', 30000, 20000)]
SELS = ['A', 'T', 'SOC', 'MPC']
_cache = {}

def setup(mix, shift, n):
    key = (mix, shift, n)
    if key in _cache: return _cache[key]
    from smt_replay import make_thread, run, pdn_predictor
    from replay2 import droop
    thr = [make_thread(w, s + shift, n) for w, s in MIXES[mix]]
    mb, k = MINB[mix]; Qb = FRAC * mb
    Lp, Cp, Rp = 507e-12 / k, 5e-9 * k, 79.6e-3 / k
    base = dict(TH=100.0 * k, Qb=Qb, Ir=25.0 * k, record_ready=False, est=dict(quant=True, lat=LAT))
    RA = run(thr, 'A', pdn=pdn_predictor(L=Lp, C=Cp, R=Rp), **base)
    vA = droop(RA['In'], L=Lp, C=Cp, R=Rp)
    d = dict(thr=thr, base=base, LCR=(Lp, Cp, Rp), mkA=RA['makespan'], a999=float(np.percentile(vA, 99.9)),
             aw=float(vA.max()), a50=float(np.percentile(vA, 50)))
    _cache[key] = d; return d

def evaluate(fam, g, mix, shift, n):
    from smt_replay import run, pdn_predictor
    from replay2 import droop
    S = setup(mix, shift, n); Lp, Cp, Rp = S['LCR']
    if fam == 'OURS':
        P = pdn_predictor(L=Lp, C=Cp, R=Rp, Hp=int(g['Hp']))
        R = run(S['thr'], g['sel'], pdn=P, res=g['res'], vth=g['q'] * S['a999'], hold='MPCH', MAXD=int(g['maxd']), **S['base'])
    else:
        P = pdn_predictor(L=Lp, C=Cp, R=Rp)
        gro = dict(lo=g['lo'] * S['a999'], hi=(g['hi'] * S['a999'] if g['use_hi'] else None))
        R = run(S['thr'], 'A', pdn=P, gro=gro, MAXD=int(g['maxd']), **S['base'])
    v = droop(R['In'], L=Lp, C=Cp, R=Rp)
    slow = 100 * (R['makespan'] / S['mkA'] - 1); dE = 100 * R['dumq'] / R['Itot']
    p999 = float(np.percentile(v, 99.9)); w = float(v.max())
    fit = p999 / S['a999'] + 0.5 * w / S['aw'] + 10 * max(0.0, slow - 1.0) + 0.05 * dE
    return dict(fit=fit, p999r=p999 / S['a999'] - 1, worstr=w / S['aw'] - 1, slow=slow, dE=dE, nheld=R['nheld'])

def fitness(a):
    fam, g = a
    rs = [evaluate(fam, g, *tr) for tr in TRAIN]
    return dict(g=g, fit=float(np.mean([r['fit'] for r in rs])), per=rs)

SPACE = {'OURS': dict(q=(0.55, 1.05), res=(0.0, 0.8), Hp=(10, 20), maxd=(2, 32), sel=SELS),
         'GRO': dict(lo=(0.4, 1.2), hi=(0.0, 0.6), use_hi=[0, 1], maxd=(2, 32))}

def rand(fam):
    g = {}
    for k, v in SPACE[fam].items():
        g[k] = random.choice(v) if isinstance(v, list) else random.uniform(*v)
    return g

def mutate(fam, g, rate=0.3):
    g = dict(g)
    for k, v in SPACE[fam].items():
        if random.random() < rate:
            if isinstance(v, list): g[k] = random.choice(v)
            else:
                lo, hi = v; g[k] = min(hi, max(lo, g[k] + random.gauss(0, 0.15 * (hi - lo))))
    return g

def cross(fam, a, b):
    c = {}
    for k, v in SPACE[fam].items():
        if isinstance(v, list): c[k] = random.choice([a[k], b[k]])
        else: u = random.random(); c[k] = u * a[k] + (1 - u) * b[k]      # blend
    return c

def ga(fam, pool, pop=48, gens=20, seed=1, log=None):
    random.seed(seed)
    P = [rand(fam) for _ in range(pop)]
    if fam == 'OURS': P[0] = dict(q=0.8, res=0.3, Hp=12, maxd=8, sel='A')        # hand-tuned seed
    else: P[0] = dict(lo=0.8, hi=0.2, use_hi=0, maxd=8)
    best = None
    for gen in range(gens):
        R = pool.map(fitness, [(fam, g) for g in P])
        R.sort(key=lambda r: r['fit'])
        if best is None or R[0]['fit'] < best['fit']: best = R[0]
        if log:
            log.write(json.dumps(dict(fam=fam, gen=gen, best=R[0], median=float(np.median([r['fit'] for r in R])))) + '\n'); log.flush()
        print(fam, 'gen', gen, 'best', round(R[0]['fit'], 4), R[0]['g'], [(round(p['p999r'] * 100, 1), round(p['slow'], 2)) for p in R[0]['per']], flush=True)
        elite = [r['g'] for r in R[:2]]
        def tour():
            a, b = random.sample(R, 2); return (a if a['fit'] < b['fit'] else b)['g']
        P = elite + [mutate(fam, cross(fam, tour(), tour())) for _ in range(pop - 2)]
    return best

if __name__ == '__main__':
    with Pool(24) as pool, open('ga_log.jsonl', 'w') as log:
        champs = {fam: ga(fam, pool, log=log) for fam in ('OURS', 'GRO')}
        json.dump({k: v['g'] for k, v in champs.items()}, open('ga_champions.json', 'w'), indent=1)
        # held-out validation: different slices AND unseen mixes, full length, both bank sizes
        VAL = [('smt8-b', 0, 60000), ('smt8-c', 0, 60000), ('smt8-mix', 0, 60000), ('smt4-mix', 0, 60000)]
        jobs = [(fam, champs[fam]['g'], v) for fam in champs for v in VAL]
        res = pool.starmap(evaluate, [(f, g, *v) for f, g, v in jobs])
        with open('ga_validate.jsonl', 'w') as f:
            for (fam, g, v), r in zip(jobs, res):
                f.write(json.dumps(dict(fam=fam, g=g, mix=v[0], **r)) + '\n')
                print('VAL', fam, v[0], {k: round(x, 3) if isinstance(x, float) else x for k, x in r.items()}, flush=True)
