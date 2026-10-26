"""Deferred + useful filler work, single-thread BOOM closed loop in Pulsim on the final supply (parallel, FIVR 35 MHz 4xL,
active bank). Elastic work bank = background work needed anyway (e.g. AES-CTR memory-encryption pads)."""
import json, numpy as np
from multiprocessing import Pool
import st_exp
KER = st_exp.KER; KW = st_exp.KW; SUPK = st_exp.SUP['D35L4+bank']
h = KER['D35L4']
CFG = {
 'A':                    dict(),
 'MPCH':                 dict(mp=0.68),
 'wb const r0.5':        dict(wb=dict(mode='const', rate=0.5, cap=64, low=8)),
 'wb mpc r0.5':          dict(wb=dict(mode='mpc', rate=0.5, cap=64, low=8)),
 'wb mpc r0.5 +MPCH':    dict(mp=0.68, wb=dict(mode='mpc', rate=0.5, cap=64, low=8)),
 'wb const r1.0 +MPCH':  dict(mp=0.68, wb=dict(mode='const', rate=1.0, cap=64, low=8)),
 'wb mpc r1.0 +MPCH':    dict(mp=0.68, wb=dict(mode='mpc', rate=1.0, cap=64, low=8)),
 'redund vfill5 +MPCH':  dict(mp=0.68, wb=dict(mode='redund', vfill=5.0)),
 'redund vfill15 +MPCH': dict(mp=0.68, wb=dict(mode='redund', vfill=15.0)),
 'dummy fill +MPCH':     dict(mp=0.68, dummy='fill'),
}
def job(a):
    wl, name = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    c = CFG[name]
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    kw = dict(mp_h=h)
    if 'mp' in c: kw.update(mp_vth=c['mp'] * st_exp.a999(wl), mp_maxd=19)
    if 'wb' in c: kw['wb'] = c['wb']
    if 'dummy' in c: kw.update(dummy='fill', THF=200.0)
    g = run_gen2(inst, prods, cons, age, pin, 'A', 100, 0, 25.0, **kw)
    o = simulate_chain(None, sched=g, nmax=n, warmup=3000, **KW, **SUPK)
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']; m = R['makespan']
    real = sum(float(np.sum(s)) for s in [])  # placeholder
    w = R.get('wb') or {}
    Itot = float(R['I'][:m].sum())
    extra = (w.get('q', 0.0) + R['dummyq']) / Itot
    return dict(wl=wl, name=name, slow=100 * (m / R['treal'].max() - 1), p999=float(np.percentile(v, 99.9)), worst=float(v.max()),
                over=float(1e3 * (o['vmax'].max() - 0.9)), extra_pct=100 * extra, starve=w.get('starve', 0.0), bank_ops=w.get('n', 0),
                cov=w.get('cov', 0), nint=sum(1 for x in inst if x[0] == 1), dummyq=R['dummyq'])
if __name__ == '__main__':
    jobs = [(wl, nm) for wl in ('qsort', 'dhrystone', 'mm') for nm in CFG]
    with Pool(23) as p, open('wb_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
    R = [json.loads(l) for l in open('wb_exp.jsonl')]; ref = {r['wl']: r for r in R if r['name'] == 'A'}
    for r in sorted(R, key=lambda r: (r['wl'], list(CFG).index(r['name']))):
        b = ref[r['wl']]
        print(f"{r['wl']:9s} {r['name']:22s} | slow {r['slow']:+5.2f}% | p99.9 {r['p999']:5.1f} ({100*(r['p999']/b['p999']-1):+5.1f}%) worst {r['worst']:5.1f} "
              f"({100*(r['worst']/b['worst']-1):+5.1f}%) over {r['over']:5.1f} ({100*(r['over']/b['over']-1):+5.1f}%) | filler charge {r['extra_pct']:4.1f}% "
              f"starve {r['starve']:.0f} cov {100*r['cov']/r['nint']:.0f}% of int ops")
