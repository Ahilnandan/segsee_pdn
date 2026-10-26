"""Single-thread BOOM (real ports/deps, PaperMatchBoomV3Config traces), CLOSED LOOP inside the Pulsim circuit:
parallel supply (FIVR 35 MHz, 4x L) with/without active bank; oldest-first vs reactive holds (REGH) vs
predictive holds (MPCH, identified kernel)."""
import json, numpy as np
from multiprocessing import Pool
import par_exp
KER = json.load(open('kernels.json'))
KW = dict(case='D', fsw_fivr=35e6, t3=par_exp.dpar(35e6, False, 448e-9), **dict(par_exp.CF, L_fivr=448e-9))
SUP = {'D35L4': dict(cdie=5e-9), 'D35L4+bank': dict(cdie=4.5e-9, bank=dict(cb=0.5, ir=25.0))}
_A = {}
def a999(wl):
    if wl in _A: return _A[wl]
    import replay2
    from sweep import setup
    inst, prods, cons, age, pin = setup(wl)
    R = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0); m = R['makespan']; I = R['I'][:m]
    _A[wl] = float(np.percentile(-np.convolve(I - I.mean(), KER['D35L4'])[:m], 99.9)); return _A[wl]
def job(a):
    wl, sup, pol, q = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    kw = dict(mp_h=KER['D35L4'], mp_vth=q * a999(wl), mp_maxd=19) if pol == 'MPCH' else {}
    g = run_gen2(inst, prods, cons, age, pin, 'A' if pol == 'MPCH' else pol, 100, 0, 25.0, **kw)
    o = simulate_chain(None, sched=g, nmax=n, warmup=3000, **KW, **SUP[sup])
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']
    return dict(wl=wl, sup=sup, pol=pol, q=q, slow=100 * (R['makespan'] / R['treal'].max() - 1), nhold=R['nhold'],
                p999=float(np.percentile(v, 99.9)), worst=float(v.max()), over=float(1e3 * (o['vmax'].max() - 0.9)))
if __name__ == '__main__':
    jobs = [(wl, s, p, q) for wl in ('qsort', 'dhrystone', 'mm') for s in SUP
            for p, q in [('A', 0), ('REGH', 0), ('MPCH', 0.68), ('MPCH', 0.8), ('MPCH', 0.9)]]
    with Pool(23) as p, open('st_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
    R = [json.loads(l) for l in open('st_exp.jsonl')]
    ref = {r['wl']: r for r in R if r['pol'] == 'A' and r['sup'] == 'D35L4'}
    for r in sorted(R, key=lambda r: (r['wl'], r['sup'], r['pol'], r['q'])):
        b = ref[r['wl']]
        print(f"{r['wl']:9s} {r['sup']:11s} {r['pol']:4s} {r['q']:4.2f} | slow {r['slow']:+5.2f}% holds {r['nhold']:6d} | p99.9 {r['p999']:6.1f} ({100*(r['p999']/b['p999']-1):+5.1f}%) "
              f"worst {r['worst']:6.1f} ({100*(r['worst']/b['worst']-1):+5.1f}%) over {r['over']:6.1f}")
