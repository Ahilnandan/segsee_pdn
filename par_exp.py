"""Direct path + FIVR path in parallel (case D) vs direct only (A) vs FIVR only (C), full traces, closed loop."""
import json, numpy as np
from multiprocessing import Pool
import t3design
CF = dict(nph_fivr=16, L_fivr=112e-9, ron_f=0.2, dcr_f=0.2, pwm='tri', ctl='t3')
def dpar(fsw, fast, L_ph=112e-9):
    d = t3design.best_par(2 * fsw if fast else fsw / 1.8, L_ph=L_ph)
    return dict(wI=d['wI'], fz=(d['fz'], d['fz']), fp=(600e6, 1.2e9), Td=0.6e-9)
CFG = []
for fsw in (140, 70, 35):
    fc = 78e6 if fsw == 140 else fsw * 1e6 / 1.8
    CFG.append((f'C FIVR-only fc{fc/1e6:.0f}', dict(case='B', ideal_in=True, fsw_fivr=fsw * 1e6, t3=dict(fc=fc), **CF)))
    CFG.append((f'D par conservative', dict(case='D', fsw_fivr=fsw * 1e6, t3=dpar(fsw * 1e6, False), **CF)))
    if fsw >= 70: CFG.append((f'D par fast-loop', dict(case='D', fsw_fivr=fsw * 1e6, t3=dpar(fsw * 1e6, True), **CF)))
CFG.insert(0, ('A direct-only', dict(case='A')))
def job(a):
    wl, (name, kw), pol = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    g = run_gen2(inst, prods, cons, age, pin, pol, 100, 0, 25.0)
    o = simulate_chain(None, sched=g, nmax=n, cdie=5e-9, warmup=3000, **kw)
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']
    Iload = R['I'][:R['makespan']]
    out = dict(wl=wl, name=name, fsw=kw.get('fsw_fivr', 0) / 1e6, pol=pol, p999=float(np.percentile(v, 99.9)), worst=float(v.max()),
               top20=float(np.sort(v)[-20:].mean()), over=float(1e3 * (o['vmax'].max() - 0.9)), dc=float(1e3 * (o['vavg'].mean() - 0.9)),
               slow=100 * (R['makespan'] / R['treal'].max() - 1), Iload_mean=float(Iload.mean()))
    if 'ifivr' in o: out.update(ifivr_mean=float(o['ifivr'].mean() * 1e3), ifivr_rms=o['ifivr_rms'] * 1e3)
    return out
if __name__ == '__main__':
    jobs = [(wl, c, pol) for wl in ('qsort', 'dhrystone', 'mm') for c in CFG for pol in ('A', 'REGH')]
    print(len(jobs), 'jobs', [c[0] for c in CFG], flush=True)
    with Pool(23) as p, open('par_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
