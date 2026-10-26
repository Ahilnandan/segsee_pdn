"""Attack the open gaps against their prior-art counterparts (see docs/prior-art-report.md)."""
import json, sys, numpy as np
from multiprocessing import Pool
B = dict(cb=0.5, ir=25.0)
def cfg(name):
    # L1: overshoot handling at equal area (4.5 nF passive + 0.5 nF bank), same trigger
    if name == 'L1-none':  return dict(bank=dict(B))
    if name == 'L1-dissip': return dict(bank=dict(B, sink_dissipate=True))     # prior art (P14-P16, A9)
    if name == 'L1-store':  return dict(bank=dict(B, harvest_v=True))           # ours
    if name == 'L1-hybrid': return dict(bank=dict(B, hybrid=True))              # store until full, then dissipate
    raise KeyError(name)
def job(a):
    wl, exp, var, pol, fsw = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    kw = dict(case='B', nmax=n, fsw_fivr=fsw * 1e6, L_fivr=2e-9 * 140 / fsw, cdie=4.5e-9)
    ir = 25.0
    if exp == 'L1':
        kw.update(cfg(var)); kw.update(ctl='vm')
    elif exp == 'L2':                       # same actuator (REGH holds) + storing bank, different sensing
        kw.update(bank=dict(B, harvest_v=True), ctl='vm', kp_vm=0.0, fint=0.5e6, sig_mode=var)
    elif exp == 'L3':                       # TD hold gated by modeled tank (P3-style) vs sensed reservoir
        b = dict(B)
        if var != 'sensed': b['tok_model'] = dict(Qb=float(var), Ir=ir)
        kw.update(bank=b, ctl='vm')
    g = run_gen2(inst, prods, cons, age, pin, pol, 100, 0, ir)
    o = simulate_chain(None, sched=g, **kw)
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']
    return dict(wl=wl, exp=exp, var=var, pol=pol, fsw=fsw, droop=float(v.max()), p999=float(np.percentile(v, 99.9)),
                over=float(1e3 * (o['vmax'].max() - 0.9)), slow=100 * (R['makespan'] / R['treal'].max() - 1),
                Ediss_nJ=o['E']['Ediss'] * 1e9, Ein_nJ=o['E']['Ein'] * 1e9, Estored_nJ=o['E']['Estored'] * 1e9,
                Eload_nJ=float((R['I'][:R['makespan']] * 1e-3 * 0.9 * 1e-9).sum() * 1e9))
if __name__ == '__main__':
    jobs = []
    if len(sys.argv) > 1 and sys.argv[1] == 'hybrid':
        jobs = [(wl, 'L1', 'L1-hybrid', pol, fsw) for wl in ('qsort', 'dhrystone', 'mm') for fsw in (140, 70) for pol in ('A', 'REGH')]
        with Pool(12) as p, open('gaps_hybrid.jsonl', 'w') as f:
            for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
        sys.exit(0)
    for wl in ('qsort', 'dhrystone', 'mm'):
        for fsw in (140, 70):
            for var in ('L1-none', 'L1-dissip', 'L1-store'):
                for pol in ('A', 'REGH'): jobs.append((wl, 'L1', var, pol, fsw))
            for var in ('v', 'i', 'both'): jobs.append((wl, 'L2', var, 'REGH', fsw))
        for var in ('sensed', '250', '125', '60', '30'): jobs.append((wl, 'L3', var, 'TD', 140))
    with Pool(12) as p, open('gaps.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
