"""Supply-induced clock jitter: baseline (FIVR-only 140 MHz, normal issue) vs our design (parallel 35 MHz + bank
+ predictive holds), single-thread BOOM closed loop in Pulsim. Jitter = clock-tree delay sensitivity x supply change
seen between consecutive edges; we report the supply statistics and convert with S = 0.35 ps/mV (300 ps tree, ~0.12%/mV)."""
import json, numpy as np
from multiprocessing import Pool
import par_exp, st_exp
KER = st_exp.KER
CFG = {'base: FIVR-only 140': (dict(case='B', ideal_in=True, fsw_fivr=140e6, t3=dict(fc=78e6), cdie=5e-9, **par_exp.CF), 'A', 0),
       'ours: par35+bank+MPCH': (dict(st_exp.KW, **st_exp.SUP['D35L4+bank']), 'MPCH', 0.68),
       'par35+bank, no holds': (dict(st_exp.KW, **st_exp.SUP['D35L4+bank']), 'A', 0)}
def job(a):
    wl, name = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    kw, pol, q = CFG[name]
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    mk = dict(mp_h=KER['D35L4'], mp_vth=q * st_exp.a999(wl), mp_maxd=19) if pol == 'MPCH' else {}
    g = run_gen2(inst, prods, cons, age, pin, 'A', 100, 0, 25.0, **mk)
    o = simulate_chain(None, sched=g, nmax=n, warmup=3000, **kw)
    va, vn, vx = o['vavg'] * 1e3, o['vmin'] * 1e3, o['vmax'] * 1e3
    d1 = np.abs(np.diff(va)); d4 = np.abs(va[4:] - va[:-4]); rip = vx - vn
    P = lambda x, q: float(np.percentile(x, q))
    return dict(wl=wl, name=name, dv1_p99=P(d1, 99), dv1_p999=P(d1, 99.9), dv1_max=float(d1.max()), dv1_rms=float(np.sqrt((d1**2).mean())),
                dv4_p999=P(d4, 99.9), rip_p50=P(rip, 50), rip_p999=P(rip, 99.9), droop_p999=P(900 - vn, 99.9))
if __name__ == '__main__':
    with Pool(9) as p: R = p.map(job, [(wl, nm) for wl in ('qsort', 'dhrystone', 'mm') for nm in CFG])
    json.dump(R, open('jitter_exp.json', 'w'), indent=1)
    S = 0.35
    for r in R:
        print(f"{r['wl']:9s} {r['name']:24s} | cycle-to-cycle dV p99 {r['dv1_p99']:5.2f} p99.9 {r['dv1_p999']:5.2f} max {r['dv1_max']:5.2f} mV "
              f"(-> {S*r['dv1_p999']:.2f} ps p99.9) | 4-cycle dV p99.9 {r['dv4_p999']:5.1f} | in-cycle ripple p50 {r['rip_p50']:4.2f} p99.9 {r['rip_p999']:5.2f} mV | droop p99.9 {r['droop_p999']:5.1f}")
