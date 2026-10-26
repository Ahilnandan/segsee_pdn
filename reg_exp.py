"""Can the CPU regulate its own supply? FIVR with its feedback loop REMOVED (fixed duty), the
scheduler's REG policy (hold on deficit, dummy-shunt on surplus) as the only actuator."""
import json, sys, numpy as np
from multiprocessing import Pool
def job(a):
    wl, fsw, ctl, pol, gain = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    g = run_gen2(inst, prods, cons, age, pin, pol, 100, 0, 0)
    o = simulate_chain(None, case='B', sched=g, nmax=n, fsw_fivr=fsw * 1e6, L_fivr=2e-9 * 140 / fsw, ctl=ctl, sig_gain=gain)
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']
    return dict(wl=wl, fsw=fsw, ctl=ctl, pol=pol, gain=gain, droop=float(v.max()), p999=float(np.percentile(v, 99.9)),
                over=float(1e3 * (o['vmax'].max() - 0.9)), vavg=float(o['vavg'].mean()), vstd=float(1e3 * o['vavg'].std()),
                slow=100 * (R['makespan'] / R['treal'].max() - 1), dummyE=100 * R['dummyq'] / (R['I'].sum() - R['dummyq']))
if __name__ == '__main__':
    jobs = [(wl, fsw, ctl, pol, 1.0) for wl in ('qsort', 'dhrystone') for fsw in (140, 70)
            for ctl, pol in (('vm', 'A'), ('open', 'A'), ('open', 'REG'), ('vm', 'REG'))]
    with Pool(6) as p, open('reg_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush(); print(r, flush=True)
