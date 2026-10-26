"""Re-test the CPU-as-regulator / bank claims against a FAITHFUL Intel-class FIVR (case C preset):
16 phases, natural-sampled double-edge PWM, type-3 compensator UGF 78 MHz / PM 43 deg incl. 0.6 ns loop delay.
Full traces, closed loop (scheduler generator inside the switching sim), off-chip 12->1.8 V VR included."""
import json, sys, numpy as np
from multiprocessing import Pool
CASEC = dict(case='B', nph_fivr=16, L_fivr=112e-9, ron_f=0.2, dcr_f=0.2, pwm='tri')
BANK = dict(cb=0.5, ir=25.0)
def cfgs():
    out = []
    for fsw in (140, 70, 35):
        bw = (78e6,) if fsw == 140 else (78e6, fsw * 1e6 / 1.8)       # fsw/1.8 = Intel's 78/140 ratio
        for fc in bw:
            t3 = dict(t3=dict(fc=fc), ctl='t3')
            out += [(f't3 fc{fc/1e6:.0f}', fsw, 'A', t3, None),
                    (f't3 fc{fc/1e6:.0f}', fsw, 'A', t3, 'bidir'),
                    (f't3 fc{fc/1e6:.0f}', fsw, 'A', t3, 'hybrid'),
                    (f't3 fc{fc/1e6:.0f}', fsw, 'REGH', t3, 'bidir')]
        slow = dict(ctl='vm', kp_vm=0.0, fint=0.5e6)
        out += [('slowint', fsw, 'REGH', slow, 'bidir'), ('slowint', fsw, 'REG', slow, 'bidir')]
    return out
def job(a):
    wl, (name, fsw, pol, ctl, bk) = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    kw = dict(CASEC, nmax=n, fsw_fivr=fsw * 1e6, **ctl)
    if bk == 'bidir': kw.update(bank=dict(BANK, harvest_v=True), cdie=4.5e-9)
    elif bk == 'hybrid': kw.update(bank=dict(BANK, hybrid=True), cdie=4.5e-9)
    else: kw.update(cdie=5e-9)
    g = run_gen2(inst, prods, cons, age, pin, pol, 100, 0, 25.0)
    o = simulate_chain(None, sched=g, **kw)
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']
    return dict(wl=wl, name=name, fsw=fsw, pol=pol, bank=bk, droop=float(v.max()), p999=float(np.percentile(v, 99.9)),
                top20=float(np.sort(v)[-20:].mean()), over=float(1e3 * (o['vmax'].max() - 0.9)),
                dc=float(1e3 * (o['vavg'].mean() - 0.9)), slow=100 * (R['makespan'] / R['treal'].max() - 1),
                ndum=int(R.get('ndum', 0) or 0), dummyq=float(R.get('dummyq', 0.0)), Iq=float(R['I'][:R['makespan']].sum()), Ediss_nJ=o['E']['Ediss'] * 1e9, Ein_nJ=o['E']['Ein'] * 1e9,
                Eload_nJ=float((R['I'][:R['makespan']] * 1e-3 * 0.9 * 1e-9).sum() * 1e9))
if __name__ == '__main__':
    jobs = [(wl, c) for c in cfgs() for wl in ('qsort', 'dhrystone', 'mm')]
    print(len(jobs), 'jobs', flush=True)
    with Pool(20) as p, open('fivr_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush(); print(r, flush=True)
