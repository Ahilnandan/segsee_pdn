"""Equal inductor + switch budget: FIVR phases traded for inductance (parallel supply).
Total inductance fixed (16 x 112 nH = 1.79 uH); DCR per phase proportional to L per phase (same total copper);
switch Ron per phase = 0.2 ohm * 16/N (same total switch area) -> switching energy per phase x16/N, so FIVR
switching power = 16 * fsw * E_cyc(16-phase calibration) independent of N. Single-thread BOOM, closed loop, policy A."""
import json, numpy as np
from multiprocessing import Pool
import t3design, par_exp
ECYC16 = 2.72e-12          # J per phase-cycle at 16 phases (calibrated: FIVR-only@140 = 90% efficient)
def mk(N, Lph, fsw, case='D', dcr_scale=True):
    ron = 0.2 * 16 / N; dcr = 0.2 * (Lph / 112e-9 if dcr_scale else 1.0)
    base = dict(nph_fivr=N, L_fivr=Lph, ron_f=ron, dcr_f=dcr, pwm='tri', ctl='t3', fsw_fivr=fsw)
    if case == 'D':
        d = t3design.best_par(fsw / 1.8, L_ph=Lph, nph=N, ron=ron, dcr=dcr)
        base['t3'] = dict(wI=d['wI'], fz=(d['fz'], d['fz']), fp=(600e6, 1.2e9), Td=0.6e-9)
        return dict(case='D', **base), d
    return dict(case='B', ideal_in=True, t3=dict(fc=78e6), **base), None
CFG = [('FIVR-only 16x112nH @140', 16, 112e-9, 140e6, 'B', True),
       ('par 16x112nH @140', 16, 112e-9, 140e6, 'D', True),
       ('par 16x112nH @35 (same parts)', 16, 112e-9, 35e6, 'D', True),
       ('par 8x224nH @70 (iso-budget)', 8, 224e-9, 70e6, 'D', True),
       ('par 4x448nH @35 (iso-budget)', 4, 448e-9, 35e6, 'D', True),
       ('par 16x448nH @35 (4x budget, DCR x4)', 16, 448e-9, 35e6, 'D', True),
       ('par 16x448nH @35 (4x L, DCR fixed: old optimistic)', 16, 448e-9, 35e6, 'D', False)]
def job(a):
    wl, (name, N, Lph, fsw, case, dsc) = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    kw, d = mk(N, Lph, fsw, case, dsc)
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    g = run_gen2(inst, prods, cons, age, pin, 'A', 100, 0, 25.0)
    o = simulate_chain(None, sched=g, nmax=n, cdie=5e-9, warmup=3000, **kw)
    E = o['EN']; T = o['T_meas']; v = 1e3 * (0.9 - o['vmin']); rip = 1e3 * (o['vmax'] - o['vmin'])
    E12, E18, El = -E['E12'], -E['E18'], E['Eload']
    psw = 16 * fsw * ECYC16 * T
    Ein = (E12 / 0.95 if case == 'D' else 0.0) + (E18 + psw) / 0.90
    return dict(wl=wl, name=name, N=N, Lph=Lph, fsw=fsw, p999=float(np.percentile(v, 99.9)), worst=float(v.max()),
                rip50=float(np.median(rip)), eta=El / Ein, Psw_mW=psw / T * 1e3, ifr=o.get('ifivr_rms', 0) * 1e3,
                chatter=E['nsw_f'] / max(1.0, 2 * N * fsw * T), loop=None if d is None else {k: float(v) for k, v in d.items()})
if __name__ == '__main__':
    for c in CFG[1:]:
        _, d = mk(c[1], c[2], c[3], c[4], c[5]); print(c[0], 'loop:', None if d is None else f"Zpk {d['Zpk']*1e3:.0f} mOhm fcmax {d['fcmax']/1e6:.1f} MHz GM {d['gm']:.1f}", flush=True)
    jobs = [(wl, c) for wl in ('qsort', 'dhrystone', 'mm') for c in CFG]
    with Pool(21) as p, open('isoL_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
    R = [json.loads(l) for l in open('isoL_exp.jsonl')]
    for r in sorted(R, key=lambda r: (r['wl'], [c[0] for c in CFG].index(r['name']))):
        print(f"{r['wl']:9s} {r['name']:46s} | p99.9 {r['p999']:6.1f} worst {r['worst']:6.1f} ripple {r['rip50']:5.2f} mV | eta {100*r['eta']:5.1f}% | "
              f"FIVR sw {r['Psw_mW']:4.2f} mW, I_rms {r['ifr']:5.1f} mA, chatter x{r['chatter']:.2f}")
