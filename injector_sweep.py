import json, numpy as np
from multiprocessing import Pool
import active_fix as A, active_fix3 as A3
from injector_sw import track
TN = '+5 nF damped (0.42 mm2)'
def job(a):
    wl, Ron, band = a
    import replay2
    from sweep import setup
    T = A.TIERS[TN]; inst, prods, cons, age, pin = setup(wl)
    RA = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0); m = RA['makespan']; I = RA['I'][:m]
    jr = A3.jref_of(A.quant_hat(RA, inst)[:m], 0.30)
    j, ntr, Ec, Esw = track(jr, R=Ron + 0.02, band=band, E_tr=32.5e-12 * 0.1 / Ron, Ron=Ron)
    Ts = m * 1e-9; d = A.droop(I - j * 1e3, T)
    return dict(wl=wl, Ron=Ron, band=band, base=A.droop(I, T)[0], p999=d[0], mx=d[1], P_mW=(Ec + Esw) / Ts * 1e3, sw_MHz=ntr / Ts / 1e6)
if __name__ == '__main__':
    jobs = [(wl, Ron, b) for wl in A.WL for Ron in (0.1, 0.3, 1.0) for b in (0.02, 0.05, 0.1)]
    with Pool(24) as p: R = p.map(job, jobs)
    json.dump(R, open('injector_sweep.json', 'w'), indent=1)
    for r in sorted(R, key=lambda r: (r['wl'], r['Ron'], r['band'])):
        print(f"{r['wl']:9s} Ron {r['Ron']:.1f} ohm  band +-{r['band']*1e3:3.0f} mA | p99.9 {r['base']:5.1f} -> {r['p999']:5.1f} (max {r['mx']:5.1f}) mV | loss {r['P_mW']:6.1f} mW = {100*r['P_mW']/180:5.1f}% of core | transitions {r['sw_MHz']:5.0f} MHz")
