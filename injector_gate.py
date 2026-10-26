"""Gated injector: track the band-pass reference only while the predictor (high-passed kernel of THIS network, on the
token-estimated current) sees a predicted dip/overshoot above a threshold within the next 13 cycles (window dilated +-8 cycles).
Switches 1 ohm, hysteresis +-100 mA, 2 nH. Tier: +5 nF damped decap."""
import json, numpy as np
from multiprocessing import Pool
from scipy import ndimage
import active_fix as A, active_fix3 as A3
from injector_sw import track
TN = '+5 nF damped (0.42 mm2)'
def job(a):
    wl, gate_q = a
    import replay2
    from sweep import setup
    T = A.TIERS[TN]; inst, prods, cons, age, pin = setup(wl)
    RA = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0); m = RA['makespan']; I = RA['I'][:m]
    Ih = A.quant_hat(RA, inst)[:m]; jr = A3.jref_of(Ih, 0.30)
    if gate_q is not None:
        h = A.kernel(T); vp = np.convolve(Ih - Ih.mean(), h)[:m]          # predicted voltage deviation (mV, negative = dip)
        look = ndimage.maximum_filter1d(np.abs(vp), 13, origin=-6)          # worst |deviation| over the next 13 cycles
        thr = np.percentile(np.abs(vp), gate_q); mask = ndimage.binary_dilation(look > thr, iterations=8)
        jr = jr * mask; duty = float(mask.mean())
    else: duty = 1.0
    j, ntr, Ec, Esw = track(jr, R=1.02, band=0.1, E_tr=3.25e-12, Ron=1.0)
    Ts = m * 1e-9; d = A.droop(I - j * 1e3, T)
    return dict(wl=wl, gate=gate_q, duty=duty, base=A.droop(I, T), inj=d, P_mW=(Ec + Esw) / Ts * 1e3)
if __name__ == '__main__':
    jobs = [(wl, g) for wl in A.WL for g in (None, 80, 90, 95, 98)]
    with Pool(15) as p: R = p.map(job, jobs)
    json.dump(R, open('injector_gate.json', 'w'), indent=1)
    for r in R:
        print(f"{r['wl']:9s} gate {'always' if r['gate'] is None else 'top ' + str(100 - r['gate']) + '%':8s} active {100*r['duty']:5.1f}% | p99.9 {r['base'][0]:5.1f} -> {r['inj'][0]:5.1f}, max {r['base'][1]:5.1f} -> {r['inj'][1]:5.1f} mV | loss {r['P_mW']:5.1f} mW ({100*r['P_mW']/180:.1f}%)")
