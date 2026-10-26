"""Step 2c: injector with switching-level tracking (injector_sw.track): real transition counts and conduction loss."""
import json, numpy as np, time
from multiprocessing import Pool
import active_fix as A
from injector_sw import track
from scipy import signal
def jref_of(I_hat, Jmax=0.30, band=(50e6, 200e6)):
    b, a = signal.butter(2, [band[0] / 5e8, band[1] / 5e8], btype='bandpass')
    return np.clip(signal.filtfilt(b, a, (I_hat - I_hat.mean()) * 1e-3), -Jmax, Jmax)
def job(a):
    wl, tn, Jmax, bandw = a
    import replay2
    from sweep import setup
    T = A.TIERS[tn]; inst, prods, cons, age, pin = setup(wl)
    RA = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0); m = RA['makespan']; I = RA['I'][:m]
    jr = jref_of(A.quant_hat(RA, inst)[:m], Jmax, bandw)
    t0 = time.time(); j, ntr, Ec, Esw = track(jr); T_s = m * 1e-9
    P = (Ec + Esw) / T_s
    d = A.droop(I - j * 1e3, T); base = A.droop(I, T)
    track_err = float(np.sqrt(((j - jr) ** 2).mean()) / max(np.sqrt((jr ** 2).mean()), 1e-9))
    return dict(wl=wl, tier=tn, Jmax=Jmax, band=bandw, base=base, inj=d, P_mW=P * 1e3, Pcond_mW=Ec / T_s * 1e3, Psw_mW=Esw / T_s * 1e3,
                sw_rate_MHz=ntr / T_s / 1e6, track_err=track_err, t=time.time() - t0)
if __name__ == '__main__':
    jobs = [(wl, tn, Jm, bd) for wl in A.WL for tn in A.TIERS for Jm, bd in ((0.30, (50e6, 200e6)), (0.15, (50e6, 200e6)), (0.30, (60e6, 140e6)))]
    with Pool(24) as p: R = p.map(job, jobs)
    json.dump(R, open('active_fix3.json', 'w'), indent=1)
    for r in R:
        print(f"{r['wl']:9s} {r['tier']:38s} Jmax {r['Jmax']:.2f} band {r['band'][0]/1e6:.0f}-{r['band'][1]/1e6:.0f} | passive {r['base'][0]:5.1f} -> {r['inj'][0]:5.1f} / max {r['inj'][1]:5.1f} mV | "
              f"loss {r['P_mW']:5.1f} mW ({100*r['P_mW']/180:4.1f}% of core: cond {r['Pcond_mW']:.1f} + sw {r['Psw_mW']:.1f}), switching {r['sw_rate_MHz']:.0f} MHz, tracking err {100*r['track_err']:.0f}%")
