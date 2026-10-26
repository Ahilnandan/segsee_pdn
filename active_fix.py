"""Step 2: active fixes on top of each passive tier, real units, linear PDN (exact for the passive network).
 - holds: predictive holds (replay2 mp_h) with the kernel computed from THIS network's impedance (no identification needed)
 - injector: feed-forward. The issue stage knows each op's token ~13 cycles ahead, so the predicted current I_hat is known
   ahead of time. The injector sources/sinks j = -(50-200 MHz part of I_hat), zero-phase filtered within the lead window,
   limited by slew (V_in2 - V_DD)/L = 0.45 A/ns for 2 nH and by |j| <= J_max. Evaluated on the TRUE current (token error included:
   I_hat uses 4-bit-quantised token rectangles). Loss: conduction j_rms^2 R + switching 65 pJ per on/off transition pair.
Droop = -IFFT(Z(f) * FFT(I - j))."""
import json, numpy as np
from scipy import signal
import pbuck_freq as M, replay2
from sweep import setup
WL = ('qsort', 'dhrystone', 'mm')
TIERS = {'free (50% of 5 nF damped, 450 mOhm)': dict(p=dict(M.P, C_die=2.5e-9), Cd=2.5e-9, Rd=0.45),
         '+5 nF damped (0.42 mm2)': dict(p=dict(M.P), Cd=5e-9, Rd=0.2),
         '+10 nF damped, 2x bumps (0.83 mm2)': dict(p=dict(M.P, L_p=250e-12, R_p=M.P['R_p'] * 250 / 507), Cd=10e-9, Rd=0.2)}
def Zf(ff, T):
    s = 2j * np.pi * np.maximum(ff, 1.0); s0 = M.s; M.s = s
    Y = 1 / M.Zdie(T['p']) + 1 / M.Z_direct(T['p']) + 1 / (1 / (s * T['Cd']) + T['Rd'])
    M.s = s0; return 1 / Y
def kernel(T, K=256, n=8192, fhp=10e6):
    x = np.zeros(n); x[0] = 1e-3                          # 1 mA for one cycle
    ff = np.fft.rfftfreq(n, 1e-9)
    hp = (1j * ff / fhp) ** 2 / (1 + 1j * ff / fhp) ** 2  # 2nd-order high-pass: the issue stage only targets fast (>10 MHz) droop;
                                                          # slow IR droop through the 80 mOhm package is the board VR / guardband's job
    v = -np.fft.irfft(np.fft.rfft(x) * Zf(ff, T) * hp, n) * 1e3   # mV droop (positive = down)
    return -v[:K] / 1.0                                   # replay2 convention: kernel in mV per mA, voltage deviation (negative = droop)
def droop(I, T):
    x = (I - I.mean()) * 1e-3; n = len(x); ff = np.fft.rfftfreq(n, 1e-9)
    v = -np.fft.irfft(np.fft.rfft(x) * Zf(ff, T), n) * 1e3
    return float(np.percentile(v, 99.9)), float(v.max())
def quant_hat(R, inst):
    """token-estimated current: each issued op as a 4-bit-amplitude rectangle over its duration (what the issue stage knows)"""
    from curmodel import SHAPE, FLOOR
    Ih = np.full(len(R['I']), FLOOR); step = 80.0 / 15
    for i, t in enumerate(R['tnew']):
        sh = SHAPE.get(inst[i][1], SHAPE[0]); d = len(sh); a = min(15, round(sh.sum() / d / step)) * step
        Ih[t:t + d] += a
    return Ih
def injector(I_true, I_hat, L=2e-9, Jmax=0.30, band=(50e6, 200e6), lead=13):
    b, a = signal.butter(2, [band[0] / 5e8, band[1] / 5e8], btype='bandpass')
    hp = signal.filtfilt(b, a, (I_hat - I_hat.mean()) * 1e-3)          # zero-phase: needs future samples -> within the lead window
    j = np.clip(hp, -Jmax, Jmax)
    slew = (1.8 - 0.9) / L * 1e-9                                       # A per cycle
    jj = np.zeros_like(j)
    for k in range(1, len(j)): jj[k] = jj[k - 1] + np.clip(j[k] - jj[k - 1], -slew, slew)
    on = np.abs(jj) > 0.01; trans = int(np.sum(on[1:] != on[:-1]))
    Ploss = (jj ** 2).mean() * (0.1 + 0.02) + trans * 65e-12 / (len(jj) * 1e-9)
    return I_true - jj * 1e3, Ploss
if __name__ == '__main__':
    R = {}
    for wl in WL:
        inst, prods, cons, age, pin = setup(wl)
        RA = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0); mA = RA['makespan']; IA = RA['I'][:mA]
        IhA = quant_hat(RA, inst)[:mA]
        for tn, T in TIERS.items():
            h = kernel(T); base = droop(IA, T)
            a999 = float(np.percentile(-np.convolve(IA - IA.mean(), h)[:mA], 99.9))
            RH = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0, mp_h=h, mp_vth=0.68 * a999, mp_maxd=19)
            mH = RH['makespan']; IH = RH['I'][:mH]; hold = droop(IH, T); slowH = 100 * (mH / mA - 1)
            Ii, Pi = injector(IA, IhA); inj = droop(Ii, T)
            IHh = quant_hat(RH, inst)[:mH]; Ib, Pb = injector(IH, IHh); both = droop(Ib, T)
            R[f'{wl}|{tn}'] = dict(base=base, holds=hold, slow_holds=slowH, injector=inj, inj_loss_mW=Pi * 1e3, both=both, both_loss_mW=Pb * 1e3)
            print(f"{wl:9s} {tn:38s} | passive {base[0]:5.1f}/{base[1]:5.1f} | +holds {hold[0]:5.1f}/{hold[1]:5.1f} ({slowH:+.2f}%) | "
                  f"+injector {inj[0]:5.1f}/{inj[1]:5.1f} ({Pi*1e3:.1f} mW) | +both {both[0]:5.1f}/{both[1]:5.1f}", flush=True)
    json.dump(R, open('active_fix.json', 'w'), indent=1)
