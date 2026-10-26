"""Identify a per-cycle impulse-response kernel of a supply chain in Pulsim (for the predictive scheduler).
Differential method: (DC + sparse 1-cycle pulses) minus (DC only), averaged over pulses at varied PWM phase."""
import numpy as np, sys, json
from pdnchain import simulate_chain
import par_exp
def kernel(kw, K=256, base=150.0, amp=100.0, npulse=40, gap=1500, seed=0):
    rng = np.random.default_rng(seed)
    pos = np.cumsum(gap + rng.integers(0, 37, npulse)) + 2000
    n = int(pos[-1] + K + 200)
    I0 = np.full(n, base); I1 = I0.copy(); I1[pos] += amp
    o0 = simulate_chain(I0, warmup=3000, **kw); o1 = simulate_chain(I1, warmup=3000, **kw)
    d = (o1['vavg'] - o0['vavg']) * 1e3          # mV
    R = np.array([d[p:p + K] for p in pos])
    return R.mean(0) / amp, R.std(0) / amp       # mV per mA per 1-cycle pulse
if __name__ == '__main__':
    CF = par_exp.CF; out = {}
    for name, kw in [('D35', dict(case='D', fsw_fivr=35e6, t3=par_exp.dpar(35e6, False), **CF)),
                     ('D140', dict(case='D', fsw_fivr=140e6, t3=par_exp.dpar(140e6, False), **CF)),
                     ('C140', dict(case='B', ideal_in=True, fsw_fivr=140e6, t3=dict(fc=78e6), **CF))]:
        h, sd = kernel(kw); out[name] = h.tolist()
        step = np.cumsum(h)
        print(f"{name}: peak {h.min():.4f} mV/mA @ cycle {h.argmin()}, step-response min {step.min():.3f} mV/mA @ {step.argmin()}, "
              f"final {step[-1]:.4f}, sum|h| {np.abs(h).sum():.3f}, noise (sd/|peak|) {sd[:20].mean()/abs(h.min()):.2f}", flush=True)
    json.dump(out, open('kernels.json', 'w'))
