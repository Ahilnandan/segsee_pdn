"""B1 per-band droop readout and B3 jitter / timing model, applied to a recorded worst-tile voltage trace (one sample per
clock cycle, 333 MHz). Assumptions are stated here, not hidden:
  bands : droop x = (VNOM - v) split EXACTLY into band1 > 20 MHz (first droop), band2 0.5-20 MHz (second droop) and
          band3 < 0.5 MHz (buck loop + DC) with zero-phase low-pass filters; shares are reported at the p99.9 events and
          at the max.
  jitter: clock-tree insertion delay TAU_INS = 400 ps scaled by the alpha-power delay law (timing6.d); edge t lands
          e_t = TAU_INS (d(v_t)/d(VNOM) - 1) late; period jitter J_t = e_t - e_{t-1}; cycle-to-cycle jitter = J_t - J_{t-1}.
  timing: a chip-critical path (one full period at nominal) launched at edge t-1 and captured at edge t is late by
          T_CYC (d(v_t)/d(VNOM) - 1) - J_t (a late capture edge gives the data more time). Net timing loss = that, ps."""
import numpy as np
from scipy import signal
import timing6 as TM
VNOM = 1.1; T_CYC = 3e-9; FS = 1 / T_CYC; TAU_INS = 400e-12
_lo20 = signal.butter(4, 20e6 / (FS / 2), 'low', output='sos'); _lo05 = signal.butter(2, 0.5e6 / (FS / 2), 'low', output='sos')
def bands(v, warm=600):
    x = (VNOM - np.asarray(v, float)[warm:]) * 1e3
    a20 = signal.sosfiltfilt(_lo20, x); a05 = signal.sosfiltfilt(_lo05, x)
    b1, b2, b3 = x - a20, a20 - a05, a05
    thr = np.percentile(x, 99.9); ev = np.where(x >= thr)[0]; k = int(np.argmax(x))
    sh = lambda idx: [float(np.mean(bb[idx]) / max(np.mean(x[idx]), 1e-9)) for bb in (b1, b2, b3)]
    return dict(p999_shares=sh(ev), max_shares=sh([k]), p999=float(thr), max=float(x[k]),
                band_p999=[float(np.percentile(bb, 99.9)) for bb in (b1, b2, b3)])
def jitter(v, warm=600):
    v = np.asarray(v, float)[warm:]; rel = TM.d(v) / TM.d(VNOM) - 1
    e = TAU_INS * rel; J = np.diff(e); c2c = np.diff(J); late = T_CYC * rel[1:] - J
    q = lambda a: (float(np.percentile(np.abs(a), 99.9) * 1e12), float(np.abs(a).max() * 1e12))
    return dict(period_jitter_ps=q(J), c2c_jitter_ps=q(c2c), net_late_ps=(float(np.percentile(late, 99.9) * 1e12), float(late.max() * 1e12)),
                droop_only_late_ps=(float(np.percentile(T_CYC * rel, 99.9) * 1e12), float((T_CYC * rel).max() * 1e12)))
