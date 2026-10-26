"""DVFS reaction time: die-voltage response to reference steps 1.0 -> 1.1 -> 1.0 V while qsort runs, for
(a) off-chip VR only (slew-limited, 20 and 50 mV/us), (b) die fed only by on-die regulating bucks (FIVR style),
(c) hybrid: off-chip path + on-die regulating bucks in parallel. Settle time = first time after the step from which the
cycle-averaged mean die voltage stays within +-1 % of the new reference for 1 us."""
import os, sys, json, numpy as np
os.environ['PDN3_LIB'] = 'libpdn3b.so'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn3
from multiprocessing import Pool
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
T = pdn3.T_CYC; S1, S2, N = 3000, 13000, 23000
def settle(v, ref, k0, k1):
    v = np.convolve(v, np.ones(100) / 100, 'same')          # 300 ns moving average removes workload noise
    tol = 0.01 * ref[k0 + 1]; ok = np.abs(v[k0:k1] - ref[k0 + 1]) < tol; need = int(1e-6 / T)
    for k in range(len(ok) - need):
        if ok[k:k + need].all(): return k * T * 1e9
    return float('nan')
def job(case):
    name, p, spec = case[:3]; Ax = case[3] if len(case) > 3 else 0.4
    I = np.load(B + 'qsort_rep_A.npz')['I'][60000:60000 + N].astype(np.float32)
    vref = np.full(N, 1.0); vref[S1:S2] = 1.1
    m = pdn3.simulate(I, dict(dict(A_x=[Ax * s for s in pdn3.SHARE]), **p), spec, vref=vref, want_v=True, warm=4000)
    v = m['vavg'].mean(1)
    up = settle(v, vref, S1, S2); dn = settle(v, vref, S2, N)
    ovs = (v[S1:S2].max() - 1.1) * 1e3; und = (1.0 - v[S2:].min()) * 1e3
    return name, up, dn, ovs, und, m['p999'], m['loss_buck']
reg = lambda L, fsw, fc, Ron, no_int=0, Ax=0.4: [dict(reg=dict(L=L, fsw=fsw, fc=fc, Ron=Ron, no_int=no_int, Cr=(Ax * 12e-9 + 3e-9) * s)) for s in pdn3.SHARE]
cases = [('off-chip VR only, slew 20 mV/us', dict(slew_vr=20e3), None),
         ('off-chip VR only, slew 50 mV/us', dict(slew_vr=50e3), None),
         ('off-chip VR only, no slew limit (loop 100 kHz)', dict(slew_vr=1e9), None),
         ('on-die FIVR only, 5 nH, 100 MHz, fc 10 MHz', dict(fivr_only=True, slew_vr=20e3), reg(5e-9, 100e6, 10e6, 1.0)),
         ('on-die FIVR only, 5 nH, 100 MHz, fc 3 MHz', dict(fivr_only=True, slew_vr=20e3), reg(5e-9, 100e6, 3e6, 1.0)),
         ('on-die FIVR only, 2 nH, 140 MHz, fc 14 MHz', dict(fivr_only=True, slew_vr=20e3), reg(2e-9, 140e6, 14e6, 1.0)),
         ('on-die FIVR only, 2 nH, 140 MHz, fc 14 MHz, decap 2 mm2', dict(fivr_only=True, slew_vr=20e3), reg(2e-9, 140e6, 14e6, 1.0, 0, 2.0), 2.0),
         ('on-die FIVR only, 1 nH, 500 MHz, fc 50 MHz, decap 1 mm2', dict(fivr_only=True, slew_vr=20e3), reg(1e-9, 500e6, 50e6, 0.5, 0, 1.0), 1.0),
         ('hybrid: off-chip 20 mV/us + on-die P-only bucks fc 10 MHz', dict(slew_vr=20e3), reg(5e-9, 100e6, 10e6, 1.0, 1)),
         ('hybrid: off-chip 20 mV/us + on-die PI bucks fc 10 MHz', dict(slew_vr=20e3), reg(5e-9, 100e6, 10e6, 1.0, 0))]
if __name__ == '__main__':
    with Pool(len(cases)) as P: res = P.map(job, cases)
    print('%-52s %10s %10s %9s %9s %8s %8s' % ('architecture', 'settle up', 'settle dn', 'overshoot', 'undersh.', 'p99.9', 'reg loss'))
    for r in res: print('%-52s %8.0f ns %8.0f ns %7.1f mV %7.1f mV %6.1f mV %7.1f%%' % (r[0], r[1], r[2], r[3], r[4], r[5], 100 * r[6]))
    json.dump(res, open('dvfs3.json', 'w'), indent=1)
