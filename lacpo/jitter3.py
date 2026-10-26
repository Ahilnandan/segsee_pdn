"""Supply-induced clock jitter and timing guardband from simulated die voltages.

Clock tree: insertion delay T_ins (0.5 ns, swept), half in a global trunk that sees the area-weighted mean die voltage and
half in local buffers that see their own region's voltage. Buffer delay vs V: alpha-power law d(V) ~ V / (V - Vt)^alpha
(Vt = 0.45 V, alpha = 1.3, Nangate-45-class), normalised to d(VNOM) = 1. Edge n arrives at region r at
a_r[n] = T_ins * (0.5 f(Vmean[n]) + 0.5 f(V_r[n])), using the cycle-averaged voltage the edge propagates through.
Period jitter J_p = a[n+1] - a[n]; cycle-to-cycle J_cc = J_p[n] - J_p[n-1].
Timing guardband: a register-to-register path of length T_path = 0.8 T_clk (launch and capture in the same region, the
capture edge one period later) needs T_path * f(Vmin_r[n]) <= T_clk + (a_r[n+1] - a_r[n]); the guardband is the worst
(p99.99) extra delay as % of the period that the clock must absorb -> frequency headroom lost to supply noise."""
from __future__ import print_function
import numpy as np, pdn3
VT, ALPHA = 0.45, 1.3
def f(V):
    V = np.asarray(V, float); d = V / (V - VT) ** ALPHA; return d / (pdn3.VNOM / (pdn3.VNOM - VT) ** ALPHA)

def jitter(vavg, vmin, T_ins=0.5e-9, T_clk=pdn3.T_CYC, path=0.8):
    vm = (vavg * pdn3.SHARE[None, :]).sum(1)
    a = T_ins * (0.5 * f(vm)[:, None] + 0.5 * f(vavg))                     # arrival per region (s)
    Jp = np.diff(a, axis=0) * 1e12                                          # ps
    Jcc = np.diff(Jp, axis=0)
    need = path * T_clk * f(vmin[:-1]) - (T_clk + np.diff(a, axis=0))       # >0 means the path misses timing
    gb = 100 * (path * T_clk * f(vmin[:-1]) - path * T_clk - np.diff(a, axis=0)) / T_clk
    return dict(Jp_rms=float(Jp.std()), Jp_p999=float(np.percentile(np.abs(Jp), 99.9)), Jp_max=float(np.abs(Jp).max()),
                Jcc_rms=float(Jcc.std()), Jcc_p999=float(np.percentile(np.abs(Jcc), 99.9)),
                guardband_p9999=float(np.percentile(gb.max(1), 99.99)), guardband_max=float(gb.max()))

if __name__ == '__main__':
    import sys
    B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
    I = np.load(B + 'qsort_rep_A.npz')['I'][60000:120000].astype(np.float32)
    for Ax in (0.2, 0.4, 0.8):
        m = pdn3.simulate(I, dict(A_x=[Ax * s for s in pdn3.SHARE]), want_v=True)
        for T_ins in (0.3e-9, 0.5e-9, 1.0e-9):
            j = jitter(m['vavg'], m['vmin'], T_ins)
            print('decap %.1f mm2 p99.9 droop %5.1f mV | T_ins %.1f ns: period jitter rms %.2f ps p99.9 %.2f ps max %.2f | c2c p99.9 %.2f ps | guardband p99.99 %.2f %% max %.2f %%'
                  % (Ax, m['p999'], T_ins * 1e9, j['Jp_rms'], j['Jp_p999'], j['Jp_max'], j['Jcc_p999'], j['guardband_p9999'], j['guardband_max']))
