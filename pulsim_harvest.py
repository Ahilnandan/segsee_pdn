"""Circuit-level check of the EDGE-HARVESTING bank (replay2 recyc) in Pulsim.

Same circuit as pulsim_chargepump.py (deviation-referenced PDN, 0.5 nF reservoir precharged
+0.25 V above nominal, comparator-gated release switch, current-limited self-gated recharge
from an external rail), plus a harvest path:
  on a FALLING edge the scheduler already knows about (load current it issued drops more than
  THF below its own recent 5-cycle max -- the same causal rule as the replay), the bank draws
  sink = (max_recent - THF) - I_load from vdie (capped at SMAX) and pushes it into the
  reservoir through a boost transfer: I_bank = eta * (0.9+v_die)/(0.9+v_bank) * sink
  (energy-correct, eta = 0.85), only while the reservoir is below its full precharge.
Harvest on/off is compared on the SAME load trace, so only the circuit physics differs.
"""
import sys, numpy as np
import pulsim as p
from pulsim_chargepump import build_circuit, V_HEAD, TH_RELEASE_V
from replay2 import droop

def run(I_mA, ir_mA=25.0, harvest=False, THF=150.0, eta=0.85, SMAX=150.0, cyc_s=1e-9, sub=4):
    dt = cyc_s / sub
    I_A = np.asarray(I_mA, float) * 1e-3; n = len(I_A)
    b = build_circuit(I_A[0])
    vid, bid = b.node_id_of('vdie'), b.node_id_of('bank')
    nunk, nsw = b.graph.num_branches, b.graph.num_switches
    st = {'vd': 0.0, 'vb': V_HEAD, 'hq': 0.0}
    DR = np.full(n, -9.0)   # per-cycle worst droop (V, positive = below nominal)
    rmax = np.array([I_mA[max(0, k - 5):k].max() if k > 0 else I_mA[0] for k in range(n)])   # causal recent max
    def obs(t, x):
        st['vd'] = float(x[vid]); st['vb'] = float(x[bid])
        k = int(t / cyc_s + 1e-9)
        if k < n and -st['vd'] > DR[k]: DR[k] = -st['vd']
    def sw(t):
        m = p.SwitchStateMask(nsw)
        m.set(0, (st['vd'] < -TH_RELEASE_V) and (st['vb'] > st['vd']))
        return m
    def bex(t):
        k = min(n - 1, int(t / cyc_s + 1e-9))
        out = [0.0] * nunk
        out[vid] = I_A[k]
        rech = (ir_mA * 1e-3) if st['vb'] < V_HEAD else 0.0
        out[bid] = -rech
        if harvest and st['vb'] < V_HEAD:
            sink = min(SMAX, max(0.0, rmax[k] - THF - I_mA[k])) * 1e-3
            if sink > 0:
                gain = eta * (0.9 + st['vd']) / (0.9 + st['vb'])
                out[vid] += sink; out[bid] -= gain * sink; st['hq'] += gain * sink * dt
        return out
    p.simulate(b, t_end=n * cyc_s, dt=dt, b_extra_fn=bex, switch_fn=sw, step_observer=obs, start_from_dc_op=False, store_every=100000)
    return DR * 1e3, st['hq'] * 1e12   # per-cycle worst droop (mV), harvested charge (pC)

if __name__ == '__main__':
    from sweep import setup
    from replay2 import run as replay
    wl = sys.argv[1] if len(sys.argv) > 1 else 'qsort'
    inst, prods, cons, age, pin = setup(wl)
    R = replay(inst, prods, cons, age, pin, 'A', 1e9, 0, 0)
    I = R['I']
    for ir in (25, 10):
        d0, _ = run(I, ir_mA=ir, harvest=False)
        d1, hq = run(I, ir_mA=ir, harvest=True)
        print(f'{wl} raw HW current, real bank circuit Ir={ir}mA: worst droop no-harvest {d0.max():.1f} mV -> harvest {d1.max():.1f} mV; '
              f'p99.9 {np.percentile(d0,99.9):.1f} -> {np.percentile(d1,99.9):.1f}; harvested {hq/1e3:.0f} nC', flush=True)
