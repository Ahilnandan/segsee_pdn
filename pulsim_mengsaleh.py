"""Meng & Saleh (JSSC 2009) switched-decap active decoupling circuit, built
from the topology validated by looking at Fig 1(b)/(c) directly (not the
OCR'd text) and confirmed by first-principles charge conservation -- see
pdn_work/docs/reorder-strategy-and-chargepump-benefits.md section 6.

Mechanism (charge REDISTRIBUTION, no external recharge energy -- different
from our own bank in pulsim_chargepump.py, which adds new charge from a
current-limited external rail):
  standby (parallel): cap_L (a1-b1) and cap_R (a2-gnd) each independently
    across vdie-gnd -- ordinary 2*Cdecap passive decap.
  boost (series, comparator-triggered on droop): b1 and a2 merge (the
    cross-coupled Mn2/Mp2 switch pair); a1 (cap_L's outer plate) is left
    connected to vdie only through the SAME finite on-resistance its
    parallel-mode switch (Mp1) always has -- this is the one place we
    approximate: the paper models Mp1's "off" state as subthreshold leakage
    via a fitted f(Ron) factor from their own HSPICE (Fig 2), which we don't
    have; we approximate it as one fixed linear R_on (10-20 Ohm range they
    quote in the text) always present between a1 and vdie. By charge
    conservation, a1 (now isolated except through R_on) is being driven
    toward ~2x its pre-boost value by the merge -- the excess above vdie
    drains back into the rail through R_on, delivering more charge than a
    plain decap would over the same droop.

Two decaps totaling 0.5 nF (same silicon budget as our own Cbank in
pulsim_chargepump.py, for an apples-to-apples nF-spent comparison), R_on=15
Ohm (their quoted 10-20 Ohm design range), comparator threshold 50 mV
(their stated ~50 mV design sensitivity, Section II.B).
"""
import sys, numpy as np
import pulsim as p
from replay2 import droop

L_H, C_F, R_OHM = 507e-12, 5e-9, 79.6e-3
CDECAP_F = 0.25e-9         # each of the two decaps; 0.5nF total, matches our own bank's silicon budget
RON = 15.0                  # ohm, mid of their quoted 10-20 Ohm design range
TH_DROOP_V = 0.050          # their stated ~50 mV comparator sensitivity

def build_circuit(i0_A):
    b = p.CircuitBuilder()
    b.add_voltage_source('Vrm', 'vrm', 'gnd', 0.0)
    b.add_inductor('L1', 'vrm', 'vmid', L_H, i0=i0_A)
    b.add_resistor('R1', 'vmid', 'vdie', R_OHM)
    b.add_capacitor('C1', 'vdie', 'gnd', C_F)
    b.set_initial('C1', -R_OHM * i0_A)
    b.add_current_source('Iload', 'vdie', 'gnd', 0.0)   # PWL via b_extra

    # cap_L: a1(top, outer/floating-during-boost) -- C1 -- b1(bottom, merge node)
    b.add_capacitor('Cdecap_L', 'a1', 'b1', CDECAP_F)
    b.set_initial('Cdecap_L', 0.0)   # charged to vdie's nominal (0 deviation) at rest
    # cap_R: a2(top, merge node) -- C2 -- gnd(bottom, always fixed)
    b.add_capacitor('Cdecap_R', 'a2', 'gnd', CDECAP_F)
    b.set_initial('Cdecap_R', 0.0)

    b.add_resistor('Ron_a1', 'a1', 'vdie', RON)          # always-present delivery path (Mp1's finite Ron)
    b.add_switch('SW_b1_gnd', 'b1', 'gnd', g_on=1e4, g_off=1e-9)   # Mn1: parallel-mode only
    b.add_switch('SW_a2_vdie', 'a2', 'vdie', g_on=1e4, g_off=1e-9) # ties cap_R to vdie: parallel-mode only
    b.add_switch('SW_merge', 'b1', 'a2', g_on=1e4, g_off=1e-9)     # Mn2/Mp2: series/boost-mode only
    return b

def run_mengsaleh(I_mA, cyc_s=1e-9, sub=4):
    dt_s = cyc_s / sub
    i0_A = I_mA[0] * 1e-3
    b = build_circuit(i0_A)
    vdie_idx = b.node_id_of('vdie')
    nunk = b.graph.num_branches
    nsw = b.graph.num_switches
    n = len(I_mA)
    I_A = np.asarray(I_mA, float) * 1e-3

    state = {'v_die': 0.0, 'boost': False}
    def step_observer(t, x):
        state['v_die'] = float(x[vdie_idx])
        # simple hysteresis, mirroring their latch-comparator design intent
        if not state['boost'] and state['v_die'] < -TH_DROOP_V:
            state['boost'] = True
        elif state['boost'] and state['v_die'] > -TH_DROOP_V * 0.3:
            state['boost'] = False

    def switch_fn(t):
        m = p.SwitchStateMask(nsw)
        boost = state['boost']
        m.set(0, not boost)  # SW_b1_gnd: ON in parallel mode
        m.set(1, not boost)  # SW_a2_vdie: ON in parallel mode
        m.set(2, boost)      # SW_merge: ON in series/boost mode
        return m

    def b_extra(t):
        k = int(t / cyc_s + 1e-9)
        if k >= n: k = n - 1
        out = [0.0] * nunk
        out[vdie_idx] = I_A[k]
        return out

    res = p.simulate(b, t_end=n * cyc_s, dt=dt_s, b_extra_fn=b_extra, switch_fn=switch_fn,
                      step_observer=step_observer, start_from_dc_op=False)
    states = np.asarray(res.states)
    vdie = states[::sub][:n, vdie_idx]
    return -vdie * 1e3  # droop mV

def eq_decap_nF(wl, target_mV, eqdecap_json='eqdecap.json'):
    import json
    rows = [r for r in json.load(open(eqdecap_json)) if r['wl'] == wl]
    rows.sort(key=lambda r: r['C'])
    Cs = np.array([r['C'] for r in rows]); Ms = np.array([r['max'] for r in rows])
    if target_mV >= Ms.max(): return float('inf')
    if target_mV <= Ms.min(): return float(Cs.min())
    return float(np.interp(target_mV, Ms[::-1], Cs[::-1]))

if __name__ == '__main__':
    wl = sys.argv[1] if len(sys.argv) > 1 else 'dhrystone'
    n_cap = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    unit = sys.argv[3] if len(sys.argv) > 3 else None

    if unit:
        from sweep import setup
        from replay2 import run as run_sched
        from perunit import per_unit_current
        inst, prods, cons, age, pin = setup(wl)
        R = run_sched(inst, prods, cons, age, pin, 'A', 1e9, 0, 0)
        per, t0 = per_unit_current(inst, R['tnew'], n=len(R['I']))
        I = per[unit][:n_cap]
        label = f'{wl}/{unit}'
    else:
        from sweep import setup
        from replay2 import run as run_sched
        inst, prods, cons, age, pin = setup(wl)
        R = run_sched(inst, prods, cons, age, pin, 'A', 1e9, 0, 0)
        I = R['I'][:n_cap]
        label = wl

    print(f'{label}: {len(I)} cycles, min={I.min():.1f} max={I.max():.1f} mA')
    d_nobank = droop(I)
    print(f'no active help (actual 5nF):        max droop = {d_nobank.max():6.2f} mV')

    d_meng = run_mengsaleh(I, sub=4)
    print(f'Meng & Saleh switched-decap (0.5nF): max droop = {d_meng.max():6.2f} mV')

    eq_nobank = eq_decap_nF(wl, d_nobank.max())
    eq_meng = eq_decap_nF(wl, d_meng.max())
    print(f'Equivalent passive decap: no-help {eq_nobank:.1f} nF  ->  Meng&Saleh circuit {eq_meng:.1f} nF')
    if np.isfinite(eq_meng) and np.isfinite(eq_nobank):
        benefit = eq_meng - eq_nobank
        print(f'Benefit: {benefit:.2f} nF for 0.5 nF of physical switched decap '
              f'-- {benefit/0.5:.1f}x multiplier')

    np.savez(f'pulsim_mengsaleh_{label.replace("/","_")}.npz', I=I, d_nobank=d_nobank, d_meng=d_meng)
