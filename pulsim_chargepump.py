"""The REAL switched-charge-pump / active-decap circuit, simulated (not just
token-bucket-approximated), cross-checked against the raw baseline and the
idealized python token-bucket bank, using the same equivalent-passive-decap
metric as eqdecap.py.

What existed before this file: a PYTHON TOKEN-BUCKET ABSTRACTION (replay2.py's
Qb/Ir bank inside run()) and a mirrored RTL ACCOUNTING REGISTER (core.scala's
BACT_TOK) -- both just track a charge NUMBER, assumed to deliver current
instantly and without its own voltage sagging. Neither is a circuit. This file
is the first actual circuit: a real capacitor that sags as it's drawn down,
a real diode that self-gates release, and a recharge path that's current-
limited and self-shuts-off -- built on the Pulsim engine already cross-
validated against replay2.droop() (see pulsim-pdn-cross-validation memory).

Topology:
  Vrm(0V) --L(507pH)--R(79.6mOhm)-- vdie --C(5nF)-- gnd      [core PDN, as before]
  vdie <--SW_release(Ron)-- bank                               [controlled switch:
                                                                 the "current-controlled
                                                                 release PMOS" from the
                                                                 architecture note]
  bank --Cbank(0.5nF)-- gnd, precharged to +V_HEAD
  recharge: current-limited injection into `bank`, gated OFF once bank reaches
  +V_HEAD again (python-side control via step_observer, mirroring a real
  current-limited charger IC's regulation loop)

FIRST ATTEMPT (kept as a note, not code): a plain diode instead of a controlled
switch for release. Wrong -- an always-forward-biased ideal/near-ideal diode
between a reservoir held ABOVE the rail and the rail itself has nothing to hold
it up: Cbank*Ron ~ 0.05ns means it instantly bleeds down to the rail and stays
there, turning "reservoir with triggered burst release" into "small always-on
current source", which is not the architecture and not what a passive decap
comparison should be measured against. Switched to a python-side comparator-
gated switch (release ON only while vdie droops past TH_RELEASE_MV AND bank
still has headroom over vdie) -- this is what actually reproduces "sits idle,
releases current in a burst that lowers total average draw at TH".

Cb=0.5nF / V_HEAD=0.25V (0.5V full swing) => 250pC usable charge, matching the
architecture note and the Qb=250 pC token-bucket runs. Ir default 25 mA,
matching the best token-bucket sweep point.
"""
import sys, numpy as np
import pulsim as p
from sweep import setup
from replay2 import run, droop

L_H, C_F, R_OHM = 507e-12, 5e-9, 79.6e-3
CB_F = 0.5e-9
V_HEAD = 0.25            # bank precharge above nominal Vdd (0V ref) -> 0.5V full swing -> 250pC
RON_RELEASE = 0.05       # release-FET on-resistance, ohm, while the comparator has it ON
TH_RELEASE_V = 0.020     # comparator threshold: engage release once vdie droops past this

def build_circuit(i0_A):
    b = p.CircuitBuilder()
    b.add_voltage_source('Vrm', 'vrm', 'gnd', 0.0)
    b.add_inductor('L1', 'vrm', 'vmid', L_H, i0=i0_A)
    b.add_resistor('R1', 'vmid', 'vdie', R_OHM)
    b.add_capacitor('C1', 'vdie', 'gnd', C_F)
    b.set_initial('C1', -R_OHM * i0_A)
    b.add_current_source('Iload', 'vdie', 'gnd', 0.0)       # PWL via b_extra

    b.add_capacitor('Cbank', 'bank', 'gnd', CB_F)
    b.set_initial('Cbank', V_HEAD)
    b.add_switch('SW_release', 'bank', 'vdie', g_on=1.0 / RON_RELEASE, g_off=1e-9)  # comparator-driven
    b.add_current_source('Irechg', 'bank', 'gnd', 0.0)      # PWL via b_extra, gated in step_observer
    return b

def run_chargepump(I_mA, ir_recharge_mA=25.0, cyc_s=1e-9, sub=4):
    dt_s = cyc_s / sub
    i0_A = I_mA[0] * 1e-3
    ir_A = ir_recharge_mA * 1e-3
    b = build_circuit(i0_A)
    vdie_idx = b.node_id_of('vdie')
    bank_idx = b.node_id_of('bank')
    nunk = b.graph.num_branches
    nsw = b.graph.num_switches
    n = len(I_mA)
    I_A = np.asarray(I_mA, float) * 1e-3

    state = {'v_bank': V_HEAD, 'v_die': 0.0}
    def step_observer(t, x):
        state['v_bank'] = float(x[bank_idx])
        state['v_die'] = float(x[vdie_idx])

    def switch_fn(t):
        m = p.SwitchStateMask(nsw)
        # release ON only while actually drooping past threshold AND bank still
        # has headroom over the rail (never let current flow backwards)
        release_on = (state['v_die'] < -TH_RELEASE_V) and (state['v_bank'] > state['v_die'])
        m.set(0, release_on)
        return m

    def b_extra(t):
        k = int(t / cyc_s + 1e-9)
        if k >= n: k = n - 1
        out = [0.0] * nunk
        out[vdie_idx] = I_A[k]                                # load draw
        recharging = state['v_bank'] < V_HEAD
        out[bank_idx] = -ir_A if recharging else 0.0          # current-limited recharge, self-gated
        return out

    res = p.simulate(b, t_end=n * cyc_s, dt=dt_s, b_extra_fn=b_extra, switch_fn=switch_fn,
                      step_observer=step_observer, start_from_dc_op=False)
    states = np.asarray(res.states)
    vdie = states[::sub][:n, vdie_idx]
    vbank = states[::sub][:n, bank_idx]
    return -vdie * 1e3, vbank * 1e3   # droop mV, bank voltage mV

def eq_decap_nF(wl, target_mV, eqdecap_json='eqdecap.json'):
    import json
    rows = [r for r in json.load(open(eqdecap_json)) if r['wl'] == wl]
    rows.sort(key=lambda r: r['C'])
    Cs = np.array([r['C'] for r in rows]); Ms = np.array([r['max'] for r in rows])
    if target_mV >= Ms.max():
        return float('inf')  # even the largest tested C doesn't reach this low
    if target_mV <= Ms.min():
        return Cs.min()      # better than the smallest tested C -- extrapolate down cautiously
    # Ms is decreasing in C; interpolate C at the target droop
    return float(np.interp(target_mV, Ms[::-1], Cs[::-1]))

if __name__ == '__main__':
    wl = sys.argv[1] if len(sys.argv) > 1 else 'dhrystone'
    n_cap = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    ir_mA = float(sys.argv[3]) if len(sys.argv) > 3 else 25.0

    inst, prods, cons, age, pin = setup(wl)
    R = run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0)   # raw HW-scheduled current, no bank, no reorder
    I_raw = R['I'][:n_cap].copy()
    print(f'{wl}: {len(I_raw)} cycles raw current, min={I_raw.min():.1f} max={I_raw.max():.1f} mA')

    d_nobank = droop(I_raw)
    print(f'NO bank (actual 5nF):        max droop = {d_nobank.max():6.2f} mV')

    d_pump, v_bank = run_chargepump(I_raw, ir_recharge_mA=ir_mA, sub=4)
    print(f'REAL charge-pump circuit:    max droop = {d_pump.max():6.2f} mV   '
          f'(bank voltage range {v_bank.min():.1f}..{v_bank.max():.1f} mV, precharge {1e3*V_HEAD:.0f} mV)')

    eq_nobank = eq_decap_nF(wl, d_nobank.max())
    eq_pump = eq_decap_nF(wl, d_pump.max())
    print(f'Equivalent passive decap:    no-bank {eq_nobank:.1f} nF  ->  charge-pump {eq_pump:.1f} nF '
          f'(actual physical Cbank = {CB_F*1e9:.2f} nF)')
    if np.isfinite(eq_pump) and np.isfinite(eq_nobank):
        print(f'Benefit: {eq_pump - eq_nobank:.1f} nF of EFFECTIVE decap added by a {CB_F*1e9:.2f} nF '
              f'physical bank -- a {(eq_pump-eq_nobank)/(CB_F*1e9):.1f}x multiplier vs adding {CB_F*1e9:.2f} nF '
              f'of plain passive decap in the same area.')

    np.savez(f'pulsim_chargepump_{wl}.npz', I_raw=I_raw, d_nobank=d_nobank, d_pump=d_pump, v_bank=v_bank)
