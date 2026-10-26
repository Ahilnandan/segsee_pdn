"""Charge-pump benefit PER EXECUTION UNIT, not lumped-tile.

A single pooled current trace hides which physical block actually needs a
local reservoir: dhrystone's CSR unit has a 334 mA peak (a single serializing
uop) riding on a 4 mA mean, which is invisible once averaged into a 138 mA-mean
ALU trace. The trigger and the bank sizing both have to be conditioned on
"which unit is this uop, and how much is ALREADY committed on that unit's
local rail" -- exactly per-instance opcode -> fu_code -> physical unit,
already captured in the trace (perunit.py).

Method: per unit, build its own no-bank equivalent-decap curve (droop(I_unit,C)
swept over C, since a per-unit trace has different amplitude than the lumped
one eqdecap.json was built from -- that curve does NOT apply here), then
compare against the SAME validated token-bucket bank math (perunit.apply_bank,
extracted from replay2.run(), TH/Qb/Ir as before) applied to that unit's own
current. The dominant unit is spot-checked with the REAL Pulsim circuit
(pulsim_chargepump.py) as physical validation, not just the closed-form model.
"""
import sys, numpy as np
from sweep import setup
from replay2 import run, droop
from perunit import per_unit_current, apply_bank, UNIT_GROUPS

def eqdecap_curve(I, Cs_nF):
    return np.array([droop(I, C=C * 1e-9).max() for C in Cs_nF])

def eq_decap_for(target_mV, Cs_nF, droops):
    order = np.argsort(Cs_nF)
    Cs, Ds = Cs_nF[order], droops[order]  # Ds decreasing in C
    if target_mV >= Ds.max(): return float('inf')
    if target_mV <= Ds.min(): return float(Cs.min())
    return float(np.interp(target_mV, Ds[::-1], Cs[::-1]))

if __name__ == '__main__':
    wl = sys.argv[1] if len(sys.argv) > 1 else 'dhrystone'
    n_cap = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    TH, QB, IR = 100, 250, 25.0   # same validated sweep point as bank-token-scheduling-result

    inst, prods, cons, age, pin = setup(wl)
    R = run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0)
    per, t0 = per_unit_current(inst, R['tnew'], n=len(R['I']))

    print(f'{wl}: per-execution-unit charge-pump analysis, {n_cap} cycles, '
          f'bank Qb={QB}pC Ir={IR}mA TH={TH}mA\n')
    print(f'{"unit":5s} {"mean mA":>8s} {"peak mA":>8s} {"no-bank eqC(nF)":>16s} '
          f'{"banked eqC(nF)":>15s} {"benefit(nF)":>12s} {"events":>7s} {"served%":>8s}')
    results = {}
    Cs = np.array([0.5, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 80])
    for unit in UNIT_GROUPS:
        I = per[unit][:n_cap]
        if I.max() - I.min() < 1.0:
            continue  # dead unit on this workload, skip
        d_nobank = droop(I)
        Ds = eqdecap_curve(I, Cs)
        eq_nobank = eq_decap_for(d_nobank.max(), Cs, Ds)

        In, stats = apply_bank(I, TH, QB, IR)
        d_banked = droop(In)
        eq_banked = eq_decap_for(d_banked.max(), Cs, Ds)
        served = 100 * stats['evfull'] / stats['events'] if stats['events'] else 100.0

        results[unit] = dict(I=I, d_nobank=d_nobank, d_banked=d_banked,
                              eq_nobank=eq_nobank, eq_banked=eq_banked)
        benefit = (eq_banked - eq_nobank) if np.isfinite(eq_banked) and np.isfinite(eq_nobank) else float('nan')
        print(f'{unit:5s} {I.mean():8.1f} {I.max():8.1f} {eq_nobank:16.1f} '
              f'{eq_banked:15.1f} {benefit:12.1f} {stats["events"]:7d} {served:7.1f}%')

    # physical validation: run the REAL Pulsim circuit on the dominant unit
    dom = max(results, key=lambda u: results[u]['I'].sum())
    print(f'\nSpot-checking the dominant unit ({dom}, {100*results[dom]["I"].sum()/sum(r["I"].sum() for r in results.values()):.0f}% '
          f'of decomposed energy) against the real Pulsim switched-capacitor circuit...')
    from pulsim_chargepump import run_chargepump, eq_decap_nF
    I_dom = results[dom]['I'][:6000]
    d_pump, v_bank = run_chargepump(I_dom, ir_recharge_mA=IR, sub=4)
    d_nobank_dom = droop(I_dom)
    print(f'  closed-form no-bank max droop:  {d_nobank_dom.max():.2f} mV')
    print(f'  real circuit (Cbank=0.5nF) max droop: {d_pump.max():.2f} mV  '
          f'(bank swung {v_bank.min():.0f}..{v_bank.max():.0f} mV)')
    print(f'  closed-form TOKEN-BUCKET (same TH/Qb/Ir) max droop: '
          f'{droop(apply_bank(I_dom, TH, QB, IR)[0]).max():.2f} mV  <- idealized-vs-real agreement check')
