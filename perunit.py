"""Per-EXECUTION-UNIT current decomposition, from the real fu_code each uop
instance issued on (not a synthetic uopc-class grouping). fu_code is BOOM's
own opcode->hardware-unit mapping (decode.scala / execution-units.scala),
already captured per instance in BACT_CAND and validated in
replay-must-use-real-ports-and-deps.

I[n] (lumped tile current) sums uniformly over every issued uop regardless of
which physical block executed it. That's wrong for anything asking "does THIS
unit's local rail need its own reservoir" -- an ALU burst and a MUL burst are
different physical nodes on the die with different average duty cycles, so a
single shared pump sized off the lumped total either over- or under-provisions
each one. This module attributes each instance's SHAPE current to its actual
fu_code instead of summing everything into one trace.
"""
import numpy as np
from curmodel import SHAPE, FLOOR

FU_NAMES = {1: 'ALU', 2: 'JMP', 4: 'MEM', 8: 'MUL', 16: 'DIV',
            32: 'CSR', 64: 'FPU', 128: 'FDV', 256: 'I2F', 512: 'F2I'}

# physical grouping: ALU+JMP+I2F share the int ALU datapath across ports p0-p3
# (execution-units.scala); CSR is its own pipe on p1; MUL/DIV are dedicated
# units on p2/p3; MEM is the shared AGU/LSU; FPU/FDV/F2I are the FP side.
UNIT_GROUPS = {
    'ALU': (1, 2, 256),      # FU_ALU | FU_JMP | FU_I2F
    'CSR': (32,),            # FU_CSR
    'MUL': (8,),             # FU_MUL
    'DIV': (16,),            # FU_DIV
    'MEM': (4,),             # FU_MEM
    'FPU': (64, 512),        # FU_FPU | FU_F2I
    'FDV': (128,),           # FU_FDV
}
FU_TO_UNIT = {fu: name for name, fus in UNIT_GROUPS.items() for fu in fus}

def per_unit_current(inst, tnew, n=None):
    """inst: candparse.parse() instances (iq,u,fu,rob,pd,prs,r,t).
    tnew: per-instance issue cycle (use inst[i][7]==treal for the unscheduled
    baseline, or run()'s R['tnew'] for a scheduled policy -- both share t0
    normalization with run(), so pass whichever this is being compared against).
    Returns {unit_name: I(t) array}, each on the SAME timeline/floor split
    UNIT_GROUPS.keys() ways so they sum back to the lumped current_from_times().
    """
    t0 = int(min(tnew))
    n = n or (int(max(tnew)) - t0 + 20)
    nunits = len(UNIT_GROUPS)
    floor_each = FLOOR / nunits   # split the always-on floor evenly across units
    out = {name: np.full(n, floor_each) for name in UNIT_GROUPS}
    for i in range(len(inst)):
        fu = inst[i][2]
        unit = FU_TO_UNIT.get(fu)
        if unit is None:
            continue  # shouldn't happen for real BOOM fu_code, but don't crash silently elsewhere
        u = inst[i][1]
        k = int(tnew[i]) - t0
        s = SHAPE.get(u)
        if s is None or k < 0:
            continue
        end = min(k + len(s), n)
        if end > k:
            out[unit][k:end] += s[:end - k]
    return out, t0

def apply_bank(I, TH, Qb, Ir, W=5):
    """Standalone extraction of replay2.run()'s bank/slew-clamp math (same
    formula, validated there), applied directly to an arbitrary current
    array -- lets us bank-assist a PER-UNIT current trace without re-running
    the full instruction-scheduling replay for each unit."""
    n = len(I); In = np.array(I, float); tok = float(Qb)
    ev = evfull = 0; unc = drawn = 0.0; inev = False; short = False
    for t in range(n):
        tok = min(float(Qb), tok + Ir)
        lo = In[max(0, t - W):t].min() if t > 0 else I[0]
        ex = max(0.0, I[t] - (lo + TH)); b = 0.0
        if ex > 0:
            b = min(ex, tok); tok -= b; drawn += b; unc += ex - b
            if not inev: ev += 1; inev = True; short = False
            if b < ex - 1e-9: short = True
        elif inev:
            evfull += (not short); inev = False
        In[t] = I[t] - b
    return In, dict(events=ev, evfull=evfull, unc_nC=unc / 1e3, drawn_nC=drawn / 1e3)

def verify_conservation(inst, tnew, I_lumped, t0_lumped):
    """Sanity check: sum of per-unit arrays must equal the lumped current
    array exactly (same instances, same shapes, just re-bucketed)."""
    per, t0 = per_unit_current(inst, tnew, n=len(I_lumped))
    assert t0 == t0_lumped, f'timeline mismatch: {t0} vs {t0_lumped}'
    total = sum(per.values())
    err = np.abs(total - I_lumped)
    return float(err.max()), float(err.mean())

if __name__ == '__main__':
    import sys
    from sweep import setup
    from replay2 import run
    from curmodel import current_from_times
    wl = sys.argv[1] if len(sys.argv) > 1 else 'dhrystone'
    inst, prods, cons, age, pin = setup(wl)
    R = run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0)
    per, t0 = per_unit_current(inst, R['tnew'], n=len(R['I']))
    maxerr, meanerr = verify_conservation(inst, R['tnew'], R['I'], t0)
    print(f'{wl}: conservation check max|err|={maxerr:.6f} mA mean|err|={meanerr:.6f} mA '
          f'(should be ~0 -- per-unit split must sum back to the lumped current)')
    for name, I in per.items():
        print(f'  {name:4s}: mean={I.mean():7.2f} mA  peak={I.max():7.2f} mA  '
              f'share of total energy={100*I.sum()/R["I"].sum():5.1f}%')
