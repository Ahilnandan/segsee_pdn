"""A1: timing-error model for droop-induced datapath errors (circuit-state-aware replay, Phase 5).
All numbers below are ASSUMPTIONS, stated so they can be swept; there are no per-path timing reports for this core.
  delay law   : alpha-power d(V) ~ V / (V - VT)^ALPHA, VT = 0.40 V, ALPHA = 1.3 (45 nm class)
  clocking    : the clock is set so the chip-critical path just meets the period at V = VNOM - G (G = guardband, mV)
  unit paths  : R[u] = unit's longest path / chip-critical path at nominal V
  sensitising : P[u] = probability one execution exercises that longest path (data-dependent)
  detection   : shadow sample at (1 + W) x period; a late path within it is a DETECTED error (replayable),
                beyond it an UNDETECTED error (silent corruption) -> the floor the PDN must hold
  timing      : a uop issued at cycle t executes in cycles t+2 .. t+1+LAT[u]; it fails if the worst-tile voltage in any of
                those cycles is low enough. Datapath units only (control state and SRAM are out of scope)."""
import numpy as np
VT, ALPHA, VNOM, W = 0.40, 1.3, 1.1, 0.25
FU = {1: 'ALU', 2: 'BRU', 4: 'MEM', 8: 'MUL', 16: 'DIV', 32: 'CSR', 64: 'FPU', 128: 'FDV', 256: 'I2F', 512: 'F2I'}
R = dict(ALU=0.85, BRU=0.80, MEM=1.00, MUL=1.00, DIV=0.70, CSR=0.60, FPU=1.00, FDV=0.80, I2F=0.80, F2I=0.80)
P = dict(ALU=0.02, BRU=0.02, MEM=0.05, MUL=0.10, DIV=0.02, CSR=0.02, FPU=0.10, FDV=0.05, I2F=0.05, F2I=0.05)
LAT = dict(ALU=1, BRU=1, MEM=1, MUL=3, DIV=1, CSR=1, FPU=4, FDV=1, I2F=2, F2I=2)
def d(V): return V / (V - VT) ** ALPHA
def v_fail(r, G_mV, w=0.0):
    """Voltage below which a path of relative length r misses (1+w) x period, for guardband G (bisection)."""
    target = d(VNOM - G_mV * 1e-3) * (1 + w) / r; lo, hi = VT + 0.05, VNOM + 0.2
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if d(m) > target: lo = m
        else: hi = m
    return 0.5 * (lo + hi)
def errors(tnew, fu, vworst, G_mV, offset=2, R=R, P=P, LAT=LAT, w=W):
    """Expected detected / undetected error counts per unit for one schedule. vworst: worst-tile V per cycle."""
    out = {}; n = len(vworst)
    for code, u in FU.items():
        idx = np.where(fu == code)[0]
        if not len(idx): continue
        t0 = np.clip(tnew[idx] + offset, 0, n - 1); vmin = vworst[t0].copy()
        for k in range(1, LAT[u]): vmin = np.minimum(vmin, vworst[np.clip(t0 + k, 0, n - 1)])
        vf, vu = v_fail(R[u], G_mV), v_fail(R[u], G_mV, w)
        late = vmin < vf; silent = vmin < vu
        out[u] = dict(n=len(idx), late=int(late.sum()), det=float(P[u] * (late & ~silent).sum()), undet=float(P[u] * silent.sum()))
    return out
