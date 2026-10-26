"""Switching-level injector model (50 ps steps): a synchronous buck from 1.8 V with L = 2 nH and 0.1-ohm switches,
hysteretic current control tracking the predicted reference j_ref (A, per cycle, zero-order hold).
States: HS on (di/dt = (1.8-0.9)/L), LS on (di/dt = -0.9/L), OFF (tri-state; current decays through the LS body diode,
blocked at 0). Idle rule: if |j_ref| < I_idle and |i_L| < I_idle -> OFF (no switching).
Loss = sum i^2 (Ron + DCR) dt + E_sw per HS<->LS/OFF transition (65 pJ per pair -> 32.5 pJ per transition)."""
import numpy as np


def track(jref, L=2e-9, R=0.12, band=0.02, i_idle=0.01, sub=20, E_tr=32.5e-12, Ron=0.1):
    n = len(jref); dt = 1e-9 / sub
    i = 0.0; state = 0; ntr = 0; Ec = 0.0
    out = np.zeros(n)
    for k in range(n):
        r = jref[k]; acc = 0.0
        for q in range(sub):
            if abs(r) < i_idle and abs(i) < i_idle:
                new = 0
            elif i < r - band:
                new = 1
            elif i > r + band:
                new = 2
            else:
                new = state if state != 0 else (1 if i < r else 2)
            if new != state: ntr += 1; state = new
            if state == 1: i += (0.9 - i * Ron) / L * dt            # switch IR drop eats into the 0.9 V headroom
            elif state == 2: i += (-0.9 - i * Ron) / L * dt
            else:
                if i > 0: i = max(0.0, i - 0.9 / L * dt)
                elif i < 0: i = min(0.0, i + 0.9 / L * dt)
            Ec += i * i * R * dt; acc += i
        out[k] = acc / sub
    return out, ntr, Ec, ntr * E_tr
