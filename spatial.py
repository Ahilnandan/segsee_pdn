"""Spatial port steering: BOOM's greedy grant puts 51-60% of integer EX charge on ALU port 0.
Re-assigning ALU ops among the four identical ALUs keeps the SAME issued set at the SAME cycles
(zero performance cost) and only moves where the current is drawn.

Spatial PDN (linear, exact ZOH like replay2.droop):
  Vrm(0) -- L_pkg -- R_pkg -- common node (C_c) --R_g-- island_w (C_l), w = 0..3 ALU ports
  each int op draws f_loc of its current at its ALU island, the rest at the common node;
  MEM/FP ops draw at the common node. Deviation coordinates (droop positive = below nominal).
Unknowns for this layout (R_g, f_loc, C_l) are swept, not assumed.
Policies: greedy (BOOM: age order, lowest free eligible port), rr (rotating start port),
          steer (most-constrained ops first, then heaviest ALU op -> island with highest voltage now).
"""
import numpy as np
from scipy.linalg import expm
from replay2 import PORTS
from curmodel import SHAPE, FLOOR

def build(Lp=507e-12, Rp=79.6e-3, Ctot=5e-9, Cl=0.25e-9, Rg=0.2, dt=0.25e-9):
    Cc = Ctot - 4 * Cl
    # states: iL, vc, v0..v3 ; inputs: Ic (common load), I0..I3 (island loads)
    n = 6; A = np.zeros((n, n)); B = np.zeros((n, 5))
    A[0, 0] = -Rp / Lp; A[0, 1] = -1 / Lp
    A[1, 0] = 1 / Cc; B[1, 0] = -1 / Cc
    for w in range(4):
        A[1, 1] -= 1 / (Rg * Cc); A[1, 2 + w] += 1 / (Rg * Cc)
        A[2 + w, 1] += 1 / (Rg * Cl); A[2 + w, 2 + w] -= 1 / (Rg * Cl); B[2 + w, 1 + w] = -1 / Cl
    M = np.zeros((n + 5, n + 5)); M[:n, :n] = A * dt; M[:n, n:] = B * dt
    E = expm(M); return E[:n, :n], E[:n, n:]

def run(inst, tnew, policy='greedy', f_loc=0.3, Rg=0.2, Cl=0.25e-9, sub=4):
    Ad, Bd = build(Rg=Rg, Cl=Cl, dt=1e-9 / sub)
    t0 = int(tnew.min()); T = int(tnew.max()) - t0 + 20
    by = {}
    for i, t in enumerate(tnew):
        by.setdefault(int(t) - t0, []).append(i)
    Icom = np.full(T, FLOOR); Iisl = np.zeros((4, T))
    x = np.zeros(6); i0 = FLOOR * 1e-3; x[0] = i0; x[1:] = -79.6e-3 * i0
    worst = np.zeros(T); rr = 0
    for t in range(T):
        ops = by.get(t, [])
        ints = [i for i in ops if inst[i][0] == 1]
        for i in ops:
            if inst[i][0] != 1:
                s = SHAPE.get(inst[i][1], SHAPE[0]); m = min(len(s), T - t); Icom[t:t + m] += s[:m]
        used = [False] * 4; assign = {}
        if policy == 'greedy':
            ints.sort(key=lambda i: inst[i][3])
            for i in ints:
                for w in range(4):
                    if not used[w] and inst[i][2] & PORTS[1][w]: used[w] = True; assign[i] = w; break
        else:
            ints.sort(key=lambda i: (sum(1 for w in range(4) if inst[i][2] & PORTS[1][w]), -SHAPE.get(inst[i][1], SHAPE[0]).sum()))
            if policy == 'rr': pref = [(rr + k) % 4 for k in range(4)]; rr = (rr + 1) % 4
            else: pref = list(np.argsort(-x[2:6]))            # steer: highest island voltage first
            for i in ints:
                for w in pref:
                    if not used[w] and inst[i][2] & PORTS[1][w]: used[w] = True; assign[i] = w; break
        for i, w in assign.items():
            s = SHAPE.get(inst[i][1], SHAPE[0]); m = min(len(s), T - t)
            Iisl[w, t:t + m] += f_loc * s[:m]; Icom[t:t + m] += (1 - f_loc) * s[:m]
        u = np.concatenate([[Icom[t]], Iisl[:, t]]) * 1e-3
        mn = 0.0
        for _ in range(sub):
            x = Ad @ x + Bd @ u; mn = min(mn, x[1:].min())
        worst[t] = -mn * 1e3
    return worst

def job(a):
    wl, pol, fl, rg = a
    from sweep import setup
    from replay2 import run as replay
    inst, prods, cons, age, pin = setup(wl)
    R = replay(inst, prods, cons, age, pin, 'A', 1e9, 0, 0)
    w = run(inst, R['tnew'], pol, fl, rg)
    return dict(wl=wl, pol=pol, f_loc=fl, Rg=rg, droop=float(w.max()), p999=float(np.percentile(w, 99.9)))

if __name__ == '__main__':
    import sys, json
    from multiprocessing import Pool
    jobs = [(wl, pol, fl, rg) for wl in ('qsort', 'dhrystone', 'mm') for fl in (0.3, 0.6) for rg in (0.05, 0.2, 0.5) for pol in ('greedy', 'rr', 'steer')]
    with Pool(6) as p, open('spatial.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
