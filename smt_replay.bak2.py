"""SMT-4 issue replay on a Fu et al. (ICPP'08, Table 2)-class machine, built from real BOOM traces.

Why: on single-thread PaperMatchBoomV3Config only 9-47% of cycles have a real issue choice (mean
choice set 2-3.5), so reordering had nothing to work with. Fu et al. report a ready queue of ~26
(max 73), <9 ready on only 10% of cycles, on an SMT-4 8-wide machine. This model asks: given THAT
much choice, how much does WORK-CONSERVING reordering (same number of ops issued every cycle as
oldest-first -> zero throughput cost) reduce droop and active-decap demand?

Machine (Fu Table 2): 8-wide issue, unified 96-entry IQ shared by 4 threads, ICOUNT dispatch,
8 I-ALU, 4 I-MUL/DIV, 4 LD/ST, 12 FP units.  Dispatch width 8/cycle.

Threads: contiguous program-order slices of the real BOOM traces (qsort/dhrystone/mm), each with
its real dependency graph, real ready times and per-uop current shapes (curmodel.py).
Timing model per instruction i of thread j:
  d_i  = (min over younger-or-equal uops of r) - LEAD   dispatch proxy, monotone in program order
  a_i  = admission into the IQ: program order per thread, >= d_i + slip_j, IQ capacity, 8/cycle
  slip_j  = how far thread j has fallen behind its solo timeline (monotone)
  ready_i = r_i + max(producer delay, a_i - d_i)        (delays only; never earlier than solo)
Producers outside a slice are treated as satisfied.

Policies (all issue EXACTLY the oldest-first count each cycle, per-class unit caps respected):
  A     oldest-first (admission order, then age)
  T     minimise the charge an active bank must supply this cycle (slew clamp lo+TH over W)
  FLAT  minimise squared deviation of projected current from its slow running mean (EMA)
A uop passed over MAXD times is forced (no starvation)."""
import sys, json, numpy as np
from collections import defaultdict
from curmodel import SHAPE, FLOOR
from replay2 import FU_ALU, FU_JMP, FU_MEM, FU_MUL, FU_DIV, FU_CSR, FU_FPU, FU_FDV, FU_I2F, FU_F2I, droop

H = 10; W = 5
CAPS = {'ialu': 8, 'imul': 4, 'ls': 4, 'fp': 12}; ISSUE_W = 8; DISP_W = 8; IQ = 96

def fclass(fu):
    if fu & FU_MEM: return 'ls'
    if fu & (FU_MUL | FU_DIV): return 'imul'
    if fu & (FU_FPU | FU_FDV | FU_F2I): return 'fp'
    return 'ialu'

def make_thread(wl, start, n):
    from sweep import setup
    inst, prods, cons, age, pin = setup(wl)
    order = np.argsort(age, kind='stable')[start:start + n]
    keep = {int(g): k for k, g in enumerate(order)}
    r = np.array([inst[g][6] for g in order], float); tr = np.array([inst[g][7] for g in order], float)
    r = np.where(pin[order], tr, r)
    t0 = r.min(); r -= t0; tr -= t0
    # a real producer is always OLDER in program order; preg-reuse matching in the trace yields some
    # younger 'producers' (false deps) that can deadlock a shared IQ -> drop them
    P = [[keep[p] for p in set(prods[g]) if p in keep and keep[p] < k] for k, g in enumerate(order)]
    sh = [SHAPE.get(inst[g][1], SHAPE[0]) for g in order]
    cl = [fclass(inst[g][2]) for g in order]
    # monotone dispatch proxy: suffix minimum of r in program order
    d = np.minimum.accumulate(r[::-1])[::-1]
    return dict(wl=wl, r=r, tr=tr, P=P, sh=sh, cl=cl, d=d, n=n)

def run(threads, pol='A', LEAD=2, MAXD=8, TH=100.0, ema_tau=64.0, iq=IQ, Qb=None, Ir=25.0, record_ready=True, res=0.0):
    # flatten
    off = np.cumsum([0] + [th['n'] for th in threads]); N = int(off[-1])
    tid = np.concatenate([np.full(th['n'], j) for j, th in enumerate(threads)])
    r = np.concatenate([th['r'] for th in threads]); d = np.concatenate([th['d'] for th in threads]) - LEAD
    TR = np.concatenate([th['tr'] for th in threads])
    sh = [s for th in threads for s in th['sh']]; cl = [c for th in threads for c in th['cl']]
    SH = np.zeros((N, H))
    for i, s in enumerate(sh): m = min(len(s), H); SH[i, :m] = s[:m]
    csum = SH.sum(1)
    prods = [[p + off[j] for p in th['P'][k]] for j, th in enumerate(threads) for k in range(th['n'])]
    cons = [[] for _ in range(N)]
    for i, ps in enumerate(prods):
        for p in ps: cons[p].append(i)
    pend = np.array([len(ps) for ps in prods]); dmax = np.zeros(N)
    nxt = [int(off[j]) for j in range(len(threads))]; slip = [0.0] * len(threads)
    icount = [0] * len(threads); inq = 0
    adm = np.full(N, -1.0); tnew = np.full(N, -1, np.int64); held = np.zeros(N, np.int64)
    wake = defaultdict(list); ready = []
    L = int(max(th['tr'].max() for th in threads) * 8) + 10000
    I = np.full(L, FLOOR); In = np.full(L, FLOOR)
    tok = float(Qb) if Qb is not None else 0.0; drawn = 0.0; unc = 0.0
    ema = None; nshort = 0; SOCt = np.zeros(L); nready = []; nchg = 0; done = 0; t = 0
    def make_ready(i, t):
        rt = r[i] + max(dmax[i], adm[i] - d[i])
        wake[int(max(np.ceil(rt), t + 1))].append(i)
    while done < N:
        if t >= L - 20: raise RuntimeError('stall')
        # --- dispatch (ICOUNT: fewest IQ entries first), program order per thread
        nd = 0
        while nd < DISP_W and inq < iq:
            cand = [j for j in range(len(threads)) if nxt[j] < off[j + 1] and d[nxt[j]] + slip[j] <= t]
            if not cand: break
            j = min(cand, key=lambda j: icount[j]); i = nxt[j]; nxt[j] += 1
            slip[j] = max(slip[j], t - d[i]); adm[i] = t; inq += 1; icount[j] += 1; nd += 1
            if pend[i] == 0: make_ready(i, t)
        if t in wake: ready.extend(wake.pop(t))
        if record_ready: nready.append(len(ready))
        pick = []
        if ready:
            ready.sort(key=lambda i: (adm[i], i))
            cnt = defaultdict(int); base = []
            for i in ready:
                if len(base) >= ISSUE_W: break
                if cnt[cl[i]] < CAPS[cl[i]]: cnt[cl[i]] += 1; base.append(i)
            pick = base
            if pol != 'A' and len(ready) > len(base):
                forced = [i for i in ready if held[i] >= MAXD]
                sel = []; cnt = defaultdict(int)
                for i in forced:
                    if len(sel) < len(base) and cnt[cl[i]] < CAPS[cl[i]]: sel.append(i); cnt[cl[i]] += 1
                proj = I[t:t + H].copy()
                for i in sel: proj += SH[i]
                pool = [i for i in ready if i not in sel]
                PA = np.array(pool); PS = SH[PA]; avail = np.ones(len(pool), bool)
                lo = In[max(0, t - W):t].min() if t > 0 else FLOOR
                tgt = ema if ema is not None else I[t]
                while len(sel) < len(base):
                    ok = avail & np.array([cnt[cl[i]] < CAPS[cl[i]] for i in pool])
                    if not ok.any(): break
                    P2 = proj[None, :] + PS
                    if pol == 'T':
                        cost = np.maximum(0.0, P2[:, :W] - (lo + TH)).sum(1)
                    elif pol in ('SOC', 'RES'):
                        # STATE-OF-CHARGE-AWARE selection. Budget = bank charge above a reserve.
                        # Among picks whose bank need fits the budget: SOC takes the HEAVIEST (spend the
                        # heavy backlog while the bank can afford it -> current climbs to the slew limit,
                        # no further: controlled ramps); RES takes the oldest. If nothing fits: min need.
                        need = np.maximum(0.0, P2[:, :W] - (lo + TH)).sum(1)
                        B = (tok + Ir - res * Qb) if Qb is not None else 0.0
                        fit = need <= B + 1e-9
                        if (fit & ok).any():
                            cost = np.where(fit, (-csum[PA] if pol == 'SOC' else 0.0), np.inf)
                        else:
                            cost = need.copy()
                    else:   # FLAT
                        cost = ((P2[:, :4] - tgt) ** 2).sum(1)
                    cost = cost + 1e-6 * np.arange(len(pool))      # tie-break: oldest
                    cost[~ok] = np.inf
                    k = int(np.argmin(cost)); i = pool[k]
                    sel.append(i); cnt[cl[i]] += 1; proj += SH[i]; avail[k] = False
                if len(sel) == len(base):
                    if set(sel) != set(base): nchg += 1
                    for i in ready:
                        if i not in sel: held[i] += 1
                    pick = sel
        for i in pick:
            ready.remove(i); tnew[i] = t; I[t:t + H] += SH[i]; done += 1; inq -= 1; icount[tid[i]] -= 1
            for c in cons[i]:
                # producer's delay vs its solo issue time propagates to the consumer (as replay2)
                pend[c] -= 1; dmax[c] = max(dmax[c], t - TR[i])
                if pend[c] == 0 and adm[c] >= 0: make_ready(c, t)
        ema = I[t] if ema is None else ema + (I[t] - ema) / ema_tau
        if Qb is not None:
            tok = min(float(Qb), tok + Ir)
            lo = In[max(0, t - W):t].min() if t > 0 else FLOOR
            ex = max(0.0, I[t] - (lo + TH)); b = min(ex, tok); tok -= b; drawn += b; unc += ex - b
            if ex - b > 1e-9: nshort += 1
            SOCt[t] = tok
            In[t] = I[t] - b
        else:
            In[t] = I[t]
        t += 1
    end = int(tnew.max()) + 20
    fin = [int(tnew[off[j]:off[j + 1]].max()) for j in range(len(threads))]
    return dict(I=I[:end], In=In[:end], makespan=int(tnew.max()), fin=fin, nready=np.array(nready),
                nchg=nchg, drawn=drawn, unc=unc, tnew=tnew, nshort=nshort, soc=SOCt[:end])
