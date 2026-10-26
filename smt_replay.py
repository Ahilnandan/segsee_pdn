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

def pdn_predictor(L=507e-12, C=5e-9, R=79.6e-3, Hp=12, sub=4, dt=1e-9):
    """Exact-ZOH series-RLC (same model as replay2.droop). Returns per-cycle transition (Ad1, Bd1) and
    horizon maps: v[substeps] = Fv @ x0 + G @ u[0:Hp]  (v = die-voltage deviation, droop = -v)."""
    from scipy.linalg import expm
    A = np.array([[-R / L, -1 / L], [1 / C, 0.0]]); Bv = np.array([0.0, -1 / C]); h = dt / sub
    M = np.zeros((3, 3)); M[:2, :2] = A * h; M[:2, 2] = Bv * h; E = expm(M); Ad = E[:2, :2]; Bd = E[:2, 2]
    def sim(x0, u):
        x = np.array(x0, float); out = []
        for k in range(Hp):
            for _ in range(sub): x = Ad @ x + Bd * u[k]; out.append(x[1])
        return np.array(out)
    Fv = np.stack([sim([1, 0], np.zeros(Hp)), sim([0, 1], np.zeros(Hp))], 1)
    G = np.stack([sim([0, 0], np.eye(Hp)[j]) for j in range(Hp)], 1)
    Ad1 = np.linalg.matrix_power(Ad, sub); Bd1 = sum(np.linalg.matrix_power(Ad, j) @ Bd for j in range(sub))
    return dict(Fv=Fv, G=G, Ad1=Ad1, Bd1=Bd1, Hp=Hp, R=R)

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

def run(threads, pol='A', LEAD=2, MAXD=8, TH=100.0, ema_tau=64.0, iq=IQ, Qb=None, Ir=25.0, record_ready=True, res=0.0, pdn=None, vth=None, hold=None, hth=1.0, ramp=None, est=None, gro=None):
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
    # --- what the HARDWARE predictor knows (est): 4-bit-amp x 4-bit-dur power tokens (rectangle of the
    # token's quantised mean amplitude, no intra-op timing), a possibly mismatched PDN model, and a
    # state sensed `lat` cycles late and propagated forward with the token-estimated current.
    est = est or {}
    if est.get('quant'):
        amax = 80.0; step = amax / 15
        SHe = np.zeros_like(SH)
        for i in range(N):
            nz = np.nonzero(SH[i])[0]
            if len(nz) == 0: continue
            dur = min(H, int(nz[-1]) + 1); amp = min(15, round(SH[i].sum() / dur / step)) * step
            SHe[i, :dur] = amp
    else:
        SHe = SH
    Ie = np.full(L, FLOOR); Ine = np.full(L, FLOOR)
    FIR = pdn is not None and pdn.get('kind') == 'fir'
    if FIR:
        # identified per-cycle impulse response of the REAL supply (Pulsim): v_dev[t] = sum_j h[j] * (I[t-j] - dc)
        hk = np.asarray(pdn['h'], float); K = len(hk); Hp = pdn['Hp']
        Hm = np.zeros((Hp, H))                         # effect of a shape issued at t on cycles t..t+Hp-1
        for c in range(Hp):
            for m in range(min(H, c + 1)): Hm[c, m] = hk[c - m]
        dcI = FLOOR; dc_tau = 2000.0
        def fir_base(tt, upto=None):
            # predicted deviation at cycles tt..tt+Hp-1 from history (< tt) + committed tails (>= tt)
            lo_ = max(0, tt - K + 1)
            x = np.concatenate([Ine[lo_:tt], Ie[tt:tt + Hp]]) - dcI
            if upto is not None: x[(upto - lo_):] = 0.0             # Grochowski latency: nothing known after upto
            x = np.concatenate([np.zeros(K - 1 - (tt - lo_)), x])
            return np.array([hk[::-1] @ x[c:c + K] for c in range(Hp)])
        pdn = dict(pdn, R=0.0)
    if pdn is not None:
        Hp = pdn['Hp']; i00 = FLOOR * 1e-3; xs = np.array([i00, -pdn['R'] * i00]) if not FIR else np.zeros(2)
        pe = est.get('pdn', pdn); lat = int(est.get('lat', 0)); XH = []
        xh = xs.copy()
    gro = gro or None; dumq = 0.0
    def predJ(sel):
        if FIR:
            sh_ = np.zeros(H)
            for i in sel: sh_ += SHe[i]
            return -float((fir_base(t) + Hm @ sh_).min())
        u = Ie[t:t + Hp].copy()
        for i in sel: u[:H] += SHe[i, :min(H, Hp)]
        bud = (tok - res * Qb) if Qb is not None else 0.0
        for c in range(Hp):
            if c > 0 and Qb is not None: bud += Ir
            ex = max(0.0, u[c] - (lo + TH)); b = min(ex, max(bud, 0.0)); u[c] -= b; bud -= b
        return -float((pe['Fv'] @ xh + pe['G'] @ (u * 1e-3)).min()) * 1e3
    nheld = 0; prevcnt = 0
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
        if pdn is not None and not FIR:
            XH.append(xs.copy())
            if lat == 0 and pe is pdn and SHe is SH: xh = xs
            else:
                xh = XH[max(0, t - lat)].copy()
                for j in range(max(0, t - lat), t): xh = pe['Ad1'] @ xh + pe['Bd1'] * (Ine[j] * 1e-3)
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
                lo = Ine[max(0, t - W):t].min() if t > 0 else FLOOR
                tgt = ema if ema is not None else I[t]
                while len(sel) < len(base):
                    ok = avail & np.array([cnt[cl[i]] < CAPS[cl[i]] for i in pool])
                    if not ok.any(): break
                    P2 = proj[None, :] + PS
                    if pol == 'T':
                        cost = np.maximum(0.0, P2[:, :W] - (lo + TH)).sum(1)
                    elif pol in ('MPC', 'MPCmin'):
                        # MODEL-PREDICTIVE selection: predicted die-voltage over Hp cycles from the PDN state,
                        # committed current (in-flight op tails), the candidate, and the bank covering the slew
                        # excess with its charge minus a reserve.
                        Pe = Ie[t:t + H].copy()
                        for i in sel: Pe += SHe[i]
                        U = np.repeat(Ie[None, t:t + Hp].copy(), len(pool), 0); U[:, :H] = Pe[None, :] + SHe[PA]
                        bud = np.full(len(pool), (tok - res * Qb) if Qb is not None else 0.0)
                        for c in range(Hp):
                            bud += Ir if (Qb is not None and c > 0) else 0.0
                            ex = np.maximum(0.0, U[:, c] - (lo + TH)); b = np.minimum(ex, np.maximum(bud, 0.0))
                            U[:, c] -= b; bud -= b
                        V = (pe['Fv'] @ xh)[None, :] + (U * 1e-3) @ pe['G'].T
                        J = -V.min(1) * 1e3                       # predicted worst droop (mV) over the horizon
                        if pol == 'MPC' and vth is not None and ((J <= vth) & ok).any():
                            cost = np.where(J <= vth, -csum[PA], np.inf)   # safe: spend the heavy backlog
                        else:
                            cost = J
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
        if gro is not None:
            # GROCHOWSKI et al. HPCA'02 di/dt controller: voltage computed by convolving (estimated) current with
            # the PDN impulse response (== the state model), `lat` cycles late; below the lower threshold stop
            # issue (all non-forced picks); above the upper threshold force units on (dummy ALU ops on free slots).
            vd = (-float(fir_base(t, upto=t - int(gro.get('lat', 0)))[0]) if FIR else -float(xh[1]) * 1e3)
            if vd > gro['lo']:
                kept = [i for i in pick if held[i] >= MAXD]
                if len(kept) < len(pick):
                    nheld += 1
                    for i in pick:
                        if i not in kept: held[i] += 1
                pick = kept
            elif gro.get('hi') is not None and vd < gro['hi']:
                nd = ISSUE_W - len(pick)
                if nd > 0:
                    DS = np.array([10., 20., 20., 10.]); I[t:t + 4] += nd * DS; Ie[t:t + 4] += nd * DS; dumq += nd * DS.sum()
        if hold and pick:
            # NON-work-conserving: hold the heaviest non-forced picks while the criterion says so
            lo = Ine[max(0, t - W):t].min() if t > 0 else FLOOR
            forced = [i for i in pick if held[i] >= MAXD]
            rest = sorted([i for i in pick if held[i] < MAXD], key=lambda i: csum[i])   # light first
            if hold == 'MPCH':
                kept = forced + rest
                while len(kept) > len(forced) and predJ(kept) > vth: kept.pop()          # drop heaviest
            elif hold == 'TDH':
                def needf(sel):
                    proj = I[t:t + H].copy()
                    for i in sel: proj += SH[i]
                    return float(np.maximum(0.0, proj[:W] - (lo + TH * hth)).sum())
                kept = forced + rest
                while len(kept) > len(forced) and needf(kept) > tok + Ir: kept.pop()
            elif hold == 'RAMP':
                kept = forced + sorted([i for i in pick if held[i] < MAXD], key=lambda i: (adm[i], i))
                kept = kept[:max(len(forced), prevcnt + ramp)]
            if len(kept) < len(pick):
                nheld += 1
                for i in pick:
                    if i not in kept: held[i] += 1
            pick = kept
        prevcnt = len(pick)
        for i in pick:
            ready.remove(i); tnew[i] = t; I[t:t + H] += SH[i]; Ie[t:t + H] += SHe[i]; done += 1; inq -= 1; icount[tid[i]] -= 1
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
            In[t] = I[t] - b; Ine[t] = Ie[t] - b
        else:
            In[t] = I[t]; Ine[t] = Ie[t]
        if pdn is not None and not FIR: xs = pdn['Ad1'] @ xs + pdn['Bd1'] * (In[t] * 1e-3)
        if FIR: dcI += (Ine[t] - dcI) / dc_tau
        t += 1
    end = int(tnew.max()) + 20
    fin = [int(tnew[off[j]:off[j + 1]].max()) for j in range(len(threads))]
    return dict(I=I[:end], In=In[:end], makespan=int(tnew.max()), fin=fin, nready=np.array(nready),
                nchg=nchg, drawn=drawn, unc=unc, tnew=tnew, nshort=nshort, soc=SOCt[:end], nheld=nheld, dumq=dumq, Itot=float(I[:end].sum()))
