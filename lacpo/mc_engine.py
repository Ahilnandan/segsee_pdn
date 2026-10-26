"""mc_engine: closed-loop MULTI-CORE cross-core issue coordinator on the shared pdn7 die (2x2 cores, 256 tiles).

Core   : per-core issue stepper that reproduces policy6.Engine.run(design, 'A', prm) issue behaviour exactly (same ready
         list, wake / dependency handling, replay3.port_fill, ROB window + in-order commit, BG time-warp via slip/src,
         I_t = BG[src] + ring[t].dot(WB)). Split into cand(t) -> greedy age-order port pick, and issue(t, pick) -> I_t (64).
         Ihat holds each core's predicted issue-coupled current (hp kernel per uop class) for the committed uops.
run_mc : 4 Cores + one pdn7.Sim stepped together each cycle (per-core I_t placed into the die).
Coordinator (coord dict, None = no coordination):
  every cycle, S[t+j] = sum over cores of committed predicted issue current (j = 0..2); the candidate picks of all cores
  are admitted greedily (forced uops first; then cores ordered most-lagged-first with a rotating tie-break, within a core
  lightest-first as in policy6's hold loop) while  max_j S[t+j] - min(S[t-w .. t-1]) <= budget_eff  (rise of the
  predicted TOTAL current over 1-3 cycles). Picks that do not fit are held one cycle.
  Caps: per-uop hold cap maxd (held >= maxd -> forced), per-core throughput cap (EWMA, 64 cycles, of cycles in which the
  core had a pick held; above cap -> all that core's picks forced this cycle).
  Optional sensor: budget_eff = budget - ks * max(0, sensed droop - v0) using the worst-tile voltage delayed 2 cycles and
  quantised to 5 mV.
"""
from __future__ import print_function
import os, sys, collections, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay3, timing6 as TM, pdn7
K = replay3.K


class Core(object):
    def __init__(self, eng, a, b, rob=80):
        rp = eng.rp; self.eng = eng; self.a = a; self.b = b
        sel = np.where((rp.t >= a) & (rp.t < b))[0]; loc = {g_: k for k, g_ in enumerate(sel)}; m = len(sel)
        self.sel, self.m = sel, m
        self.iq, self.fu = rp.iq[sel], rp.fu[sel]; self.cls = rp.cls[sel]; self.age = rp.age[sel]; self.tt = rp.t[sel] - a
        r = (np.where(rp.pin[sel], rp.t[sel], rp.r0[sel]) - a).astype(np.int64); self.r = np.maximum(r, 0)
        self.prods = [[loc[p] for p in set(rp.prods[g_]) if p in loc] for g_ in sel]; self.cons = [[] for _ in range(m)]
        for i, ps in enumerate(self.prods):
            for p_ in ps: self.cons[p_].append(i)
        self.pend = np.array([len(ps) for ps in self.prods]); self.delay = np.zeros(m, np.int64)
        self.tnew = np.full(m, -1, np.int64); self.held = np.zeros(m, np.int64)
        self.wake = collections.defaultdict(list)
        for i in range(m):
            if self.pend[i] == 0: self.wake[int(self.r[i])].append(i)
        L = 3 * (b - a) + 20000; self.L = L
        self.ring = np.zeros((L + K, eng.A5.shape[1])); self.Ihat = np.zeros(L + K); self.BG = eng.BG[a:b]
        uname = [TM.FU.get(int(c), 'ALU') for c in self.fu]; self.lat = np.array([TM.LAT[u] for u in uname])
        self.ROBN = int(rob); self.order = np.argsort(self.age, kind='stable'); self.rank = np.empty(m, np.int64)
        self.rank[self.order] = np.arange(m)
        self.commit = 0; self.rob_wait = []; self.ready = []; self.done = 0; self.slip = 0.0; self.src = 0
        self.kw = eng.hp.sum(1)[self.cls]                      # predicted charge per uop (lightness key)
        self.tend = None

    def active(self, t):
        return self.done < self.m or t < (self.b - self.a)

    def cand(self, t):
        """ROB commit / admission, wakeups, then the greedy age-order port pick (policy6 'A')."""
        if t >= self.L - 2: raise RuntimeError('overflow')
        tnew, rank, ready = self.tnew, self.rank, self.ready
        if self.ROBN:
            while self.commit < self.m:
                i_ = self.order[self.commit]
                if tnew[i_] >= 0 and tnew[i_] + 2 + self.lat[i_] <= t: self.commit += 1
                else: break
            if self.rob_wait:
                lim = self.commit + self.ROBN
                ok_ = [i for i in self.rob_wait if rank[i] < lim]
                if ok_:
                    self.rob_wait = [i for i in self.rob_wait if rank[i] >= lim]
                    ready.extend(i for i in sorted(ok_, key=lambda x: rank[x]) if self.pend[i] == 0 and tnew[i] < 0 and i not in ready)
        if t in self.wake:
            for i in self.wake.pop(t):
                if self.pend[i] == 0 and tnew[i] < 0 and i not in ready:
                    if self.ROBN and rank[i] >= self.commit + self.ROBN:
                        if i not in self.rob_wait: self.rob_wait.append(i)
                    else: ready.append(i)
        if not ready: return []
        ready.sort(key=lambda i: self.age[i]); base, _ = replay3.port_fill(ready, self.fu, self.iq)
        return base

    def issue(self, t, pick, held=()):
        eng = self.eng
        for i in pick:
            self.ready.remove(i); self.tnew[i] = t; self.done += 1
            self.ring[t:t + K] += eng.A5[self.sel[i]].T; self.Ihat[t:t + K] += eng.hp[self.cls[i]]
            d = t - self.tt[i]; self.slip += (max(d, 0) - self.slip) / 32.0
            for c in self.cons[i]:
                self.delay[c] = max(self.delay[c], d); self.pend[c] -= 1
                if self.pend[c] == 0: self.wake[int(max(self.r[c] + max(self.delay[c], 0), t + 1))].append(c)
        for i in held: self.held[i] += 1
        tgt = t - int(round(self.slip)); src = max(0, tgt) if t == 0 else max(self.src, min(self.src + 1, tgt))
        self.src = min(src, len(self.BG) - 1)
        if self.tend is None and not (self.done < self.m or t + 1 < (self.b - self.a)): self.tend = t + 1
        return self.BG[self.src] + self.ring[t].dot(eng.WB)

    def slowdown(self):
        return (int(self.tnew.max()) - int(self.tt.max())) / max(1, int(self.tt.max()))


CDEF = dict(budget=None, w=3, maxd=4, cap=0.25, ks=0.0, v0=30.0)


def coordinate(cores, cands, t, S_hist, st, prm, vsens):
    """returns per-core (pick, held). st: per-core state dict (hold ewma, lag)."""
    nc = len(cores); fut = sum(c.Ihat[t:t + 3] for c in cores)
    lo = min(S_hist[max(0, t - prm['w']):t]) if t > 0 else 0.0
    bud = prm['budget'] - prm['ks'] * max(0.0, (pdn7.VNOM - vsens) * 1e3 - prm['v0'])
    keep = [[] for _ in range(nc)]; rest = []
    for k, c in enumerate(cores):
        over = st['ew'][k] > prm['cap']
        for i in cands[k]:
            if over or c.held[i] >= prm['maxd']:
                keep[k].append(i); fut = fut + c.eng.hp[c.cls[i]][:3]
            else: rest.append((k, i))
    if rest:
        rot = t % nc
        corder = sorted(range(nc), key=lambda k: (-st['lag'][k], (k - rot) % nc))
        pri = {k: n for n, k in enumerate(corder)}
        rest.sort(key=lambda x: (pri[x[0]], cores[x[0]].kw[x[1]], cores[x[0]].age[x[1]]))
        for (k, i) in rest:
            f2 = fut + cores[k].eng.hp[cores[k].cls[i]][:3]
            if f2.max() - lo <= bud: keep[k].append(i); fut = f2
    out = []
    for k in range(nc):
        ks_ = set(keep[k]); hd = [i for i in cands[k] if i not in ks_]
        h = 1.0 if hd else 0.0; st['ew'][k] += (h - st['ew'][k]) / 64.0; st['lag'][k] += len(hd)
        out.append(([i for i in cands[k] if i in ks_], hd))
    return out


def run_mc(engs, starts, design1, W=60000, coord=None, rob=80, warm=600, q_mV=5.0, sensor_delay=2):
    nc = len(engs); assert nc == pdn7.NCORE
    cores = [Core(engs[k], starts[k], starts[k] + W, rob) for k in range(nc)]
    d4 = pdn7.place_design(design1); has_bank = bool(d4.get('bank'))
    sim = pdn7.Sim(d4, Im_A=pdn7.place([c.BG[:3000].mean(0) for c in cores]) * 1e-3)
    prm = dict(CDEF, **(coord or {})); on = coord is not None and prm['budget'] is not None
    st = dict(ew=[0.0] * nc, lag=[0] * nc); S_hist = []; vq = collections.deque([pdn7.VNOM] * sensor_delay, maxlen=sensor_delay)
    vw, itot, ipred = [], [], []; arm = np.zeros((1, pdn7.T), np.uint8); vref = np.array([pdn7.VNOM]); nheld = 0; t = 0
    while any(c.active(t) for c in cores):
        cands = [c.cand(t) for c in cores]
        if on and any(cands): pk = coordinate(cores, cands, t, S_hist, st, prm, vq[0])
        else: pk = [(cd, []) for cd in cands]
        It = [c.issue(t, p, h) for c, (p, h) in zip(cores, pk)]; nheld += sum(len(h) for _, h in pk)
        S_hist.append(sum(c.Ihat[t] for c in cores))
        Id = pdn7.place(It)
        vmin, _, _ = sim.run(Id[None, :].astype(np.float32), vref=vref, arm=arm, want_res=has_bank)
        v = float(vmin[0].min()); vw.append(v); itot.append(float(Id.sum()))
        vq.append(np.round(v / (q_mV * 1e-3)) * q_mV * 1e-3); t += 1
    tend = min(c.tend for c in cores)
    dr = (pdn7.VNOM - np.array(vw[warm:tend])) * 1e3; it = np.array(itot[warm:tend]); S = np.array(S_hist)
    rise = np.array([S[j:j + 3].max() - S[max(0, j - 3):j].min() for j in range(max(warm, 3), tend - 3)])
    return dict(p999=float(np.percentile(dr, 99.9)), max=float(dr.max()), p99=float(np.percentile(dr, 99)),
                dI999=float(np.percentile(np.abs(np.diff(it)), 99.9)), rise=[float(np.percentile(rise, q)) for q in (50, 90, 99, 99.9)],
                sd=[c.slowdown() for c in cores], exec=[int(c.tnew.max()) for c in cores], nheld=nheld, cycles=t, tend=tend)


def run_single(eng, a, b, design, rob=80, warm=600):
    """gate helper: one Core through pdn6.Sim, identical bookkeeping to policy6 'A'."""
    import pdn6
    c = Core(eng, a, b, rob); sim = pdn6.Sim(design, Im_A=c.BG[:min(3000, len(c.BG))].mean(0) * 1e-3)
    has_bank = bool(design.get('bank')); arm = np.zeros((1, pdn6.T), np.uint8); V = []; t = 0
    while c.active(t):
        base = c.cand(t); I_t = c.issue(t, base)
        vmin, _, _ = sim.run(I_t[None, :].astype(np.float32), vref=np.array([pdn6.VNOM]), arm=arm, want_res=has_bank, want_state=True)
        V.append(vmin[0]); t += 1
    worst = ((pdn6.VNOM - np.array(V[warm:])) * 1e3).max(1)
    return dict(p999=float(np.percentile(worst, 99.9)), max=float(worst.max()), slowdown=c.slowdown(), tnew=c.tnew)
