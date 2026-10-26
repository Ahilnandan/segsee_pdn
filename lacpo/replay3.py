"""Dependency-aware issue replay for MediumBoom (2019 BOOM, chipyard 717c4658) on real ready sets from the VCD,
with per-uop LACPo-attributed current (workload.py) and predictor-driven power policies.

Uop i: ready r_i, grant t_i (real). Producers: latest instance writing (rtype, preg) granted before r_i.
Counterfactual ready r'_i = r_i + max(0, max_p (t'_p - t_p)). Baseline 'A' = oldest-first per IQ with the real port map;
pinning (r_i := t_i where hardware stalled for reasons not visible in the IQ: busy divider, LSU nacks...) is calibrated
until A reproduces every real grant. Policies use ONLY what the core knows: a data-agnostic issue-stream current estimate
(NNLS kernels) convolved with a nominal die-voltage impulse response g (mV per mA) -> forecast v_hat.
  A   : baseline
  R   : reorder (work-conserving: same issue count; choose lightest/heaviest/base set minimising forecast energy)
  H   : hold heaviest picks while forecast droop < -th (each uop held at most MAXD cycles)
  D   : dummy ALU ops on free int ports while forecast overshoot > +th_d (fills current dips / damps ringing)
  RH, RHD, HD: combinations
Background (frontend, caches, leakage) is time-stretched: BG'[n'] = BG[n' - D(n')], D = running max grant delay."""
from __future__ import print_function
import os, sys, json, numpy as np
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
K = 16; NR = 4
FU_ALU, FU_BRU, FU_MEM, FU_MUL, FU_DIV, FU_CSR, FU_FPU, FU_FDV, FU_I2F, FU_F2I = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512
PORTS = {0: [FU_ALU | FU_BRU | FU_MUL | FU_I2F, FU_ALU | FU_CSR | FU_DIV], 1: [FU_MEM], 2: [FU_FPU | FU_FDV | FU_F2I]}
ROB = 64

def nominal_g(L=0.5e-9, C=8e-9, R=50e-3, Resr=20e-3, T=3e-9, n=64, sub=12):
    """die-voltage deviation (mV, cycle-averaged) at cycle k for a 1 mA load pulse during cycle 0: fixed supply ->
    series R-L -> die node with C + ESR. States x = [iL, vC]; v_die = vC + Resr*(iL - u)."""
    from scipy.linalg import expm
    A = np.array([[-(R + Resr) / L, -1 / L], [1 / C, 0]]); Bm = np.array([Resr / L, -1 / C])
    h = T / sub; Ad = expm(A * h); Bd = np.linalg.solve(A, (Ad - np.eye(2)).dot(Bm))
    x = np.zeros(2); g = np.zeros(n)
    for k in range(n):
        u = 1e-3 if k == 0 else 0.0; acc = 0.0
        for _ in range(sub):
            x = Ad.dot(x) + Bd * u; acc += x[1] + Resr * (x[0] - u)
        g[k] = acc / sub * 1e3
    return g

def unwrap_age(rob, t):
    order = np.argsort(t, kind='stable'); age = np.zeros(len(rob), np.int64); ref = None
    for i in order:
        rb = int(rob[i])
        if ref is None: v = rb
        else:
            k = (ref - rb) // ROB; v = min((rb + (k + j) * ROB for j in (-1, 0, 1, 2)), key=lambda x: abs(x - ref))
        age[i] = v; ref = v if ref is None else max(ref, v)
    return age

def producers(z):
    n = len(z['t']); order = np.argsort(z['t'], kind='stable'); hist = defaultdict(list); prods = [[] for _ in range(n)]
    for i in order:
        r = z['r'][i]
        for p, rt in ((z['prs1'][i], z['t1'][i]), (z['prs2'][i], z['t2'][i])):
            if rt not in (0, 1) or (rt == 0 and p == 0): continue
            h = hist.get((int(rt), int(p)))
            if not h: continue
            k = len(h) - 1
            while k >= 0 and h[k][0] >= r: k -= 1
            if k >= 0: prods[i].append(h[k][1])
        if z['dtype'][i] in (0, 1) and not (z['dtype'][i] == 0 and z['pdst'][i] == 0): hist[(int(z['dtype'][i]), int(z['pdst'][i]))].append((int(z['t'][i]), i))
    cons = [[] for _ in range(n)]
    for i, ps in enumerate(prods):
        for p in set(ps): cons[p].append(i)
    return prods, cons

def port_fill(order, fu, iq, used=None):
    used = used or {q: [False] * len(p) for q, p in PORTS.items()}; out = []
    for i in order:
        q = iq[i]; P = PORTS[q]
        for w in range(len(P)):
            if not used[q][w] and (fu[i] & P[w]): used[q][w] = True; out.append(i); break
    return out, used

class Replay(object):
    def __init__(self, w):
        z = dict(np.load(B + w + '_attr.npz')); self.w = w
        self.iq = z['iq'].astype(int); self.fu = z['fu'].astype(int); self.t = z['t'].astype(np.int64); self.r0 = z['r'].astype(np.int64)
        self.cls = z['cls'].astype(int); self.A = z['A']; self.BG = z['BG'].astype(float); self.n = len(self.BG)
        self.age = unwrap_age(z['rob'], self.t); self.prods, self.cons = producers(z); self.m = len(self.t)
        h = np.load(HERE + 'hpos.npy'); self.hp = h.sum(0)              # (classes, K) tile-level predictor kernel, mA
        self.charge = self.hp.sum(1)                                     # per-class estimated charge (mA*cycles)
        alu = self.cls == 0; self.dummyA = self.A[alu].mean(0) if alu.any() else np.zeros((NR, K))
        self.g = nominal_g()
        self.pin = np.zeros(self.m, bool)

    def calibrate(self, maxit=12):
        for it in range(maxit):
            res = self.run('A', record=False)
            early = res['tnew'] < self.t
            late = res['tnew'] > self.t
            if not early.any(): break
            self.pin |= early
        res = self.run('A', record=False)
        self.calib = dict(it=it, pinned=float(self.pin.mean()), exact=float((res['tnew'] == self.t).mean()), late=float((res['tnew'] > self.t).mean()))
        return self.calib

    def run(self, pol, th=6.0, th_d=4.0, maxd=4, H=8, record=True, dummy_max=2):
        m, n = self.m, self.n; iq, fu, age, hp, g = self.iq, self.fu, self.age, self.hp, self.g
        r = np.where(self.pin, self.t, self.r0).astype(np.int64)
        pend = np.array([len(set(p)) for p in self.prods]); delay = np.zeros(m, np.int64)
        tnew = np.full(m, -1, np.int64); held = np.zeros(m, np.int64)
        wake = defaultdict(list)
        for i in range(m):
            if pend[i] == 0: wake[int(r[i])].append(i)
        Lmax = n + 200000; Ihat = np.zeros(Lmax); ready = []; done = 0; t = 0; Dmax = 0
        Kg = len(g); nd = 0; nhold = 0; nreo = 0; dummies = []
        mu = 0.0; alpha = 1.0 / 256
        def forecast(extra):
            """forecast v_hat over next H cycles (mV) from estimated current incl. the candidate picks 'extra'."""
            fut = Ihat[t:t + H + 1].copy()
            for (c, k0) in extra: fut[k0:k0 + K][:H + 1 - k0] += hp[c][:H + 1 - k0]
            lo = max(0, t - Kg); past = Ihat[lo:t] - mu; seq = np.concatenate([np.zeros(Kg - (t - lo)), past, fut - mu])
            return np.array([g[::-1].dot(seq[j + 1:j + 1 + Kg]) for j in range(H)])
        while done < m:
            if t >= Lmax - 100: raise RuntimeError('replay overflow')
            if t in wake:
                ready.extend(wake.pop(t))
            pick = []
            if ready:
                ready.sort(key=lambda i: age[i])
                base, used = port_fill(ready, fu, iq); pick = base
                if pol != 'A' and base:
                    fb = forecast([(self.cls[i], 0) for i in base])
                    if 'R' in pol:
                        light = sorted(ready, key=lambda i: (self.charge[self.cls[i]], age[i]))
                        heavy = sorted(ready, key=lambda i: (-self.charge[self.cls[i]], age[i]))
                        best = (float((fb ** 2).sum()), base)
                        for ordr in (light, heavy):
                            forced = [i for i in ready if held[i] >= maxd]
                            sel, _ = port_fill(forced + [i for i in ordr if i not in forced], fu, iq)
                            if len(sel) != len(base): continue
                            f = forecast([(self.cls[i], 0) for i in sel]); s = float((f ** 2).sum())
                            if s < best[0] - 1e-9: best = (s, sel)
                        if set(best[1]) != set(base): nreo += 1
                        pick = best[1]; fb = forecast([(self.cls[i], 0) for i in pick])
                    if 'H' in pol and fb.min() < -th:
                        keep = [i for i in pick if held[i] >= maxd]
                        rest = sorted([i for i in pick if held[i] < maxd], key=lambda i: self.charge[self.cls[i]])
                        for i in rest:
                            f = forecast([(self.cls[j], 0) for j in keep + [i]])
                            if f.min() >= -th: keep.append(i)
                        if len(keep) < len(pick): nhold += 1
                        pick = keep
            if 'D' in pol:
                fd = forecast([(self.cls[i], 0) for i in pick])
                if fd.max() > th_d:
                    used = {q: [False] * len(p) for q, p in PORTS.items()}
                    port_fill(pick, fu, iq, used)
                    free = sum(1 for x in used[0] if not x)
                    for _ in range(min(free, dummy_max)):
                        dummies.append(t); Ihat[t:t + K] += hp[0]; nd += 1
            for i in pick:
                ready.remove(i); tnew[i] = t; done += 1
                Ihat[t:t + K] += hp[self.cls[i]]
                d = t - self.t[i]; Dmax = max(Dmax, d)
                for c in self.cons[i]:
                    delay[c] = max(delay[c], d); pend[c] -= 1
                    if pend[c] == 0: wake[int(max(r[c] + max(delay[c], 0), t + 1))].append(c)
            for i in ready:
                if i not in pick: held[i] += 1 if pol != 'A' else 0
            mu = (1 - alpha) * mu + alpha * Ihat[t]
            t += 1
        res = dict(tnew=tnew, nhold=nhold, nreo=nreo, ndummy=nd)
        if not record: return res
        # ---- build the region current on the new timeline
        Ln = int(max(tnew.max() + K + 1, n + max(0, int(np.percentile(tnew - self.t, 99.5)) + 64)))
        # background time-warp: pipeline slip = rolling median grant delay (last 64 grants); the background clock
        # advances at most one cycle per cycle and never goes backwards (repeats while the core is held)
        order = np.argsort(tnew, kind='stable'); dl = (tnew - self.t)[order]; tg = tnew[order]
        import pandas as _pd
        med = _pd.Series(dl).rolling(64, min_periods=1).median().values
        slip = np.zeros(Ln); slip[tg] = med; has = np.zeros(Ln, bool); has[tg] = True
        idx = np.where(has, np.arange(Ln), 0); idx = np.maximum.accumulate(idx); slip = np.maximum(slip[idx], 0)
        src = np.zeros(Ln, np.int64); m_ = 0
        for k in range(Ln):
            tgt = k - int(slip[k]); m_ = max(0, tgt) if k == 0 else max(m_, min(m_ + 1, tgt)); src[k] = min(m_, n - 1)
        D = np.arange(Ln) - src
        I = self.BG[src].copy()
        for k in range(K):
            np.add.at(I, tnew + k, self.A[:, :, k].astype(float))
        for td in dummies:
            hi = min(Ln, td + K); I[td:hi] += self.dummyA.T[:hi - td]
        Ih = np.zeros((Ln, NR))                                          # predictor's view per region (for device triggers)
        hreg = np.load(HERE + 'hpos.npy')
        for k in range(K):
            np.add.at(Ih, tnew + k, hreg[:, self.cls, k].T)
        t0 = int(np.percentile(self.t, 2)); span = max(1, int(self.t.max()) - t0)
        res.update(I=I.astype(np.float32), Ihat=Ih.astype(np.float32), slowdown=float(D[-1]) / span,
                   cycles=Ln, dummy_charge=float(self.dummyA.sum() * nd), exec_delay=int(D[-1]))
        return res

VARIANTS = {'A': ('A', {}), 'R': ('R', {}),
            'H4': ('H', dict(th=4.0)), 'H6': ('H', dict(th=6.0)), 'H10': ('H', dict(th=10.0)),
            'D3': ('D', dict(th_d=3.0)), 'D5': ('D', dict(th_d=5.0)), 'D8': ('D', dict(th_d=8.0)),
            'RH6': ('RH', dict(th=6.0)), 'HD6': ('HD', dict(th=6.0, th_d=5.0)), 'RHD6': ('RHD', dict(th=6.0, th_d=5.0)),
            'RHD4': ('RHD', dict(th=4.0, th_d=3.0))}

if __name__ == '__main__':
    w = sys.argv[1]; vs = sys.argv[2].split(',') if len(sys.argv) > 2 else list(VARIANTS)
    rp = Replay(w); cal = rp.calibrate(); print(w, 'calibration', cal, flush=True)
    out = {}
    for v in vs:
        pol, kw = VARIANTS[v]; res = rp.run(pol, **kw)
        np.savez_compressed(B + '%s_rep_%s.npz' % (w, v), I=res['I'], Ihat=res['Ihat'])
        out[v] = dict(slowdown=res['slowdown'], nhold=res['nhold'], nreo=res['nreo'], ndummy=res['ndummy'], cycles=res['cycles'],
                      dummy_charge=res['dummy_charge'], exec_delay=res['exec_delay'])
        print(w, v, out[v], flush=True)
    json.dump(dict(calib=cal, runs=out), open(B + '%s_rep.json' % w, 'w'), indent=1)
