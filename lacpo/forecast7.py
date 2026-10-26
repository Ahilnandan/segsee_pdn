"""forecast7: learned short-horizon (H=8) supply-voltage forecaster for the policy6 closed loop.

Model: per design, per horizon j, ridge regression (lam 0.03 on standardized features, chosen by leave-one-train-bench-out
CV) on causal features (feature set BEST = GROUPS['F E+cls'], ~600 inputs): measured current lags, LDO-bank current lags,
committed predicted issue current, learned per-class issue history (n<c>-<d>), predicted next two greedy selects pushed
through the design kernel (iqp), 2-cycle-old quantized sensor, cache miss/refill flags (24 lags) and all cache-pipeline
activity flags (8 lags), and the existing design-kernel forecast (phys) as inputs.  lag1 weights use cache flags of
cycles <= t-1 only; lag0 weights also use the flags of cycle t (registered pipeline-valid bits of the current cycle).
Pieces
  record(bench, design, a, b)  runs policy6.Engine (pol 'A': schedule independent of the forecast) and logs per-cycle
                               signals through non-invasive hooks (recording Sim subclass, BG-cursor ndarray subclass,
                               port_fill wrapper for the ready list). policy6.py is not modified.
  feats / target / Acc         causal feature matrix, target v[t+j] - EWMA256(v)[t] (mV), streaming ridge.
  train(design, lag0)          fit + save forecast7_w_<design>_<lag>_<train7|all10>.npz next to this file.
  Forecaster7(wfile)           run-time object: begin() / forecast() / end() per cycle, E[c] per-class candidate effect.
  verify()                     steps Forecaster7 through a recorded run; agrees with the offline path to 1e-5 mV.
Usage: python forecast7.py rec | eval [lam] [letters|-] [lag0|lag1] | train | verify
Held-out protocol: TRAIN benches (7, windows below) for fitting, TEST = qsort/mm/towers windows for scoring only."""
from __future__ import print_function
import os, sys, re, json, collections, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
sys.path.insert(0, HERE)
SCR = '/tmp/claude-1004/scratch/fc7/'
CB = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
H = 8
TEST = dict(qsort=(100000, 160000), mm=(50000, 110000), towers=(0, None))
TRAIN = dict(rsort=(100000, 160000), dhrystone=(100000, 160000), spmv=(40000, 100000), multiply=(10000, 70000),
             memcpy=(10000, 70000), median=(0, None), vvadd=(0, None))
BASE = dict(fc_mode='meas', fc_kernel='design', fc_obs=0.5, fc_log=1)

# ------------------------------------------------------------------ recording -------------------------------------------------
class _Cur(np.ndarray):
    """BG trace that logs every scalar index (the engine's src cursor)."""
    def __array_finalize__(self, obj):
        self._log = getattr(obj, '_log', None)
    def __getitem__(self, k):
        if isinstance(k, (int, np.integer)): self._log.append(int(k))
        return np.ndarray.__getitem__(self, k)

def record(bench, dn, a=0, b=None, eng=None):
    import policy6 as P6, eval6
    D, prmb = eval6.designs(); e = eng or P6.Engine(bench)
    Sim0 = P6.pdn5.Sim; log = dict(I=[], Ib=[], v=[])
    class RecSim(Sim0):
        def run(self, I_mA, *aa, **kw):
            out = Sim0.run(self, I_mA, *aa, **kw)
            log['I'].append(float(I_mA.sum())); log['Ib'].append(1e3 * float(self.state[0][2])); log['v'].append(float(out[0][0].min()))
            return out
    BG0 = e.BG; cur = BG0.view(_Cur); cur._log = []
    PF0 = P6.replay3.port_fill; rdy = []
    def PFrec(order, *aa, **kw):                         # issue-queue ready list seen by the select logic at cycle t
        rdy.append((len(log['I']), list(order))); return PF0(order, *aa, **kw)
    P6.pdn5.Sim = RecSim; e.BG = cur; P6.replay3.port_fill = PFrec
    try:
        r = e.run(D[dn], 'A', dict(prmb, **BASE), a=a, b=b, record=True)
    finally:
        P6.pdn5.Sim = Sim0; e.BG = BG0; P6.replay3.port_fill = PF0
    b = e.n if b is None else b
    rp = e.rp; sel = np.where((rp.t >= a) & (rp.t < b))[0]; tnew = r['tnew']; n = len(log['I'])
    src = np.array(cur._log[-n:], np.int64)
    assert len(src) == n and len(r['vavg_t']) == n
    K = e.hp.shape[1]
    # committed predicted issue current, by issue cycle (per uop class kernel hp) and true issue current
    cls = rp.cls[sel]; inc = np.zeros((n + K, K)); ncls = np.zeros((n + K, e.hp.shape[0]))
    np.add.at(inc, tnew, e.hp[cls]); np.add.at(ncls, (tnew, cls), 1)
    Iiss = np.zeros(n + K); A = np.einsum('ubk,b->uk', e.A5[sel], e.WB.sum(1))      # (uops, K) true mA per cycle offset
    for k in range(K): np.add.at(Iiss, tnew + k, A[:, k])
    # issue-queue lookahead (all known to the scheduler at cycle t):
    #   BL  = ready uops left in the queue after this cycle's select (assumed to issue at t+1)
    #   WK1/WK2 = uops whose wakeup is already scheduled for t+1 / t+2 by producers that issued at <= t
    m = len(sel); loc = {g_: k for k, g_ in enumerate(sel)}
    lastp = np.full(m, -1, np.int64)
    for i, g_ in enumerate(sel):
        ps = [loc[p] for p in set(rp.prods[g_]) if p in loc]
        if ps: lastp[i] = max(tnew[p] for p in ps)
    first = np.full(m, 1 << 40, np.int64); rl = {}
    for (tc, lst) in rdy:
        rl[tc] = lst
        for i in lst:
            if first[i] > tc: first[i] = tc
    w1 = collections.defaultdict(list); w2 = collections.defaultdict(list)
    for i in range(m):
        if lastp[i] < 0 or first[i] >= n + 2: continue
        if first[i] - 1 >= lastp[i]: w1[int(first[i]) - 1].append(i)
        if first[i] - 2 >= lastp[i]: w2[int(first[i]) - 2].append(i)
    # predicted next selects (greedy age-ordered port fill, as the select logic would do): S1 at t+1, S2 at t+2
    fu, iq, age = rp.fu[sel], rp.iq[sel], rp.age[sel]; PF = P6.replay3.port_fill
    BL = np.zeros((n, H + 1)); S1 = np.zeros((n, H + 1)); S2 = np.zeros((n, H + 1)); hit = [0, 0]
    for t in range(n):
        bl = [i for i in rl.get(t, ()) if tnew[i] != t]
        for i in bl: BL[t, 1:] += e.hp[cls[i]][:H]
        c1 = sorted(set(bl) | set(w1.get(t, ())), key=lambda i: age[i]); n1, _ = PF(c1, fu, iq)
        for i in n1: S1[t, 1:] += e.hp[cls[i]][:H]; hit[0] += tnew[i] == t + 1
        s1 = set(n1); c2 = sorted(set(i for i in c1 if i not in s1) | set(w2.get(t, ())), key=lambda i: age[i]); n2, _ = PF(c2, fu, iq)
        for i in n2: S2[t, 2:] += e.hp[cls[i]][:H - 1]
        hit[1] += len(n1)
    print('next-select precision %.3f, recall %.3f' % (hit[0] / max(1, hit[1]), hit[0] / max(1, int(((tnew >= 1) & (tnew < n)).sum()))), flush=True)
    z = np.load(CB + bench + '_cache.npz'); ev = {k: (z[k][a:b] > 0).astype(np.int8) for k in ('dc_mshr', 'ic_miss', 'dc_refill', 'ic_refill', 'dc_req', 'dc_wb')}
    return dict(bench=bench, design=dn, I=np.array(log['I']), Ib=np.array(log['Ib']), v=np.array(log['v']), src=src,
                inc=inc[:n], ncls=ncls[:n], BL=BL, S1=S1, S2=S2, a=a, Iiss=Iiss[:n], flog=r['flog'], vw=r['vavg_t'], g=np.asarray(e.kernel(D[dn])),
                **{'ev_' + k: x for k, x in ev.items()})

def _rec_job(args):
    bench, win = args
    import policy6 as P6
    e = P6.Engine(bench); out = []
    for dn in ('decap', 'banks'):
        f = SCR + '%s_%s.npz' % (bench, dn)
        if os.path.exists(f): continue
        rec = record(bench, dn, win[0], win[1], eng=e); np.savez(f, **rec); out.append(f)
        print('rec', bench, dn, len(rec['I']), flush=True)
    return out

def load(bench, dn):
    z = np.load(SCR + '%s_%s.npz' % (bench, dn)); return {k: z[k] for k in z.files}

# ------------------------------------------------------------------ features --------------------------------------------------
def ewma(x, n=256.0):
    y = np.empty(len(x)); m = x[0]
    for i in range(len(x)): m += (x[i] - m) / n; y[i] = m
    return y

def lag(x, k):
    """x[t - k] (zero before start); negative k = future (targets only)."""
    y = np.zeros(len(x))
    if k == 0: return np.array(x, float)
    if k > 0: y[k:] = x[:len(x) - k]
    else: y[:k] = x[-k:]
    return y

NP = 24; NB = 8; NS = 6; NE = 24; NX = 8; NC = 10; KC = 16
MISS = ('dc_mshr', 'ic_miss', 'dc_refill', 'ic_refill', 'dc_req', 'dc_wb')
_CZ = {}
def cache_sigs(bench):
    if bench not in _CZ:
        z = np.load(CB + bench + '_cache.npz'); _CZ.clear()
        _CZ[bench] = [(k, z[k].astype(np.float32) if k.endswith('ways') else (z[k] > 0).astype(np.float32)) for k in sorted(z.files)]
    return _CZ[bench]

def snapshot(inc):
    """Ihat[t + k] as known at cycle t (uops issued at <= t), k = 0..H."""
    n, K = inc.shape; snap = np.zeros((n, H + 1))
    for k in range(H + 1):
        for d in range(K - k): snap[:, k] += lag(inc[:, d + k], d)
    return snap

def streams(rec):
    """Per-cycle raw streams exactly as the run-time forecaster sees them (index = cycle t; all causal)."""
    Im = rec['I'] - rec['Ib']; snap = snapshot(rec['inc']); q = 5.0
    return dict(Im=Im, Ib=rec['Ib'], snap=snap, vs=np.round(rec['v'] * 1e3 / q) * q, ncls=rec['ncls'],
                iq=np.hstack([rec['S1'], rec['S2'], rec['BL'] - rec['S1']]), src=rec['src'])

def feats(rec, groups, ev_lag0=False):
    """Causal features at cycle t (after the select of cycle t, before its PDN step).
    cur : measured current Im = I - Ibank of cycles t-1..t-NP, minus EWMA256 (as of t-1)
    bank: measured LDO-bank current of cycles t-1..t-NB
    fut : committed predicted issue current Ihat[t..t+H] (incl. this cycle's picks) minus its EWMA; last measured background
    cls : per-class issue counts of cycles t..t-KC+1 (lets the fit learn its own per-class current/voltage kernels)
    iq  : issue-queue lookahead - predicted-current kernels of the ready backlog (issue at t+1) and of uops whose wakeup
          is already scheduled for t+1 / t+2
    sen : quantized worst-tile sensor of cycles t-2..t-NS-1 minus its EWMA
    ev  : D$/I$ miss, refill, request, writeback flags over the last NE cycles (at the engine's trace cursor)
    evx : all 35 cache-pipeline activity flags over the last NX cycles
    phys: the existing design-kernel forecast (policy6 'meas' path without observer), 8 values"""
    S = streams(rec); n = len(S['Im']); Im = S['Im']; mu1 = lag(ewma(Im), 1); snap = S['snap']; Ih = snap[:, 0]; muh = ewma(Ih)
    k0 = 0 if ev_lag0 else 1
    cols = []; names = []
    def add(c, nm): cols.append(np.asarray(c, float)); names.append(nm)
    if 'cur' in groups:
        for k in range(1, NP + 1): add(lag(Im, k) - mu1, 'Im-%d' % k)
    if 'bank' in groups:
        for k in range(1, NB + 1): add(lag(S['Ib'], k), 'Ib-%d' % k)
    if 'fut' in groups:
        for k in range(H + 1): add(snap[:, k] - muh, 'Ih+%d' % k)
        add(lag(Im, 1) - lag(Ih, 1), 'bg-1')
    if 'cls' in groups:
        for c in range(NC):
            for d in range(KC): add(lag(S['ncls'][:, c], d), 'n%d-%d' % (c, d))
    if 'iq' in groups:
        for c in range(S['iq'].shape[1]): add(S['iq'][:, c], 'iq%d' % c)
    if 'sen' in groups:
        vb2 = lag(ewma(S['vs']), 2)
        for k in range(2, 2 + NS): add(lag(S['vs'], k) - vb2, 'vs-%d' % k)
    if 'ev' in groups or 'evx' in groups:
        src = rec['src'] + int(rec['a'])
        for nm, e_ in cache_sigs(str(rec['bench'])):
            x = e_[np.minimum(src, len(e_) - 1)]
            if 'ev' in groups and nm in MISS:
                for k in range(k0, NE + 1): add(lag(x, k), '%s-%d' % (nm, k))
            elif 'evx' in groups and nm not in MISS:
                for k in range(k0, NX + 1): add(lag(x, k), '%s-%d' % (nm, k))
    if 'phys' in groups:
        g = rec['g']; bg = lag(Im, 1) - lag(Ih, 1)
        for j in range(H):
            acc = np.zeros(n)
            for m_ in range(len(g)):
                tau = j - m_
                acc += g[m_] * ((snap[:, tau] + bg - mu1) if tau >= 0 else (lag(Im, -tau) - mu1))
            add(acc, 'phys%d' % j)
    if 'iqp' in groups:                               # predicted next selects pushed through the design kernel
        g = rec['g']; Sx = S['iq'][:, :H + 1] + S['iq'][:, H + 1:2 * H + 2]
        for j in range(H):
            add(sum(g[j - tau] * Sx[:, tau] for tau in range(j + 1)), 'iqp%d' % j)
    return np.column_stack(cols).astype(np.float32), names

def target(rec):
    v = rec['vw'] * 1e3; ew = ewma(v)
    return np.column_stack([lag(v, -j) - ew for j in range(H)])

# ------------------------------------------------------------------ fit / score -----------------------------------------------
W0 = 600
def score(f, y, j):
    n = len(y) - j; p = f[W0:n]; a = y[W0:n]
    r2 = np.corrcoef(a, p)[0, 1] ** 2
    ta = a <= np.percentile(a, 1); tp = p <= np.percentile(p, 1)
    return r2, (ta & tp).sum() / max(1, ta.sum())

class Acc(object):
    """Streaming ridge: accumulate sufficient statistics bench by bench (rows W0..n-H)."""
    def __init__(self): self.N = 0; self.sx = self.sy = self.xx = self.xy = None
    def add(self, X, Y):
        X = X[W0:len(X) - H].astype(np.float64); Y = Y[W0:len(Y) - H]
        if self.sx is None: self.sx = X.sum(0); self.sy = Y.sum(0); self.xx = X.T.dot(X); self.xy = X.T.dot(Y)
        else: self.sx += X.sum(0); self.sy += Y.sum(0); self.xx += X.T.dot(X); self.xy += X.T.dot(Y)
        self.N += len(X)
    def solve(self, lam=1e-3):
        N = self.N; mx = self.sx / N; my = self.sy / N
        C = self.xx / N - np.outer(mx, mx); c = self.xy / N - np.outer(mx, my)
        s = np.sqrt(np.maximum(np.diag(C), 0)); ok = s > 1e-9; s[~ok] = 1.0
        Cz = C / np.outer(s, s); cz = c / s[:, None]
        Cz[~ok, :] = 0; Cz[:, ~ok] = 0; cz[~ok] = 0
        Wz = np.linalg.solve(Cz + lam * np.eye(len(s)), cz); W = Wz / s[:, None]; b = my - mx.dot(W)
        return W, b

GROUPS = {'phys (=baseline w/o obs, refit)': ('phys',),
          'phys+sen': ('phys', 'sen'),
          'A cur+fut+sen+phys': ('cur', 'fut', 'sen', 'phys', 'bank'),
          'B A+iqp': ('cur', 'fut', 'sen', 'phys', 'bank', 'iqp'),
          'C B+iq(raw)': ('cur', 'fut', 'sen', 'phys', 'bank', 'iqp', 'iq'),
          'D B+ev(miss)': ('cur', 'fut', 'sen', 'phys', 'bank', 'iqp', 'ev'),
          'E D+evx(all cache)': ('cur', 'fut', 'sen', 'phys', 'bank', 'iqp', 'ev', 'evx'),
          'F E+cls': ('cur', 'fut', 'sen', 'phys', 'bank', 'iqp', 'ev', 'evx', 'cls')}

def evaluate(lam=1e-3, which=None, ev_lag0=False, save=None):
    res = {}
    for dn in ('decap', 'banks'):
        te = {bn: load(bn, dn) for bn in TEST}
        for bn, r in te.items():
            y = target(r); res[(dn, bn, '0 baseline meas+kern+obs.5')] = [score(r['flog'][:, j], y[:, j], j) for j in range(H)]
        for gn, gs in GROUPS.items():
            if which and gn[0] not in which: continue
            acc = Acc(); names = None
            for bn in TRAIN:
                r = load(bn, dn); X, names = feats(r, gs, ev_lag0); acc.add(X, target(r)); del X
            W, b = acc.solve(lam)
            for bn, r in te.items():
                X, _ = feats(r, gs, ev_lag0); y = target(r); f = b + X.dot(W)
                res[(dn, bn, gn)] = [score(f[:, j], y[:, j], j) for j in range(H)]
            if save and gn.startswith(save):
                np.savez(SCR + 'fc7_%s.npz' % dn, W=W, b=b, names=np.array(names), g=r['g'], groups=np.array(gs), lam=lam)
                sd = np.sqrt(np.maximum(np.diag(acc.xx / acc.N - np.outer(acc.sx / acc.N, acc.sx / acc.N)), 0))
                imp = np.abs(W) * sd[:, None]
                for j in (0, 2, 4):
                    gi = {}
                    for i, nm in enumerate(names):
                        key = re.sub(r'\d+$', '', re.sub(r'[-+]\d+$', '', nm))
                        gi[key] = gi.get(key, 0) + imp[i, j] ** 2
                    print(dn, 'j%d' % j, 'feature-group weight*std (rss):', ', '.join('%s %.2f' % (k, np.sqrt(v)) for k, v in sorted(gi.items(), key=lambda kv: -kv[1])[:12]), flush=True)
            print('fit', dn, gn, len(names), 'features', flush=True)
    for k, v in sorted(res.items()):
        print('%-6s %-7s %-32s ' % k + '  '.join('j%d R2 %.2f rec %.2f' % (j, v[j][0], v[j][1]) for j in (0, 2, 4)), flush=True)
    json.dump({'|'.join(k): v for k, v in res.items()}, open(SCR + 'eval_%s.json' % ('lag0' if ev_lag0 else 'lag1'), 'w'))
    return res

# ------------------------------------------------------------------ training / run-time ---------------------------------------
BEST = GROUPS['F E+cls']; LAM = 0.03
def wfile(dn, lag0, tag='train7'): return HERE + 'forecast7_w_%s_%s_%s.npz' % (dn, 'lag0' if lag0 else 'lag1', tag)

def train(dn, lag0=False, benches=None, lam=LAM, tag='train7'):
    """Fit the per-design model (group set BEST) on the given benches and save the weights next to this file."""
    import policy6 as P6
    acc = Acc(); names = None; benches = list(benches or TRAIN)
    for bn in benches:
        r = load(bn, dn); X, names = feats(r, BEST, lag0); acc.add(X, target(r)); del X
    W, b = acc.solve(lam); sig = [nm for nm, _ in cache_sigs(str(r['bench']))]
    hp = np.load(HERE + 'hpos_blocks.npy').sum(0)
    np.savez(wfile(dn, lag0, tag), W=W, b=b, names=np.array(names), groups=np.array(BEST), g=r['g'], hp=hp, lag0=int(lag0),
             sig=np.array(sig), benches=np.array(benches), lam=lam)
    return wfile(dn, lag0, tag)

class Forecaster7(object):
    """Run-time forecaster with the same information the policy6 engine has at cycle t.  Per cycle:
         fc.begin(Im_prev, Ib_prev, vs_mV, flags)     start of cycle t (before the select)
              Im_prev = measured current (I - Ibank, mA) of cycle t-1, Ib_prev = measured LDO-bank current (mA) of t-1
              (both None at t = 0); vs_mV = quantized worst-tile sensor sample of cycle t-2 in mV (None for t < 2);
              flags = cache activity flags at the engine's trace cursor for cycle t (vector in fc.sig order; 'ways'
              signals are counts). With lag0 weights the cycle-t flags are used, otherwise only t-1 and older.
         f = fc.forecast(fut, ncls, S12)               after the select: fut = Ihat[t..t+H] (committed predicted issue
              current incl. this cycle's picks, H+1 values), ncls = per-class count of uops picked at t (10), S12 = predicted
              current of the next two greedy selects (H+1 values: Ihat contribution of select t+1 shifted by 1 plus select
              t+2 shifted by 2; zeros if unavailable). Returns 8 values: forecast (mV) of v[t+j] - slow average, j = 0..7.
         f_pre = fc.forecast(fut_pre, zeros, S12) and f_pre + fc.E[cls].sum(0) over candidate uops: candidate evaluation
              (exact: the model is linear in ncls / fut).
         fc.end(Ih_t, ncls)                            end of cycle t: Ih_t = Ihat[t] after the picks."""
    def __init__(self, wf):
        z = np.load(wf); self.W = z['W']; self.b = z['b']; self.names = [str(s) for s in z['names']]; self.g = z['g']
        self.hp = z['hp']; self.lag0 = bool(z['lag0']); self.sig = [str(s) for s in z['sig']]; self.groups = [str(s) for s in z['groups']]
        self.idx = {n: i for i, n in enumerate(self.names)}; self.Kg = len(self.g)
        self.Im = collections.deque([0.0] * max(self.Kg, NP + 1), maxlen=max(self.Kg, NP + 1))     # [-1] = cycle t-1
        self.Ib = collections.deque([0.0] * NB, maxlen=NB); self.vs = collections.deque([0.0] * (NS + 2), maxlen=NS + 2)
        self.ncls = collections.deque([np.zeros(NC)] * KC, maxlen=KC); self.fl = collections.deque([np.zeros(len(self.sig))] * (NE + 1), maxlen=NE + 1)
        self.mu_m = None; self.mu1 = 0.0; self.vbar = None; self.vb2 = 0.0; self.muh = None; self.Ih_prev = 0.0; self.t = 0
        self.ev_sel = [(nm, k) for nm in self.sig for k in range(0 if self.lag0 else 1, (NE if nm in MISS else NX) + 1)]
        # per-class effect of one extra uop picked this cycle (linear model -> exact)
        self.E = np.zeros((NC, H))
        for c in range(NC):
            dx = np.zeros(len(self.names)); hc = self.hp[c]
            if 'n%d-0' % c in self.idx: dx[self.idx['n%d-0' % c]] = 1.0
            for k in range(H + 1):
                if 'Ih+%d' % k in self.idx: dx[self.idx['Ih+%d' % k]] = hc[k]
            for j in range(H):
                if 'phys%d' % j in self.idx: dx[self.idx['phys%d' % j]] = sum(self.g[j - tau] * hc[tau] for tau in range(j + 1))
            self.E[c] = dx.dot(self.W)
    def begin(self, Im_prev, Ib_prev, vs_mV, flags):
        if Im_prev is not None:
            self.Im.append(float(Im_prev)); self.Ib.append(float(Ib_prev))
            self.mu_m = float(Im_prev) if self.mu_m is None else self.mu_m + (Im_prev - self.mu_m) / 256.0
            self.mu1 = self.mu_m
        if vs_mV is not None:
            self.vs.append(float(vs_mV)); self.vbar = float(vs_mV) if self.vbar is None else self.vbar + (vs_mV - self.vbar) / 256.0
            self.vb2 = self.vbar
        else: self.vs.append(0.0)
        self.fl.append(np.asarray(flags, float))
    def _x(self, fut, ncls, S12):
        d = {}; t = self.t; Im = self.Im; L = len(Im); mu1 = self.mu1
        imk = lambda k: Im[L - k] if k <= t else 0.0                          # Im[t-k], zero before start
        for k in range(1, NP + 1): d['Im-%d' % k] = imk(k) - mu1
        for k in range(1, NB + 1): d['Ib-%d' % k] = self.Ib[NB - k] if k <= t else 0.0
        muh = float(fut[0]) if self.muh is None else self.muh + (fut[0] - self.muh) / 256.0
        for k in range(H + 1): d['Ih+%d' % k] = fut[k] - muh
        bg = (imk(1) - self.Ih_prev) if t >= 1 else 0.0; d['bg-1'] = bg
        nh = list(self.ncls)
        for c in range(NC):
            d['n%d-0' % c] = ncls[c]
            for dd in range(1, KC): d['n%d-%d' % (c, dd)] = nh[KC - dd][c] if dd <= t else 0.0
        for k in range(2, 2 + NS): d['vs-%d' % k] = (self.vs[NS + 2 - 1 - (k - 2)] if k <= t else 0.0) - self.vb2
        fl = self.fl; si = {nm: i for i, nm in enumerate(self.sig)}
        for nm, k in self.ev_sel: d['%s-%d' % (nm, k)] = fl[NE - k][si[nm]] if k <= t else 0.0
        g = self.g
        for j in range(H):
            acc = 0.0
            for m_ in range(self.Kg):
                tau = j - m_
                acc += g[m_] * ((fut[tau] + bg - mu1) if tau >= 0 else (imk(-tau) - mu1))
            d['phys%d' % j] = acc
        for j in range(H): d['iqp%d' % j] = sum(g[j - tau] * S12[tau] for tau in range(j + 1))
        return np.array([d[n] for n in self.names])
    def forecast(self, fut, ncls, S12):
        return self.b + self._x(fut, ncls, S12).dot(self.W)
    def end(self, Ih_t, ncls):
        self.muh = float(Ih_t) if self.muh is None else self.muh + (Ih_t - self.muh) / 256.0
        self.Ih_prev = float(Ih_t); self.ncls.append(np.asarray(ncls, float)); self.t += 1

def verify(bench, dn, lag0, n=3000):
    """Step Forecaster7 through a recorded run and compare with the offline feature path (must agree to rounding)."""
    r = load(bench, dn); wf = wfile(dn, lag0); fc = Forecaster7(wf); z = np.load(wf)
    X, names = feats(r, BEST, lag0); assert names == fc.names
    S = streams(r); sig = dict(cache_sigs(bench)); src = r['src'] + int(r['a'])
    F_off = z['b'] + X[:n].dot(z['W']); err = 0.0; errE = 0.0
    S12 = S['iq'][:, :H + 1] + S['iq'][:, H + 1:2 * H + 2]
    for t in range(n):
        flags = [sig[nm][min(src[t], len(sig[nm]) - 1)] for nm in fc.sig]
        fc.begin(S['Im'][t - 1] if t else None, S['Ib'][t - 1] if t else None, S['vs'][t - 2] if t >= 2 else None, flags)
        f = fc.forecast(S['snap'][t], S['ncls'][t], S12[t]); err = max(err, np.abs(f - F_off[t]).max())
        # candidate path: pre-pick state + per-class effect vectors
        pre = S['snap'][t] - sum(S['ncls'][t][c] * fc.hp[c][:H + 1] for c in range(NC))
        f2 = fc.forecast(pre, np.zeros(NC), S12[t]) + S['ncls'][t].dot(fc.E); errE = max(errE, np.abs(f2 - f).max())
        fc.end(S['snap'][t][0], S['ncls'][t])
    return err, errE

if __name__ == '__main__':
    if not os.path.isdir(SCR): os.makedirs(SCR)
    if sys.argv[1] == 'rec':
        import multiprocessing as mp
        jobs = sorted(list(TEST.items()) + list(TRAIN.items()), key=lambda kv: -((kv[1][1] or 60000) - kv[1][0]))
        with mp.get_context('spawn').Pool(4) as P: P.map(_rec_job, jobs, chunksize=1)
    elif sys.argv[1] == 'eval':     # eval [lam] [group letters] [lag0] [save-prefix]
        evaluate(float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3, sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != '-' else None,
                 len(sys.argv) > 4 and sys.argv[4] == 'lag0', sys.argv[5] if len(sys.argv) > 5 else None)
    elif sys.argv[1] == 'train':    # train: weights for both designs, lag1 and lag0, on TRAIN (held-out) and on all 10 benches
        for dn in ('decap', 'banks'):
            for l0 in (False, True):
                print(train(dn, l0), train(dn, l0, list(TRAIN) + list(TEST), tag='all10'), flush=True)
    elif sys.argv[1] == 'verify':
        for dn in ('decap', 'banks'):
            for l0 in (False, True): print('verify', dn, 'lag0' if l0 else 'lag1', 'max |runtime - offline| = %.2e mV, candidate path err %.2e mV' % verify('towers', dn, l0), flush=True)
