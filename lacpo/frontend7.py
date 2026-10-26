"""B6: frontend supply-current model vs fetch activity, and the fetch-throttle lever for the closed-loop engine.

Frontend = LACPo blocks bpd (gshare counter table), btb (BTB+BIM), fetch (fetch_controller: FTQ, fetch buffer, branch
checker, RVC expand) + the calibrated I$ event model (cache_model.py). Decode is core-side and reported, not modelled.
Core = BOOM v2 MediumBoom (chipyard 717c4658): WithRVC doubles fetchWidth 2->4 halfword slots = 8 B/packet,
16-entry fetch buffer (enq needs room for a full packet), decode/dispatch 2/cycle, issue 4/cycle (mem1, int2, fp1).

Measured I$ pipeline relations (exact on traces): s1[t+1] = s0[t];  s2[t+1] = s1[t] & ~kill[t];  kill ~= 1 when s1 = 0;
data-out toggles only when s1 valid; tag/way reads ~= s0.  A fetch gate at s0 for cycles [t, t+g) therefore removes
s0/tag/way activity on [t,t+g), s1 + data-out on [t+1,t+g+1), s2 (packet into fetch controller / fetch buffer) on
[t+2,t+g+2).  Model: per-block linear ridge on these signals with lags 0..4 (+ issue count as non-gated covariate),
fit leave-one-bench-out.  HTIF load phase stripped (same detector as tilecurrents.py).

Run:  OPENBLAS_NUM_THREADS=1 <lacpo37 python> frontend7.py     -> prints report, writes frontend7_model.npz
"""
from __future__ import print_function
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, '/home/ahilnandan/Documents/PMIC/pdn_work/lacpo')
import tilecurrents as TC
B, HERE = TC.B, TC.HERE
WL = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
FE = ['icache', 'bpd', 'btb', 'fetch']
SIG = ['s0', 's1', 's2', 'kill', 'tag', 'ways', 'hd', 'refill', 'wen', 'miss', 'dbus', 'iss']
GATED = ['s0', 's1', 's2', 'kill', 'tag', 'ways', 'hd']       # signals a fetch gate changes
LAGS = list(range(0, 5))                                      # current[t] sees signal[t-lag] (causal; leads gave no gain)
FB_ENTRIES, DEC_W, FE_LAT = 16, 2, 3                          # fetch buffer, decode width, s0 -> fetch-buffer latency
RIDGE = 1e-3

def start_of_exec(tot):
    n = len(tot)
    return int(np.argmax(pd.Series(tot).rolling(100).mean().fillna(0).values * (np.arange(n) > 3000) > 150))

def base_signals(E):
    s1 = E['ic_s1'] > 0
    return {'s0': (E['ic_s0'] > 0) * 1.0, 's1': s1 * 1.0, 's2': (E['ic_s2'] > 0) * 1.0, 'kill': (E['ic_kill'] > 0) * 1.0,
            'tag': (E['ic_tag'] > 0) * 1.0, 'ways': (E['ic_w0'] + E['ic_w1'] + E['ic_w2'] + E['ic_w3']) / 4.0,
            'hd': (E['ic_d0_hd'] + E['ic_d1_hd'] + E['ic_d2_hd'] + E['ic_d3_hd']) / 256.0, 'refill': E['ic_refill'] * 1.0,
            'wen': E['ic_wen'] * 1.0, 'miss': E['ic_miss'] * 1.0,
            'dbus': np.maximum(E['dc_dbus_hd'], E['ic_dbus_hd']) / 64.0, 'iss': E['iss'] * 1.0}

def design(S, a=None, b=None):
    """rows a..b of the lagged design matrix (const first)."""
    n = len(S['s0']); a = 0 if a is None else a; b = n if b is None else b
    cols = [np.ones(b - a)]
    for k in SIG:
        x = S[k]
        for L in LAGS:
            idx = np.arange(a, b) - L; ok = (idx >= 0) & (idx < n); c = np.zeros(b - a); c[ok] = x[idx[ok]]; cols.append(c)
    return np.stack(cols, 1)

def load(w):
    blk, n = TC.block_currents(w)
    tileI = np.load(B + w + '_tileI.npy', mmap_mode='r'); n = min(n, len(tileI))
    tot = np.asarray(tileI[:n]).sum(1).astype(float); s = start_of_exec(tot)
    E = {k: v[:n].astype(float) for k, v in np.load(B + w + '_cache.npz').items()}
    z = np.load(B + w + '_attr5.npz'); t = z['t']; t = t[(t >= s) & (t < n)] - s
    BG = z['BG'][s:n].astype(float)
    E = {k: v[s:] for k, v in E.items()}; E['iss'] = np.bincount(t, minlength=n - s)[:n - s]
    Y = np.stack([blk[b][s:n] for b in FE], 1)
    fe_tiles = np.where(sum(TC.W[b] for b in FE) > 0)[0]
    return dict(w=w, s=s, S=base_signals(E), Y=Y, dec=blk['decode'][s:n], tot=tot[s:], icpl=tot[s:] - BG.sum(1),
                fe_tiles=fe_tiles, tile_fe=np.asarray(tileI[s:n])[:, fe_tiles].astype(float))

def fit(XtX, XtY):
    A = XtX + RIDGE * np.trace(XtX) / len(XtX) * np.eye(len(XtX)); A[0, 0] = XtX[0, 0]
    return np.linalg.solve(A, XtY)

def predict(C, S, a=None, b=None):
    return design(S, a, b).dot(C)

def gate_signals(S, t, g):
    """copy of signals on window [t-2, t+g+8) with fetch gated at s0 for [t, t+g)."""
    a, b = max(0, t - 6), min(len(S['s0']), t + g + 10); G = {k: S[k][a:b].copy() for k in SIG}; o = t - a
    for k in ('s0', 'tag', 'ways'): G[k][o:o + g] = 0
    for k in ('s1', 'hd'): G[k][o + 1:o + g + 1] = 0
    G['kill'][o + 1:o + g + 1] = 1
    G['s2'][o + 2:o + g + 2] = 0
    return G, a, b

def delta_current(C, S, t, g):
    """per-block current removed (mA, >0 = reduction) on cycles t-2 .. t+g+5 by gating fetch for g cycles at t."""
    G, a, b = gate_signals(S, t, g)
    O = {k: S[k][a:b] for k in SIG}
    d = (predict(C, O) - predict(C, G))[:, :len(FE)]; o = t - a      # drop the fitted-total column
    return d[max(0, o - 2):o + g + 6]

def slack_cycles(S, t, H=8):
    """backpressure-stall cycles (fetch idle, no I$ miss/refill nearby) in [t, t+H): fetch-buffer-full slack the gate can
    borrow without starving decode."""
    s0 = S['s0'][t:t + H]; m = S['miss'][max(0, t - 20):t + H + 20].sum() + S['refill'][max(0, t - 20):t + H + 20].sum()
    return 0 if m > 0 else int((s0 == 0).sum())

class FrontendThrottle(object):
    """Closed-loop engine interface.  ft = FrontendThrottle('qsort');  r = ft.gate(t, g)
    t is a cycle index into the stripped trace (same origin as tilecurrents start_of_exec / BG rows after load phase).
    Returns dict: dI_total (mA removed per cycle, array from cycle t-2 on), dI_tiles (cycles x 64 mA, via floorplan),
    g_eff (cycles fetch was actually active in the gate), delay (cycles by which all later instructions become ready
    later: g_eff - slack, >= 0), slack, occupancy_full (bool)."""
    def __init__(self, w, model=HERE + 'frontend7_model.npz'):
        z = np.load(model); self.C = z['C_all'] if w not in list(z['benches']) else z['C_lobo'][list(z['benches']).index(w)]
        self.d = load(w); self.S = self.d['S']
        self.WB = np.stack([TC.W[b] for b in FE], 1)               # (64, 4)
    def baseline(self, t):
        return float(predict(self.C, self.S, t, t + 1)[0, :len(FE)].sum())
    def gate(self, t, g):
        dI = delta_current(self.C, self.S, t, g)                   # (cycles, 4 blocks)
        g_eff = int(self.S['s0'][t:t + g].sum()); sl = slack_cycles(self.S, t + g)
        return dict(t0=t - 2, dI_total=dI.sum(1), dI_tiles=dI.dot(self.WB.T), g_eff=g_eff, slack=sl,
                    delay=max(0, g_eff - sl), occupancy_full=sl > 0)

def main():
    D = {}; acc = {}
    for w in WL:
        d = load(w); D[w] = d; X = design(d['S']); Yt = np.c_[d['Y'], d['Y'].sum(1)]
        acc[w] = (X.T.dot(X), X.T.dot(Yt)); del X
        print('loaded %-10s exec cycles %6d (load phase %d stripped)' % (w, len(d['Y']), d['s']), file=sys.stderr)
    nf = len(acc[WL[0]][0])
    tot_XtX = sum(acc[w][0] for w in WL); tot_XtY = sum(acc[w][1] for w in WL)
    C_all = fit(tot_XtX, tot_XtY)
    C_lobo = np.stack([fit(tot_XtX - acc[w][0], tot_XtY - acc[w][1]) for w in WL])
    np.savez(HERE + 'frontend7_model.npz', C_all=C_all, C_lobo=C_lobo, benches=np.array(WL), sig=np.array(SIG), lags=np.array(LAGS))

    d0 = D['qsort']
    print('== 1. Frontend blocks / tiles ==')
    print('frontend tiles (row*8+col): %s  (%d of 64)' % (list(d0['fe_tiles']), len(d0['fe_tiles'])))
    print('%-10s %7s %7s %7s %7s %7s | %6s %6s %6s | %7s %7s %7s' % ('bench', 'I$', 'bpd', 'btb', 'fetch', 'FEtot', 'FEstd', 'FEp1', 'FEp99',
                                                                    'tile', 'icplP50', 'icplP99'))
    for w in WL:
        d = D[w]; m = d['Y'].mean(0); F = d['Y'].sum(1)
        print('%-10s %7.1f %7.1f %7.1f %7.1f %7.1f | %6.1f %6.1f %6.1f | %7.1f %7.1f %7.1f' % (w, m[0], m[1], m[2], m[3], F.mean(), F.std(),
              np.percentile(F, 1), np.percentile(F, 99), d['tot'].mean(), np.percentile(d['icpl'], 50), np.percentile(d['icpl'], 99)))
    # idle level: load phase (fetch off, s0 = 0) per block
    blk, n = TC.block_currents('qsort'); idle = [np.median(blk[b][5000:40000]) for b in FE]
    print('idle (fetch off, HTIF load phase) mA: ' + '  '.join('%s %.1f' % (b, v) for b, v in zip(FE, idle)) + '  sum %.1f' % sum(idle))

    print('\n== 2. Held-out fit (leave-one-bench-out), total frontend current ==')
    print('%-10s %6s %6s %6s %6s %7s | %6s %6s | per-block R2 %s' % ('bench', 'R2', 'rms', 'bias', 'std', 'stdhat', 'dIp99', 'hat', ' '.join(FE)))
    r2s = []
    for i, w in enumerate(WL):
        d = D[w]; P = predict(C_lobo[i], d['S']); Yt = np.c_[d['Y'], d['Y'].sum(1)]; e = P - Yt
        r2 = 1 - (e ** 2).mean(0) / Yt.var(0); r2s.append(r2)
        print('%-10s %6.3f %6.2f %+6.2f %6.2f %7.2f | %6.2f %6.2f | %s' % (w, r2[4], np.sqrt((e[:, 4] ** 2).mean()), e[:, 4].mean(), Yt[:, 4].std(),
              P[:, 4].std(), np.percentile(np.abs(np.diff(Yt[:, 4])), 99), np.percentile(np.abs(np.diff(P[:, 4])), 99),
              ' '.join('%.2f' % x for x in r2[:4])))
    r2s = np.array(r2s); print('median held-out R2 total %.3f; blocks %s' % (np.median(r2s[:, 4]), ' '.join('%s %.2f' % (b, x) for b, x in zip(FE, np.median(r2s[:, :4], 0)))))

    print('\n== 3a. Gate lever (model, LOBO coefficients), random fetch-active cycles and heavy-burst cycles ==')
    rng = np.random.RandomState(0); rows = {}
    for i, w in enumerate(WL):
        d = D[w]; S = d['S']; n = len(S['s0']); hi = np.percentile(d['icpl'], 95)
        cand = np.where(S['s0'][10:n - 20] > 0)[0] + 10; burst = np.where((d['icpl'][10:n - 20] > hi) & (S['s0'][10:n - 20] > 0))[0] + 10
        for tag, pool in (('rand', cand), ('burst', burst)):
            for g in (1, 2, 3, 4):
                T = rng.choice(pool, min(400, len(pool)), replace=False)
                prof = np.array([delta_current(C_lobo[i], S, t, g).sum(1) for t in T])     # (events, cycles from t-2)
                sl = np.array([slack_cycles(S, t + g) for t in T]); ge = np.array([S['s0'][t:t + g].sum() for t in T])
                rows.setdefault((tag, g), []).append((prof.mean(0), np.maximum(0, ge - sl).mean(), (sl > 0).mean(), w))
    for (tag, g), R in sorted(rows.items()):
        P = np.mean([r[0] for r in R], 0); mx = P.max(); edge = np.abs(np.diff(np.r_[0, P])).max()
        q = P.sum()
        print('%-5s g=%d  dI per cycle from t-2: %s  peak %.1f mA  edge %.1f mA/cyc  charge %.0f mA*cyc (%.2f nC)  delay %.2f cyc  FB-full %.2f'
              % (tag, g, ' '.join('%.1f' % x for x in P[:g + 6]), mx, edge, q, q * 3e-3, np.mean([r[1] for r in R]), np.mean([r[2] for r in R])))
    print('per-bench (burst, g=4): ' + '  '.join('%s %.1fmA/d%.1f/full%.2f' % (r[3], r[0].max(), r[1], r[2]) for r in rows[('burst', 4)]))

    print('\n== 3b. Model-independent check: natural clean fetch stalls (s0 off L cycles, no I$ miss +-20) ==')
    print('actual frontend current change vs pre-level (mean of t-4..t-1), cycles t-2..t+L+5; model (LOBO) in brackets')
    for L in (1, 2, 3, 4):
        A, M, cnt = [], [], 0
        for i, w in enumerate(WL):
            d = D[w]; S = d['S']; s0 = S['s0']; n = len(s0)
            z = np.diff(np.r_[1, s0, 1]); st = np.where(z == -1)[0]; en = np.where(z == 1)[0]; ln = en - st
            mi = np.convolve(S['miss'] + S['refill'], np.ones(41), 'same') > 0
            for a in st[(ln == L)]:
                if a < 10 or a + L + 10 >= n or mi[a]: continue
                F = d['Y'][a - 4:a + L + 6].sum(1); A.append(F[:4].mean() - F[2:])
                Pm = predict(C_lobo[i], S, a - 4, a + L + 6)[:, :len(FE)].sum(1); M.append(Pm[:4].mean() - Pm[2:]); cnt += 1
        if cnt:
            A, M = np.mean(A, 0), np.mean(M, 0)
            print('L=%d n=%5d  %s' % (L, cnt, ' '.join('%.1f[%.1f]' % (x, y) for x, y in zip(A, M))))

    print('\n== 3c. Instruction supply ==')
    print('%-10s %6s %6s %6s %7s %7s %8s' % ('bench', 'IPCiss', 's0act', 'ins/pk', 'bpStall', 'cover@1', 'cover@2'))
    for w in WL:
        S = D[w]['S']; iss = S['iss']; ipk = iss.sum() / max(1, S['s2'].sum())
        mi = np.convolve(S['miss'] + S['refill'], np.ones(41), 'same') > 0; bp = ((S['s0'] == 0) & ~mi).mean()
        print('%-10s %6.2f %6.3f %6.2f %7.3f %7.1f %8.1f' % (w, iss.mean(), S['s0'].mean(), ipk, bp, FB_ENTRIES / 1.0, FB_ENTRIES / 2.0))
    print('cover@x = cycles a FULL 16-entry fetch buffer feeds decode at x inst/cycle (+%d-cycle in-flight s1/s2 packets).' % (FE_LAT - 1))

    print('\n== 4. Lever vs issue-coupled current, and interface demo ==')
    g4 = np.mean([r[0] for r in rows[('burst', 4)]], 0)
    for w in WL:
        d = D[w]; print('%-10s issue-coupled p50 %5.1f p99 %5.1f p99.9 %5.1f max %5.1f | tile mean %5.1f p99 %5.1f | g=4 peak cut %.1f mA = %.0f%% of icpl p99'
              % (w, np.percentile(d['icpl'], 50), np.percentile(d['icpl'], 99), np.percentile(d['icpl'], 99.9), d['icpl'].max(), d['tot'].mean(),
                 np.percentile(d['tot'], 99), g4.max(), 100 * g4.max() / np.percentile(d['icpl'], 99)))
    ft = FrontendThrottle('mm'); t = int(np.where(ft.S['s0'][1000:] > 0)[0][0]) + 1000
    for g in (1, 4):
        r = ft.gate(t, g)
        print('demo mm t=%d g=%d: baseline FE %.1f mA, dI_total %s, g_eff %d slack %d delay %d, hottest tile cut %.2f mA (tile %d)'
              % (t, g, ft.baseline(t), ' '.join('%.1f' % x for x in r['dI_total']), r['g_eff'], r['slack'], r['delay'],
                 r['dI_tiles'].max(), int(np.argmax(r['dI_tiles'].max(0)))))

if __name__ == '__main__':
    main()
