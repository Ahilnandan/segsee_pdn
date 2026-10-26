"""Per-uop current attribution for replay (exact in the baseline, moves with the uop when it is rescheduled).

Issue-coupled region current IC[n, r] (LACPo execution blocks + D$ hit path, see kernels.region_currents) above a
per-region floor F_r is shared among the uops granted in cycles n-K+1..n, in proportion to non-negative per-class
kernel shapes h+[r, c, k] (NNLS fit). Uop i then owns A[i, r, k] mA at cycle t_i + k. Everything else
(background + floor + unattributed excess when nothing issued) is BG[n, r], which a replay keeps and time-stretches.
Identity: BG[n] + sum_i A[i, :, n - t_i] == tile region current[n] for the real schedule."""
from __future__ import print_function
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.optimize import nnls
import kernels as KM
B = KM.B; HERE = KM.HERE; K = KM.K; NC = len(KM.CLS); NR = 4

def fit_hpos(nsub=120000, seed=0):
    rng = np.random.RandomState(seed); Xs, Ys = [], []
    for w in KM.FIT:
        IC, N, s, BG = KM.load(w); X = KM.design(N); idx = rng.choice(len(X), min(nsub // len(KM.FIT), len(X)), replace=False)
        Xs.append(X[idx]); Ys.append(IC[idx] - np.percentile(IC, 5, axis=0))
    X = np.vstack(Xs).astype(float); Y = np.vstack(Ys); h = np.zeros((NR, NC, K))
    for r in range(NR):
        if Y[:, r].std() < 1e-6: continue
        coef, _ = nnls(X, np.maximum(Y[:, r], 0)); h[r] = coef.reshape(NC, K)
    h += 1e-3                                                   # every class owns a sliver of every region (no zero-weight holes)
    return h

def build(w, h):
    IC, BGo = KM.region_currents(w, split=True); n = len(IC)
    z = np.load(B + w + '_issue.npz'); t = z['t']; keep = t < n
    inst = {k: z[k][keep] for k in z.files if k not in ('ncyc', 'nvalid', 'nreq')}
    cls = KM.uop_class(inst['fu'], inst['pdst'], inst['dtype']); t = inst['t']; m = len(t)
    F = np.percentile(IC, 5, axis=0); ex = np.maximum(IC - F, 0.0); BG = BGo + np.minimum(IC, F)
    S = np.zeros((n + K, NR)); W = h[:, cls, :].transpose(1, 0, 2)   # (m, R, K)
    for k in range(K): np.add.at(S, t + k, W[:, :, k])
    A = np.zeros((m, NR, K), np.float32)
    for k in range(K):
        nn = t + k; ok = nn < n
        A[ok, :, k] = (ex[nn[ok]] * W[ok, :, k] / np.maximum(S[nn[ok]], 1e-12))
    covered = S[:n] > 1e-9; BG = BG + np.where(covered, 0.0, ex)      # excess with no issuing uop stays background
    # check identity
    R = BG.copy()
    for k in range(K):
        nn = t + k; ok = nn < n; np.add.at(R, nn[ok], A[ok, :, k])
    err = np.abs(R - (IC + BGo)).max()
    return dict(inst=inst, cls=cls, A=A, BG=BG.astype(np.float32), n=n, err=err)

if __name__ == '__main__':
    h = fit_hpos(); np.save(HERE + 'hpos.npy', h)
    print('NNLS kernel charge per uop (pC), regions R0..R3:')
    for c, nm in enumerate(KM.CLS): print('  %-5s' % nm, ' '.join('%7.1f' % v for v in h[:, c, :].sum(1) * 3))
    for w in KM.WL:
        d = build(w, h)
        np.savez_compressed(B + w + '_attr.npz', A=d['A'], BG=d['BG'], cls=d['cls'], **d['inst'])
        movable = sum(d['A'][:, r, :].sum() for r in range(4)) / (d['BG'].sum() + d['A'].sum())
        print('%-10s instances %7d  identity err %.2e mA  movable share of tile charge %.2f' % (w, len(d['cls']), d['err'], movable))
