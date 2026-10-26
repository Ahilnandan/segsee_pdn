"""Block-level per-uop current attribution for the closed-loop tiled simulator (exact in the baseline).

Issue-coupled blocks (LACPo execution blocks + D$ hit path of the cache model): their current above a per-block floor is
shared among the uops granted in the last K cycles, in proportion to non-negative per-class kernels (NNLS, fit on FIT
benches). Uop i owns A[i, b, k] mA of block b at cycle t_i + k. Everything else is per-tile background BG[n, t].
Identity: BG[n, :] + sum_b W_b * sum_i A[i, b, n - t_i] == tile current[n, :] for the real schedule (checked)."""
from __future__ import print_function
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.optimize import nnls
import kernels as KM, tilecurrents as TC
B = KM.B; HERE = KM.HERE; K = KM.K; NC = len(KM.CLS)
BLOCKS = ['alu0', 'alu1_csr', 'iregfile', 'iregread', 'int_iq', 'mem_iq', 'lsu', 'fp', 'dc_hit']
WB = np.stack([TC.W[b] if b != 'dc_hit' else TC.W['dcache_data'] for b in BLOCKS])   # (nb, tiles)

def issue_blocks(w):
    """per-block issue-coupled currents (mA) and full tile currents, baseline schedule"""
    blk, n = TC.block_currents(w)
    cf = json.load(open(HERE + 'cache_model_coef.json'))['dc']
    E = {k: v.astype(float)[:n] for k, v in np.load(B + w + '_cache.npz').items()}
    from cache_model import feats_dc
    F = feats_dc(E); hit = sum(cf.get(k, 0) * F[k] for k in ('rd', 'wr', 'rd_wr', 'rd_data_hd', 'wr_next', 'meta_wr')) / 1.1
    IC = np.stack([blk[b] if b != 'dc_hit' else hit for b in BLOCKS], 1)                # (n, nb) mA
    return IC, n

def fit(nsub=150000, seed=0):
    rng = np.random.RandomState(seed); Xs, Ys = [], []
    for w in KM.FIT:
        IC, n = issue_blocks(w); N = KM.counts(w, n); s = KM.exec_start(IC.sum(1, keepdims=True) * 0 + KM.region_currents(w))
        X = KM.design(N[s:]); Y = IC[s:] - np.percentile(IC[s:], 5, axis=0)
        idx = rng.choice(len(X), min(nsub // len(KM.FIT), len(X)), replace=False); Xs.append(X[idx]); Ys.append(Y[idx])
    X = np.vstack(Xs).astype(float); Y = np.vstack(Ys); h = np.zeros((len(BLOCKS), NC, K))
    for b in range(len(BLOCKS)):
        if Y[:, b].std() < 1e-6: continue
        coef, _ = nnls(X, np.maximum(Y[:, b], 0)); h[b] = coef.reshape(NC, K)
    return h + 1e-3

def build(w, h):
    IC, n = issue_blocks(w); It = TC.tile_currents(w)[:n].astype(float)
    z = np.load(B + w + '_issue.npz'); keep = z['t'] < n
    inst = {k: z[k][keep] for k in z.files if k not in ('ncyc', 'nvalid', 'nreq')}
    cls = KM.uop_class(inst['fu'], inst['pdst'], inst['dtype']); t = inst['t']; m = len(t); nb = len(BLOCKS)
    F = np.percentile(IC, 5, axis=0); ex = np.maximum(IC - F, 0.0)
    Wk = h[:, cls, :].transpose(1, 0, 2)                                                  # (m, nb, K)
    S = np.zeros((n + K, nb))
    for k in range(K): np.add.at(S, t + k, Wk[:, :, k])
    A = np.zeros((m, nb, K), np.float32); attrib = np.zeros((n, nb))
    for k in range(K):
        nn = t + k; ok = nn < n
        A[ok, :, k] = ex[nn[ok]] * Wk[ok, :, k] / np.maximum(S[nn[ok]], 1e-12)
        np.add.at(attrib, nn[ok], A[ok, :, k])
    BG = It - attrib.dot(WB)                                                               # (n, tiles)
    R = BG + attrib.dot(WB); err = np.abs(R - It).max()
    return dict(inst=inst, cls=cls, A=A, BG=BG.astype(np.float32), err=err, movable=attrib.sum() / It.sum())

if __name__ == '__main__':
    h = fit(); np.save(HERE + 'hpos_blocks.npy', h)
    print('NNLS kernel charge per uop (pC) per block:'); print('       ' + ''.join('%9s' % b[:8] for b in BLOCKS))
    for c, nm in enumerate(KM.CLS): print('  %-5s' % nm + ''.join('%9.1f' % v for v in h[:, c, :].sum(1) * 3))
    for w in (sys.argv[1:] or KM.WL):
        d = build(w, h)
        np.savez_compressed(B + w + '_attr5.npz', A=d['A'], BG=d['BG'], cls=d['cls'], **d['inst'])
        print('%-10s instances %7d identity err %.1e mA movable share %.3f' % (w, len(d['cls']), d['err'], d['movable']))
