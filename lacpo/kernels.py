"""Per-region current = unschedulable background + sum over issued uop classes of a per-class current kernel.

Regions (mA at 1.1 V): R0 frontend+I$, R1 integer core, R2 LSU+D$+PTW, R3 FP. Region currents come from LACPo
per-block predictions + the calibrated cache model (tile_P_t.csv). Kernels h[r, c, k] (k = 0..K-1 cycles after grant)
are ridge-fitted on the execution phase of FIT benches and scored on HELD-OUT benches (R^2 of the region current and of
its fast part). Background B_r[n] = I_r[n] - (h * N)[n] is what replay keeps (time-stretched) when issue timing changes."""
from __future__ import print_function
import json, sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
HERE = '/home/ahilnandan/Documents/PMIC/pdn_work/lacpo/'
WL = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
FIT = ['qsort', 'mm', 'spmv', 'memcpy', 'towers']; TEST = [w for w in WL if w not in FIT]
REG = ['R0_frontend_I$', 'R1_int_core', 'R2_lsu_D$', 'R3_fp']
CLS = ['ALU', 'BRU', 'LD', 'ST', 'MUL', 'DIV', 'CSR', 'FPU', 'FDV', 'CVT']
K = 16; LAM = 30.0
FU_ALU, FU_BRU, FU_MEM, FU_MUL, FU_DIV, FU_CSR, FU_FPU, FU_FDV, FU_I2F, FU_F2I = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512

def uop_class(fu, pdst, dtype):
    """vectorised class index per instance"""
    c = np.full(len(fu), 0)
    c[(fu & FU_BRU) > 0] = 1
    mem = (fu & FU_MEM) > 0; c[mem & (dtype != 2)] = 2; c[mem & (dtype == 2)] = 3
    c[(fu & FU_MUL) > 0] = 4; c[(fu & FU_DIV) > 0] = 5; c[(fu & FU_CSR) > 0] = 6
    c[(fu & FU_FPU) > 0] = 7; c[(fu & FU_FDV) > 0] = 8; c[(fu & (FU_I2F | FU_F2I)) > 0] = 9
    return c

ISSUE_BLOCKS = ('ALUExeUnit', 'iregfile', 'iregister_read', 'int_issue', 'mem_issue', 'lsu', 'fp_pipeline', 'csr_exe_unit')
def region_currents(w, split=False):
    """region currents (mA). split=True -> (issue-coupled part, background part): issue-coupled = LACPo execution blocks
    (ALU/CSR exe, int RF + read, int/mem IQ, LSU, FP pipeline) + D$ hit-path terms of the cache model (array read/write,
    read-data toggles, meta write); background = frontend, decode, rename, ROB, glue, I$, D$ floor/refill/writeback, PTW."""
    d = pd.read_csv(B + w + '_out/composed_P_t.csv'); t = pd.read_csv(B + w + '_out/tile_P_t.csv')
    n = min(len(d), len(t)); col = lambda s: d[[c for c in d.columns if s(c)]].values[:n].sum(1)
    blk = [c for c in d.columns if c.startswith('core') or c.startswith('frontend')]
    reg = lambda c: 0 if (c.startswith('frontend') or 'decode' in c) else (3 if 'fp_' in c else (2 if ('lsu' in c or 'mem_issue' in c) else 1))
    ic = np.zeros((n, 4)); bg = np.zeros((n, 4))
    for c in blk:
        (ic if any(b in c for b in ISSUE_BLOCKS) else bg)[:, reg(c)] += d[c].values[:n]
    import json as _j
    cf = _j.load(open(HERE + 'cache_model_coef.json'))['dc']
    E = {k: v.astype(float)[:n] for k, v in np.load(B + w + '_cache.npz').items()}
    from cache_model import feats_dc
    F = feats_dc(E); hit = sum(cf.get(k, 0) * F[k] for k in ('rd', 'wr', 'rd_wr', 'rd_data_hd', 'wr_next', 'meta_wr'))
    ic[:, 2] += hit; bg[:, 2] += t.P_dcache_mW.values[:n] - hit + 7.24
    bg[:, 0] += t.P_icache_mW.values[:n]; bg[:, 1] += 1.5
    assert np.allclose((ic + bg).sum(1), t.P_tile_mW.values[:n], atol=1e-4)
    if split: return ic / 1.1, bg / 1.1
    return (ic + bg) / 1.1

def counts(w, n):
    z = np.load(B + w + '_issue.npz'); c = uop_class(z['fu'], z['pdst'], z['dtype']); t = z['t']; m = t < n
    N = np.zeros((n, len(CLS))); np.add.at(N, (t[m], c[m]), 1); return N

def exec_start(R):
    tot = R.sum(1); m = pd.Series(tot).rolling(100).mean().values; m = np.where(np.isnan(m), 0, m); m[:3000] = 0
    return max(0, int(np.argmax(m > 150)) - 100)

def design(N):
    n = len(N); X = np.zeros((n, len(CLS) * K), np.float32)
    for c in range(len(CLS)):
        for k in range(K): X[k:, c * K + k] = N[:n - k, c]
    return X

def load(w):
    IC, BG = region_currents(w, split=True); R = IC + BG; N = counts(w, len(R)); s = exec_start(R); return IC[s:], N[s:], s, BG[s:]

if __name__ == '__main__':
    data = {w: load(w) for w in WL}
    Xf = np.vstack([design(data[w][1]) for w in FIT]); Yf = np.vstack([data[w][0] for w in FIT])
    mu_x = Xf.mean(0); mu_y = Yf.mean(0); Xc = Xf - mu_x
    A = Xc.T.dot(Xc) + LAM * len(Xf) / 1e4 * np.eye(Xf.shape[1])
    H = np.linalg.solve(A, Xc.T.dot(Yf - mu_y))                  # (C*K, regions)
    h = H.T.reshape(len(REG), len(CLS), K)
    print('kernel charge per uop (pC = mA x 3 ns summed over k), per region:')
    print('%-6s' % '' + ''.join('%16s' % r for r in REG) + '%10s' % 'total')
    for c, nm in enumerate(CLS):
        q = h[:, c, :].sum(1) * 3.0; print('%-6s' % nm + ''.join('%16.1f' % v for v in q) + '%10.1f' % q.sum())
    def score(ws):
        out = []
        for w in ws:
            R, N, s, BG = data[w]; P = design(N).dot(H) + mu_y - mu_x.dot(H); Bk = R - P
            yt = R.sum(1); yp = P.sum(1); fast = lambda x: x - pd.Series(x).rolling(64, center=True, min_periods=1).mean().values
            r2 = 1 - ((yt - yt.mean()) - (yp - yp.mean())).var() / yt.var()
            r2f = 1 - (fast(yt) - fast(yp)).var() / fast(yt).var()
            tot = (R + BG).sum(1); out.append((w, len(R), r2, r2f, yt.std(), tot.std(), np.cov(yt, tot)[0, 1] / tot.var()))
        return out
    print('\nissue-coupled current (exec phase): kernel R2 | R2 of fast part (<64 cyc) | std issue-coupled / std tile mA | share of tile variance')
    for tag, ws in (('FIT', FIT), ('HELD-OUT', TEST)):
        for w, n, r2, r2f, si, st, sh in score(ws): print('  %-9s %-10s n %7d  R2 %.2f  R2fast %.2f  std %.1f / %.1f  share %.2f' % (tag, w, n, r2, r2f, si, st, sh))
    np.savez(HERE + 'kernels.npz', h=h, mu_x=mu_x, mu_y=mu_y, CLS=np.array(CLS), REG=np.array(REG))
    for w in WL:                                                  # full-length background per bench for replay
        IC, BG = region_currents(w, split=True); N = counts(w, len(IC))
        np.save(B + w + '_bg.npy', (BG + IC - design(N).dot(H)).astype(np.float32))   # = tile - kernel part (incl. mean)
    print('saved kernels.npz and <wl>_bg.npy')
