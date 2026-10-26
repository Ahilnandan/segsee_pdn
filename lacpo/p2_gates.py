"""2.5 apparatus gate (resonant vs smooth synthetic load at equal mean) and 2.1 predictor accuracy."""
import os, numpy as np
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pdn5, closedloop as CLM, kernels as KM
# 2.5: 40 mA extra on the core column tiles, either as a 4-cycle on/off pattern (~83 MHz) or smooth
n = 6000; base = np.tile(pdn5.LEAK_T * 2.0, (n, 1)).astype(np.float32); core = np.array([t for t in range(64) if t % 8 in (3, 4)])
res = {}
for name, pat in (('resonant 4-cycle', (np.arange(n) // 2) % 2), ('smooth', np.full(n, 0.5))):
    I = base.copy(); I[:, core] += (80.0 / len(core)) * pat[:, None]
    s = pdn5.Sim(dict(A_x=0.3 * pdn5.CELL_SHARE), Im_A=I[:500].mean(0) * 1e-3); vmin, _, _ = s.run(I)
    d = (pdn5.VNOM - vmin[1000:]) * 1e3; res[name] = (np.percentile(d.max(1), 99.9), d.max())
    print('2.5 %-18s mean extra 40 mA  worst-tile droop p99.9 %.1f max %.1f mV' % (name, *res[name]))
# 2.1: issue-queue predictor vs actual issue-coupled current (total), qsort execution
cl = CLM.CL('qsort'); rp = cl.rp; tt = rp.t; m = rp.m; n = cl.n
pred = np.zeros(n + 16); act = np.zeros(n + 16)
for k in range(16):
    np.add.at(pred, tt + k, cl.hp[rp.cls, k]); np.add.at(act, tt + k, cl.A5[:, :, k].sum(1))
s0 = 45000; P, A = pred[s0:n], act[s0:n]
r2 = 1 - ((A - A.mean()) - (P - P.mean())).var() / A.var()
fast = lambda x: x - np.convolve(x, np.ones(32) / 32, 'same')
r2f = 1 - (fast(A) - fast(P)).var() / fast(A).var()
ev = fast(A) > np.percentile(fast(A), 99); pe = fast(P) > np.percentile(fast(P), 99)
print('2.1 predictor (issue-coupled total current): R2 %.2f, R2 of fast part %.2f, top-1%% event recall %.2f precision %.2f'
      % (r2, r2f, (ev & pe).sum() / ev.sum(), (ev & pe).sum() / pe.sum()))
tot_act = act[s0:n] + cl.BG[s0:n].sum(1)
print('    share of TOTAL tile current variance that is issue-coupled: %.2f' % (np.cov(A, tot_act)[0, 1] / tot_act.var()))
