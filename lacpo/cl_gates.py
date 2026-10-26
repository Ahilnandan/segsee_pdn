"""Gates 1.4d/e: closed-loop baseline must reproduce the open-loop replay grants exactly, and its droop must equal
pdn5 driven open-loop with the same per-tile currents. Plus speed."""
import os, time, numpy as np
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import closedloop as CLM, replay3, pdn5
w = 'vvadd'; cl = CLM.CL(w); D = dict(A_x=0.4 * pdn5.CELL_SHARE)
t0 = time.time(); r = cl.run(D, 'A', record_t=True); t1 = time.time()
ol = cl.rp.run('A', record=False)
same = (r['tnew'] == ol['tnew'][r['sel']]).mean()
print('1.4d grants identical to open-loop replay: %.4f   (%d cycles in %.1fs, %.0f cyc/s)' % (same, cl.n, t1 - t0, cl.n / (t1 - t0)))
# open-loop pdn5 on the reconstructed baseline tile currents (replay with zero slip == original order)
It = cl.BG.copy()
z = np.load(CLM.B + w + '_attr5.npz'); A = z['A']; tt = z['t']
for k in range(CLM.K):
    ok = tt + k < len(It); np.add.at(It, tt[ok] + k, A[ok, :, k].dot(cl.WB))
s = pdn5.Sim(D, Im_A=cl.BG[:3000].mean(0) * 1e-3); vmin, _, _ = s.run(It.astype(np.float32))
d = (pdn5.VNOM - vmin[600:]) * 1e3; wo = d.max(1)
print('1.4e open-loop pdn5 on baseline currents: p99.9 %.2f max %.2f | closed-loop A: p99.9 %.2f max %.2f' % (np.percentile(wo, 99.9), wo.max(), r['p999'], r['max']))
