"""Phase 1 gates for pdn5: (a) time-domain vs exact frequency-domain, (b) lumped/old-package equivalence vs pdn4,
(c) charge conservation, plus speed."""
import os, time, numpy as np
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pdn5, pdn4
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
I = np.load(B + 'qsort_tileI.npy')[105000:113192].astype(np.float32)          # 8192 cycles of qsort execution
A_x = 0.4 * pdn5.CELL_SHARE
# (a) exact check with the off-chip loop frozen
p = dict(A_x=A_x, fc_off=1e-6)
s = pdn5.Sim(p, Im_A=I.mean(0) * 1e-3); s.run(I); vmin, vavg, _ = s.run(I)             # second pass ~ periodic steady state
net = s.net; T = pdn5.T; n = len(I) * pdn5.SUB
dyn = np.maximum(I - pdn5.LEAK_T, 0); u = (pdn5.LEAK_T[None, None, :] + dyn[:, None, :] * pdn5.SHAPE[None, :, None]).reshape(n, T) * 1e-3
U = np.fft.rfft(u - u.mean(0), axis=0); f = np.fft.rfftfreq(n, pdn5.DT)
import scipy.sparse.linalg as spl
Bl = net['B'][:, 1:1 + T]; w = pdn5.CELL_SHARE; V = np.zeros(len(f), complex)
for k in range(1, len(f)):
    M = (2j * np.pi * f[k]) * net['E'] - net['A']; x = spl.spsolve(M.tocsc(), Bl.dot(U[k])); V[k] = x[net['didx']].dot(w)
v = np.fft.irfft(V, n).reshape(len(I), pdn5.SUB).mean(1) * 1e3; a = (vavg.dot(w) - vavg.dot(w).mean()) * 1e3
print('(a) exact check: td std %.2f fd std %.2f max|diff| %.3f mV corr %.5f' % (a.std(), v.std(), np.abs(a - v).max(), np.corrcoef(a, v)[0, 1]))
# (b) old-package equivalence: tiled network with the lumped 0.5 nH / 50 mOhm and a stiff grid vs pdn4 4-region
fq = np.logspace(5, 9, 120)
pt = dict(L_plane=0.5e-9 - 0.25e-9 / 40, R_plane=50e-3 - 20e-3 / 40, R_edge=0.005, A_x=0.4 * pdn5.CELL_SHARE)
z5 = pdn5.impedance(pdn5.build(pt), fq)
P4 = dict(pdn4.DEFAULT, A_x=[0.4 * s for s in pdn4.SHARE]); Pm, Qm, vidx, ridx, vo, dcf, nx, _ = pdn4.build(P4); dt = pdn4.DT
Ac = 2 / dt * np.linalg.solve(Pm + np.eye(nx), Pm - np.eye(nx)); Bc = (np.eye(nx) - Ac * dt / 2).dot(Qm) * 2 / dt
z4 = np.array([abs(np.linalg.solve(2j * np.pi * fk * np.eye(nx) - Ac, -Bc[:, 1:5].dot(pdn4.SHARE))[vidx].dot(pdn4.SHARE)) for fk in fq])
for fk in (1e6, 3.6e6, 2e7, 5e7, 7e7, 1e8, 3e8):
    i = np.argmin(abs(fq - fk)); print('(b) %7.1f MHz  tiled %6.1f mOhm  pdn4 %6.1f mOhm  ratio %.3f' % (fq[i] / 1e6, z5[i] * 1e3, z4[i] * 1e3, z5[i] / z4[i]))
print('(b) peak: tiled %.1f mOhm @ %.1f MHz, pdn4 %.1f mOhm @ %.1f MHz' % (z5.max() * 1e3, fq[z5.argmax()] / 1e6, z4.max() * 1e3, fq[z4.argmax()] / 1e6))
# speed: full qsort execution, open loop
J = np.load(B + 'qsort_tileI.npy')[45000:].astype(np.float32)
t0 = time.time(); s2 = pdn5.Sim(dict(A_x=A_x), Im_A=J[:3000].mean(0) * 1e-3); vmin2, vavg2, _ = s2.run(J); t1 = time.time()
d = (pdn5.VNOM - vmin2) * 1e3; worst = d.max(1)
print('speed: %d cycles x %d tiles in %.1f s (%.0f cycles/s); worst-tile droop p99.9 %.1f max %.1f mV; tile of max %d'
      % (len(J), pdn5.T, t1 - t0, len(J) / (t1 - t0), np.percentile(worst, 99.9), worst.max(), d.max(0).argmax()))
