"""Independent check of the C++ trapezoidal PDN kernel: frequency-domain solve of the same linear network
(x = (jwE - A)^-1 B u per FFT bin, same intra-cycle load shape, VR loop off) vs time-domain kernel."""
import numpy as np, pdn3
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
I = np.load(B + 'vvadd_rep_A.npz')['I'][5000:13192].astype(float)          # 8192 cycles
p = dict(pdn3.DEFAULT, A_x=[0.25 * s for s in pdn3.SHARE], fc_vr=1e-3)
# time domain (VR loop effectively off); periodic steady state approximated by repeating the window twice
m = pdn3.simulate(np.vstack([I, I]), p, warm=10, want_v=True); vt = m['vavg'][len(I):]
# frequency domain on the same network
import scipy.linalg as sl
p2 = dict(p); Pm, Qm, vidx, dc, nx, _ = pdn3.build(p2)
# rebuild E, A, B (copy of build internals)
def mats(p):
    import types
    src = open(pdn3.__file__).read()
    return None
NR = pdn3.NR; sub = pdn3.SUB; n = len(I) * sub
leak = pdn3.LEAK_MA; dyn = np.maximum(I - leak, 0)
u = (leak[None, None, :] + dyn[:, None, :] * pdn3.SHAPE[None, :, None]).reshape(n, NR) * 1e-3
# continuous-time matrices from the discrete ones: trapezoid is the bilinear map, invert it exactly
Pm = np.asarray(Pm); Qm = np.asarray(Qm); dt = pdn3.DT
Ac = 2 / dt * np.linalg.solve(Pm + np.eye(nx), Pm - np.eye(nx)); Bc = np.linalg.solve(Pm + np.eye(nx), Qm) * 2 / dt * 1.0
Bc = (np.eye(nx) - Ac * dt / 2).dot(Qm) * 2 / dt
U = np.fft.rfft(u - u.mean(0), axis=0); f = np.fft.rfftfreq(n, dt); V = np.zeros((len(f), NR), complex)
for k, fk in enumerate(f):
    if k == 0: continue
    X = np.linalg.solve(2j * np.pi * fk * np.eye(nx) - Ac, Bc[:, 1:1 + NR].dot(U[k]))
    V[k] = X[vidx[:NR]]
v = np.fft.irfft(V, n, axis=0).reshape(len(I), sub, NR).mean(1)
a = (vt - vt.mean(0)) * 1e3; b = v * 1e3
print('cycle-avg die voltage, mV: time-domain std %.2f  freq-domain std %.2f  max |diff| %.2f  corr %.4f' % (a.std(), b.std(), np.abs(a - b).max(), np.corrcoef(a.ravel(), b.ravel())[0, 1]))
