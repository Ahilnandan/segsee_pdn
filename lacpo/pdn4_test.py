"""pdn4 checks: (1) linear network vs exact frequency-domain solve (off-chip loop frozen),
(2) passive design vs pdn3 on the same window, (3) each new device sane on a real window."""
import os, numpy as np
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pdn4, ga3
B = ga3.B; z = np.load(B + 'qsort_rep_A.npz'); s0 = ga3.STARTS['qsort']; I = z['I'][s0 + 60000:s0 + 68192].astype(np.float32)
# (1) exact check with constant off-chip current (fc_off ~ 0, slew 0)
p = dict(pdn4.DEFAULT, A_x=[0.4 * s for s in pdn4.SHARE], fc_off=1e-6)
m = pdn4.simulate(np.vstack([I, I]), p, warm=10, want_v=True); vt = m['vavg'][len(I):]
Pm, Qm, vidx, ridx, vo, dc, nx, _ = pdn4.build(p); dt = pdn4.DT
Ac = 2 / dt * np.linalg.solve(Pm + np.eye(nx), Pm - np.eye(nx)); Bc = (np.eye(nx) - Ac * dt / 2).dot(Qm) * 2 / dt
NR = 4; sub = pdn4.SUB; n = len(I) * sub; leak = pdn4.LEAK_MA; dyn = np.maximum(I - leak, 0)
u = (leak[None, None, :] + dyn[:, None, :] * pdn4.SHAPE[None, :, None]).reshape(n, NR) * 1e-3
U = np.fft.rfft(u - u.mean(0), axis=0); f = np.fft.rfftfreq(n, dt); V = np.zeros((len(f), NR), complex)
for k, fk in enumerate(f):
    if k == 0: continue
    X = np.linalg.solve(2j * np.pi * fk * np.eye(nx) - Ac, Bc[:, 1:1 + NR].dot(U[k])); V[k] = X[vidx]
v = np.fft.irfft(V, n, axis=0).reshape(len(I), sub, NR).mean(1); a = (vt - vt.mean(0)) * 1e3; b = v * 1e3
print('(1) exact check: td std %.2f fd std %.2f max|diff| %.3f mV corr %.5f' % (a.std(), b.std(), np.abs(a - b).max(), np.corrcoef(a.ravel(), b.ravel())[0, 1]))
# (2) vs pdn3 (ideal source) on a full exec slice
os.environ['PDN3_LIB'] = 'libpdn3b.so'; import pdn3
J = z['I'][s0:s0 + 120000].astype(np.float32)
for Ax in (0.47, 0.7):
    a3 = pdn3.simulate(J, dict(A_x=[Ax * s for s in pdn3.SHARE])); a4 = pdn4.simulate(J, dict(A_x=[Ax * s for s in pdn4.SHARE]))
    print('(2) %.2f mm2  pdn3 p99.9 %.1f max %.1f | pdn4 p99.9 %.1f max %.1f mV, mean v %.4f, off-chip cond loss %.2f%%' % (Ax, a3['p999'], a3['max'], a4['p999'], a4['max'], a4['mean_v'], 100 * a4['loss']['off']))
# (3) devices
P = dict(A_x=[0.3 * s for s in pdn4.SHARE])
def show(tag, spec=None, p=P, **kw):
    m = pdn4.simulate(J, dict(pdn4.DEFAULT, **p), spec, **kw); ar = pdn4.area(dict(pdn4.DEFAULT, **p), spec)
    print('%-62s p99.9 %5.1f max %5.1f over %5.1f | loss %s | ripple %s | area %.3f | ion %.1f mA' % (tag, m['p999'], m['max'], m['over_max'],
          {k: round(100 * v, 2) for k, v in m['loss'].items()}, [round(x, 1) for x in m['ripple_mV']], ar['total'], 1e3 * m['ion_mean']))
show('passive 0.3 mm2')
reg = lambda role, **k: [dict(reg=dict(dict(role=role, L=100e-9, on_die=False, Vin=1.8, Ron=1.0, phases=4, fsw=100e6, fc=20e6, Ipk=0.1 * s, th=0.03), **k)) for s in pdn4.SHARE]
show('A: off-chip V + on-die injector (pkg 20nH, 1.3V, Ron .5)', [dict(reg=dict(role=1, L=20e-9, on_die=False, Vin=1.3, Ron=0.5, Ipk=0.2 * s, th=0.04)) for s in pdn4.SHARE])
show('B: off-chip I-helper + on-die regulator 4ph 100nH pkg 100MHz fc20M', reg(2), role_off=np.ones(len(J), np.uint8))
show('B + predictive boost 15 mV', reg(2, boost=0.015), role_off=np.ones(len(J), np.uint8),
     pred=pdn4.predict_triggers(z['Ihat'][s0:s0 + 120000], 6, 3))
show('both regulate (no role split): off-chip V + on-die regulator', reg(2))
ldo = lambda series, **k: [dict(ldo=dict(dict(series=series, V=1.2, R_pass=0.3, gm=5.0, tau=1e-9, off=0.02), **k)) for s in pdn4.SHARE]
show('E: parallel LDO clamp, real 1.2 V rail (own pkg path, 0.47uF)', ldo(0), p=dict(A_x=[0.3 * s for s in pdn4.SHARE], ldo_rail=True))
show('E + rail decap 0.05 mm2 + predictive boost 15 mV', ldo(0, boost=0.015), p=dict(A_x=[0.3 * s for s in pdn4.SHARE], ldo_rail=True, A_rail=[0.05 * s for s in pdn4.SHARE]), pred=pdn4.predict_triggers(z['Ihat'][s0:s0 + 120000], 6, 3))
show('D: series LDO, rail +80 mV, rail decap 0.1 mm2', ldo(1, off=0.0, gm=10.0), p=dict(A_x=[0.3 * s for s in pdn4.SHARE], series_ldo=True, rail_off=0.08, A_rail=[0.1 * s for s in pdn4.SHARE]))
