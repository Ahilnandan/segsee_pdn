"""Figure: (a) tile current and (b) die voltage around qsort's worst droop event for three decap designs,
(c) impedance seen by the die vs frequency (first / second droop resonances)."""
import os, sys, numpy as np
os.environ.setdefault('PDN3_LIB', 'libpdn3b.so')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn3, ga3
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
C = ['#2a78d6', '#eb6834', '#1baf7a']
DES = [('0.47 mm² decap, undamped', dict(A_x=[0.47 * s for s in pdn3.SHARE])),
       ('0.47 mm² decap, damped (0.3 Ω)', dict(A_x=[0.47 * s for s in pdn3.SHARE], R_d=[0.3] * 4)),
       ('0.70 mm² decap, undamped (reference)', dict(A_x=[0.70 * s for s in pdn3.SHARE]))]

z = np.load(ga3.B + 'qsort_rep_A.npz'); s0 = ga3.STARTS['qsort']; I = z['I'][s0:].astype(np.float32)
m = pdn3.simulate(I, DES[0][1], want_v=True)
k = int(np.argmin(m['vmin'].min(1)))                    # worst dip, undamped 0.47 mm2
a, b = k - 25, k + 40
seg = I[max(0, a - 3000):b]; off = a - max(0, a - 3000)
V = []
for _, p in DES:
    r = pdn3.simulate(seg, p, want_v=True, warm=500); V.append((r['vmin'].min(1)[off:] - pdn3.VNOM) * 1e3)
t = (np.arange(b - a) - (k - a)) * pdn3.T_CYC * 1e9         # ns relative to the worst dip

# impedance |Z(f)| seen by the die (mean die voltage per ampere of load, loads split by region share)
def zf(p, f):
    Pm, Qm, vidx, dc, nx, _ = pdn3.build(dict(pdn3.DEFAULT, **p)); dt = pdn3.DT
    Ac = 2 / dt * np.linalg.solve(Pm + np.eye(nx), Pm - np.eye(nx)); Bc = (np.eye(nx) - Ac * dt / 2).dot(Qm) * 2 / dt
    u = np.zeros(Qm.shape[1]); u[1:1 + pdn3.NR] = pdn3.SHARE
    out = []
    for fk in f:
        X = np.linalg.solve(2j * np.pi * fk * np.eye(nx) - Ac, Bc.dot(u)); out.append(abs(X[vidx[:pdn3.NR]].dot(pdn3.SHARE)))
    return np.array(out)
f = np.logspace(5, 9, 400); Z = [zf(p, f) for _, p in DES]

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'text.color': INK, 'axes.facecolor': SURF, 'figure.facecolor': SURF})
fig = plt.figure(figsize=(11, 7.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], hspace=0.38, wspace=0.22)
ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[1, 0], sharex=ax1); ax3 = fig.add_subplot(gs[:, 1])
for ax in (ax1, ax2, ax3):
    ax.grid(True, color=GRID, lw=0.6); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax1.step(t, I[a:b].sum(1), where='post', color=INK2, lw=1.4)
ax1.set_ylabel('tile current (mA)'); ax1.set_title('(a) qsort tile current, worst droop event (3 ns cycles)', loc='left', fontsize=9.5)
for (lab, _), v, c in zip(DES, V, C): ax2.plot(t, v, color=c, lw=2, label=lab)
ax2.axhline(-55, color=INK2, lw=1, ls='--'); ax2.text(t[0] + 1, -52, '−55 mV (5 %) budget', color=INK2, fontsize=8)
ax2.set_ylabel('die voltage − 1.1 V (mV, cycle minimum)'); ax2.set_xlabel('time relative to worst dip (ns)')
ax2.set_title('(b) die voltage: dip, then ringing (~4–5 cycles ≈ 77 MHz)', loc='left', fontsize=9.5)
ax2.legend(frameon=False, fontsize=8, loc='lower right')
for (lab, _), zz, c in zip(DES, Z, C): ax3.loglog(f, zz * 1e3, color=c, lw=2, label=lab)
i1 = np.argmax(Z[0] * (f > 2e7)); i2 = np.argmax(Z[0] * (f < 5e6) * (f > 1e5))
ax3.annotate('1st droop\npackage L × die C\n%.0f MHz' % (f[i1] / 1e6), (f[i1], Z[0][i1] * 1e3), xytext=(f[i1] * 0.08, Z[0][i1] * 1e3 * 0.85),
             fontsize=8, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
ax3.annotate('2nd droop\nboard L × package + bulk C\n%.1f MHz' % (f[i2] / 1e6), (f[i2], Z[0][i2] * 1e3), xytext=(f[i2] * 0.4, Z[0][i2] * 1e3 * 1.8),
             fontsize=8, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
ax3.set_xlabel('frequency (Hz)'); ax3.set_ylabel('impedance seen by the die (mΩ)')
ax3.set_title('(c) impedance the die sees vs frequency', loc='left', fontsize=9.5); ax3.legend(frameon=False, fontsize=8, loc='lower right')
out = '/home/ahilnandan/Documents/PMIC/pdn_work/docs/figs/droop_example.png'
fig.savefig(out, dpi=140, bbox_inches='tight'); print(out, 'worst dip cycle', k, 'droops', [round(float(v.min()), 1) for v in V])
