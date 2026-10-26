"""How a load step is made of frequencies, and how each frequency passes through the PDN impedance.
Square-wave load (0 <-> +150 mA, 400 ns period) = sum of odd sine harmonics with amplitude (4/pi)(dI/2)/n."""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'
L, R, C, ESR = 0.5e-9, 50e-3, 8.6e-9, 35e-3
T = 400e-9; f1 = 1 / T; dI = 0.15; t = np.linspace(0, T, 20000, endpoint=False)
def Z(f):
    w = 2 * np.pi * f; ZL = R + 1j * w * L; ZC = ESR + 1 / (1j * w * C); return ZL * ZC / (ZL + ZC)
def build(nmax, through_Z=False):
    s = np.full_like(t, dI / 2 if not through_Z else 0.0)
    if through_Z: s += -Z(1e-3).real * dI / 2 * 0
    for n in range(1, nmax + 1, 2):
        a = (4 / np.pi) * (dI / 2) / n; f = n * f1
        if through_Z:
            z = Z(f); s += -a * abs(z) * np.sin(2 * np.pi * f * t + np.angle(z))
        else:
            s += a * np.sin(2 * np.pi * f * t)
    return s
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK, 'axes.facecolor': SURF, 'figure.facecolor': SURF})
fig = plt.figure(figsize=(11, 9)); gs = fig.add_gridspec(3, 1, hspace=0.45)
ax = fig.add_subplot(gs[0]); tn = t * 1e9
for n, c in ((1, C1), (5, C3), (31, C4), (401, C2)):
    ax.plot(tn, build(n) * 1e3, color=c, lw=1.4 if n < 401 else 2, label='sum of sines up to %d × %.1f MHz' % (n, f1 / 1e6))
ax.set_xlim(0, T * 1e9); ax.set_ylabel('load current (mA)'); ax.legend(frameon=False, fontsize=8, loc='upper right', ncol=2)
ax.set_title('(a) a step (here repeated, 0 ↔ 150 mA) = slow sines for the level + ever faster sines for the sharp edges', loc='left', fontsize=9.5)
ax = fig.add_subplot(gs[1]); n = np.arange(1, 4001, 2); f = n * f1; amp = (4 / np.pi) * (dI / 2) / n
ax.loglog(f, amp * 1e3, 'o', ms=2, color=C1, label='step content: amplitude of each sine (mA) ∝ 1/f')
ax.loglog(f, np.abs(Z(f)) * 1e3, color=C2, lw=2, label='PDN impedance |Z(f)| (mΩ)')
ax.loglog(f, amp * np.abs(Z(f)) * 1e3, color=C3, lw=2, label='voltage each sine produces (mV) = amplitude × |Z|')
ax.set_xlabel('frequency (Hz)'); ax.grid(True, which='both', color=GRID, lw=0.4); ax.legend(frameon=False, fontsize=8, loc='lower left')
ax.set_title('(b) multiply each sine by the impedance at its frequency: the ones near the 77 MHz peak get amplified', loc='left', fontsize=9.5)
ax = fig.add_subplot(gs[2]); v = build(4001, through_Z=True) * 1e3
ax.plot(tn, v - v[:50].mean() * 0, color=C3, lw=1.8, label='die voltage = sum of every sine × its impedance')
ax.set_xlim(0, T * 1e9); ax.set_xlabel('time (ns)'); ax.set_ylabel('die voltage change (mV)'); ax.grid(True, color=GRID, lw=0.4)
ax.set_title('(c) adding them back up gives the dip + ringing at each edge (fast sines) and the level shift (slow sines × R)', loc='left', fontsize=9.5)
for a in fig.axes: [a.spines[s].set_visible(False) for s in ('top', 'right')]
fig.savefig('/home/ahilnandan/Documents/PMIC/pdn_work/docs/figs/step_fourier.png', dpi=140, bbox_inches='tight')
print('dip after rising edge %.1f mV' % (v[:5000].min() - v[4000:4900].mean()))
