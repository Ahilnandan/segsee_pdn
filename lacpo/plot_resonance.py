"""Explainer figure: package L + die C resonance (values from the pdn4 model: L 0.5 nH, R 50 mOhm, C 8.6 nF, ESR 35 mOhm)."""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
L, R, C, ESR = 0.5e-9, 50e-3, 8.6e-9, 35e-3
f0 = 1 / (2 * np.pi * np.sqrt(L * C))
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.labelcolor': INK, 'axes.facecolor': SURF, 'figure.facecolor': SURF})
fig = plt.figure(figsize=(12, 8.2)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.42, wspace=0.25)

# (a) schematic
ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis('off')
ax.set_title('(a) what the die sees: supply through the package, decap on the die', loc='left', fontsize=9.5)
def wire(x, y, **k): ax.plot(x, y, color=INK, lw=1.6, **k)
wire([0.6, 1.4], [4.5, 4.5]); ax.add_patch(Rectangle((0.2, 3.9), 0.8, 1.2, fc='#dfe9f7', ec=INK, lw=1.2))
ax.text(0.6, 4.5, '1.1 V\nsupply', ha='center', va='center', fontsize=7.5)
# inductor coils
xs = np.linspace(1.4, 3.6, 200); wire(xs, 4.5 + 0.18 * np.abs(np.sin((xs - 1.4) / 2.2 * 4 * np.pi)))
ax.text(2.4, 5.05, 'package L = 0.5 nH', ha='center', fontsize=8)
wire([3.6, 4.0], [4.5, 4.5]); zz = np.linspace(4.0, 5.2, 7); wire(np.linspace(4.0, 5.2, 13), 4.5 + 0.15 * np.array([0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 0]))
ax.text(4.75, 5.05, 'R = 50 mΩ', ha='center', fontsize=8)
wire([5.2, 8.6], [4.5, 4.5]); ax.plot(6.2, 4.5, 'o', color=INK, ms=5); ax.text(6.2, 4.85, 'die node', ha='center', fontsize=8, color=INK2)
# capacitor
wire([6.2, 6.2], [4.5, 2.9]); wire([5.75, 6.65], [2.9, 2.9]); wire([5.75, 6.65], [2.6, 2.6]); wire([6.2, 6.2], [2.6, 1.0])
ax.text(6.85, 2.75, 'die\ndecap\n8.6 nF', fontsize=8, va='center')
# load
wire([8.6, 8.6], [4.5, 3.4]); ax.add_patch(plt.Circle((8.6, 2.8), 0.55, fc=SURF, ec=INK, lw=1.4)); ax.annotate('', (8.6, 2.45), (8.6, 3.15), arrowprops=dict(arrowstyle='->', color=INK, lw=1.4))
ax.text(9.25, 2.8, 'core\nload i(t)', fontsize=8, va='center'); wire([8.6, 8.6], [2.25, 1.0])
wire([0.6, 0.6], [3.9, 1.0]); wire([0.6, 8.6], [1.0, 1.0]); ax.text(4.6, 0.55, 'ground', ha='center', fontsize=8, color=INK2)

# (b) impedance
ax = fig.add_subplot(gs[0, 1]); f = np.logspace(6.5, 9, 600); w = 2 * np.pi * f
ZL = R + 1j * w * L; ZC = ESR + 1 / (1j * w * C); Zp = ZL * ZC / (ZL + ZC)
ax.loglog(f, np.abs(ZL), color=C1, lw=1.6, ls='--', label='package path |R + jωL| (rises)')
ax.loglog(f, np.abs(ZC), color=C3, lw=1.6, ls='--', label='decap |1/jωC| (falls)')
ax.loglog(f, np.abs(Zp), color=C2, lw=2.4, label='what the load sees (both in parallel)')
ax.axvline(f0, color=INK2, lw=0.8, ls=':'); ax.text(f0 * 1.06, 0.014, 'resonance %.0f MHz' % (f0 / 1e6), fontsize=8, color=INK2)
ax.axhline(0.22, color=INK2, lw=0.8); ax.text(3.5e6, 0.24, 'target 0.22 Ω (55 mV / 0.25 A)', fontsize=7.5, color=INK2)
ax.set_ylim(0.01, 10); ax.set_xlabel('frequency of the current change (Hz)'); ax.set_ylabel('impedance (Ω)')
ax.set_title('(b) below f0 the supply wins, above f0 the decap wins, at f0 neither', loc='left', fontsize=9.5)
ax.legend(frameon=False, fontsize=7.5, loc='upper left', bbox_to_anchor=(0.18, 1.0)); ax.grid(True, which='both', color=GRID, lw=0.5)
for s in ('top', 'right'): ax.spines[s].set_visible(False)

# (c) time domain: same 100 mA load wiggle at three frequencies
ax = fig.add_subplot(gs[1, :]); dt = 0.02e-9; t = np.arange(0, 300e-9, dt)
def run(fl):
    i_load = 0.05 * (1 - np.cos(2 * np.pi * fl * t)) * (t > 20e-9)          # 0 -> 100 mA swing starting at 20 ns
    iL, vC, out = 0.0, 0.0, []
    for k in range(len(t)):                                                  # node voltage deviation v = vC + ESR*(iL - load)
        v = (vC * 1 + ESR * (iL - i_load[k]))
        diL = (-v - R * iL) / L; dvC = (iL - i_load[k]) / C
        iL += diL * dt; vC += dvC * dt; out.append(v)
    return np.array(out) * 1e3
for fl, c, lab in ((f0 / 5, C1, 'load wiggles at f0/5 (%.0f MHz): supply follows' % (f0 / 5e6)),
                   (f0, C2, 'load wiggles AT f0 (%.0f MHz): swing builds up' % (f0 / 1e6)),
                   (f0 * 4, C3, 'load wiggles at 4·f0 (%.0f MHz): decap absorbs it' % (4 * f0 / 1e6))):
    ax.plot(t * 1e9, run(fl), color=c, lw=1.6 if fl != f0 else 2.2, label=lab)
ax.axhline(-55, color=INK2, lw=0.8, ls='--'); ax.text(255, -52, '−55 mV budget', fontsize=8, color=INK2)
ax.set_xlabel('time (ns)'); ax.set_ylabel('die voltage change (mV)'); ax.grid(True, color=GRID, lw=0.5)
ax.set_title('(c) same 0–100 mA load wiggle at three speeds: only the one at resonance makes a big droop', loc='left', fontsize=9.5)
ax.legend(frameon=False, fontsize=8, loc='upper center', ncol=3); ax.set_ylim(-60, 55)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
out = '/home/ahilnandan/Documents/PMIC/pdn_work/docs/figs/resonance_explained.png'
fig.savefig(out, dpi=140, bbox_inches='tight'); print(out, 'f0 %.1f MHz, Zpeak %.3f ohm' % (f0 / 1e6, np.abs(Zp).max()))
