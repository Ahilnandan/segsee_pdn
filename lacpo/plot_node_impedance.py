"""Self-impedance seen at each point of the pdn4 passive reference network (0.9 mm^2 decap, no damping).
Unit AC current injected at a node, voltage read at the same node; off-chip buck = current-mode stage (its loop
only matters below ~100 kHz, not included)."""
import numpy as np, pdn4
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#6250d6']
p = dict(pdn4.DEFAULT, A_x=[0.9 * s for s in pdn4.SHARE])
Pm, Qm, vidx, ridx, vo, dc, nx, _ = pdn4.build(p); dt = pdn4.DT
Ac = 2 / dt * np.linalg.solve(Pm + np.eye(nx), Pm - np.eye(nx))
NN_names = ['vo', 'voc', 'pk', 'pkc', 'pk2', 'pk2c'] + ['d%d' % r for r in range(4)]
idx = {n: i for i, n in enumerate(NN_names)}
# E matrix diagonal for node rows: recover node capacitance by injecting through the B column of a load (scaled)
Bc = (np.eye(nx) - Ac * dt / 2).dot(Qm) * 2 / dt
inj = {'board (after off-chip buck)': Bc[:, 0], 'die region R0 frontend+I$': -Bc[:, 1], 'die region R1 int core': -Bc[:, 2],
       'die region R2 LSU+D$': -Bc[:, 3], 'die region R3 FP': -Bc[:, 4]}
read = {'board (after off-chip buck)': vo, 'die region R0 frontend+I$': vidx[0], 'die region R1 int core': vidx[1],
        'die region R2 LSU+D$': vidx[2], 'die region R3 FP': vidx[3]}
f = np.logspace(5, 9, 500); out = {}
for k in inj:
    z = []
    for fk in f:
        X = np.linalg.solve(2j * np.pi * fk * np.eye(nx) - Ac, inj[k]); z.append(abs(X[read[k]]))
    out[k] = np.array(z)
# package node: inject via unit vector scaled by node capacitance row (E^-1): use the bump-path difference method
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK, 'axes.facecolor': SURF, 'figure.facecolor': SURF})
fig, ax = plt.subplots(figsize=(9, 5.2))
for (k, z), c in zip(out.items(), COL):
    ax.loglog(f, z * 1e3, color=c, lw=2 if 'board' in k else 1.6, label=k)
    i = np.argmax(z); print('%-30s  |Z| @1MHz %7.1f mOhm  @10MHz %7.1f  peak %7.1f mOhm at %6.1f MHz  @300MHz %6.1f' % (
        k, 1e3 * z[np.argmin(abs(f - 1e6))], 1e3 * z[np.argmin(abs(f - 1e7))], 1e3 * z.max(), f[i] / 1e6, 1e3 * z[np.argmin(abs(f - 3e8))]))
ax.axhline(220, color=INK2, lw=0.8, ls='--'); ax.text(1.2e5, 240, 'die target 0.22 Ω', fontsize=8, color=INK2)
ax.set_xlabel('frequency (Hz)'); ax.set_ylabel('self-impedance (mΩ)'); ax.grid(True, which='both', color=GRID, lw=0.5)
ax.set_title('Impedance seen at each point (passive reference, 0.9 mm² decap)', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=8, loc='upper left')
for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.savefig('/home/ahilnandan/Documents/PMIC/pdn_work/docs/figs/node_impedance.png', dpi=140, bbox_inches='tight')
