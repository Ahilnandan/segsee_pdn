"""Plain current plots: raw per-cycle core current (every cycle, full run), and a stacked per-unit breakdown."""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
TILE = ['fetch', 'decode', 'ren/ROB/RF', 'CSR', 'ALU0', 'ALU1', 'ALU2', 'ALU3', 'LSU0', 'LSU1', 'D$a', 'D$b', 'FPU0', 'FPU1', 'L2if', 'misc']
GROUPS = [('ALU0', [4]), ('ALU1', [5]), ('ALU2+3', [6, 7]), ('LSU', [8, 9]), ('D-cache', [10, 11]), ('FPU', [12, 13]), ('rename/ROB/RF', [2]), ('CSR', [3]), ('static floor (rest)', [0, 1, 14, 15])]
KERNEL = {'qsort': (40, 235), 'dhrystone': (2, 67), 'mm': (7, 77)}     # µs, read off the time-domain plots
WL = ('qsort', 'dhrystone', 'mm')
# 1) raw per-cycle current, full run
fig, ax = plt.subplots(3, 1, figsize=(12, 8))
for a, wl in zip(ax, WL):
    I = np.load(f'It_{wl}_A.npy').astype(float).sum(1); t = np.arange(len(I)) / 1e3
    a.plot(t, I, lw=0.15, color='C0', rasterized=True)
    k0, k1 = KERNEL[wl]; a.axvspan(k0, k1, color='C2', alpha=0.08, label='kernel phase')
    a.set_xlim(0, t[-1]); a.set_ylabel('mA'); a.set_title(f'{wl}: core current, every cycle ({len(I)} cycles); mean {I.mean():.0f} mA, '
                                                          f'kernel mean {I[int(k0*1e3):int(k1*1e3)].mean():.0f} mA', fontsize=9)
    a.legend(fontsize=7, loc='upper right')
ax[-1].set_xlabel('time (µs; 1 cycle = 1 ns assumed)')
plt.tight_layout(); plt.savefig('docs/figs/current_full.png', dpi=130); plt.close()
# 2) stacked per-unit current, 400 cycles in the middle of each kernel phase
fig, ax = plt.subplots(3, 1, figsize=(12, 9))
for a, wl in zip(ax, WL):
    It = np.load(f'It_{wl}_A.npy').astype(float); k0, k1 = KERNEL[wl]; c = int((k0 + k1) / 2 * 1e3); x = np.arange(c, c + 400)
    layers = [It[c:c + 400][:, idx].sum(1) for _, idx in GROUPS]
    a.stackplot(x, *layers, labels=[g for g, _ in GROUPS], step='post', lw=0)
    a.step(x, It[c:c + 400].sum(1), where='post', color='k', lw=0.6)
    a.set_xlim(x[0], x[-1]); a.set_ylabel('mA'); a.set_title(f'{wl}: current by unit, 400 cycles from the middle of the kernel phase', fontsize=9)
ax[0].legend(fontsize=7, ncol=9, loc='upper left', bbox_to_anchor=(0, 1.32)); ax[-1].set_xlabel('cycle')
plt.tight_layout(); plt.savefig('docs/figs/current_units.png', dpi=130); plt.close()
print('saved')
