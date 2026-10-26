"""Time-domain view of the per-cycle core current (policy A), full length + zoomed windows, with phase boundaries."""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
TILE = ['fetch', 'decode', 'ren/ROB/RF', 'CSR', 'ALU0', 'ALU1', 'ALU2', 'ALU3', 'LSU0', 'LSU1', 'D$a', 'D$b', 'FPU0', 'FPU1', 'L2if', 'misc']
for wl in ('qsort', 'dhrystone', 'mm'):
    It = np.load(f'It_{wl}_A.npy').astype(float); I = It.sum(1); n = len(I); t_us = np.arange(n) / 1e3
    # envelope per 200-cycle block
    B = 200; m = n // B; blk = I[:m * B].reshape(m, B); tb = (np.arange(m) * B + B / 2) / 1e3
    fig = plt.figure(figsize=(11, 8.5)); gs = fig.add_gridspec(4, 2, height_ratios=[1.3, 1, 1, 1.1])
    ax = fig.add_subplot(gs[0, :])
    ax.fill_between(tb, blk.min(1), blk.max(1), color='0.75', lw=0, label='min-max per 200 cycles')
    ax.plot(tb, blk.mean(1), lw=0.8, color='C0', label='mean per 200 cycles')
    ax.set_xlim(0, t_us[-1]); ax.set_ylabel('core current (mA)'); ax.set_xlabel('time (µs, 1 cycle = 1 ns assumed)')
    ax.set_title(f'{wl}: per-cycle core current, full run ({n} cycles)', fontsize=10); ax.legend(fontsize=7, loc='upper right')
    # zoom windows: start, middle of main phase, near end
    picks = [('start', int(0.02 * n)), ('middle', int(0.45 * n)), ('end', int(0.97 * n))]
    for c, (lab, k0) in enumerate(picks):
        k1 = min(n, k0 + 300); x = np.arange(k0, k1)
        a = fig.add_subplot(gs[1 + c // 2, c % 2]) if c < 2 else fig.add_subplot(gs[2, 0])
        a.step(x, I[k0:k1], where='post', lw=0.8); a.set_title(f'zoom ({lab}): cycles {k0}-{k1}', fontsize=8)
        a.set_xlabel('cycle'); a.set_ylabel('mA'); ax.axvspan(k0 / 1e3, k1 / 1e3, color='C3', alpha=0.25)
    a = fig.add_subplot(gs[2, 1])
    k0 = int(0.45 * n); k1 = k0 + 120
    for j, tl in enumerate((4, 5, 8, 9, 12)):
        a.step(np.arange(k0, k1), It[k0:k1, tl] + 60 * j, where='post', lw=0.7, label=TILE[tl])
    a.set_title(f'per-unit current (offset 60 mA each), cycles {k0}-{k1}', fontsize=8); a.legend(fontsize=6, ncol=5, loc='upper left'); a.set_xlabel('cycle')
    # bottom: 50-cycle autocorrelation of the middle window -> loop period
    a = fig.add_subplot(gs[3, :])
    seg = I[int(0.3 * n):int(0.6 * n)]; seg = seg - seg.mean()
    ac = np.correlate(seg[:20000], seg[:20000], 'full')[19999:19999 + 61]; ac /= ac[0]
    a.stem(np.arange(61), ac, basefmt=' '); a.set_xlabel('lag (cycles)'); a.set_ylabel('autocorr')
    a.set_title('autocorrelation of the current in the main phase (peaks = repeating loop period)', fontsize=8)
    plt.tight_layout(); plt.savefig(f'docs/figs/time_{wl}.png', dpi=130); plt.close()
    print(wl, 'saved')
