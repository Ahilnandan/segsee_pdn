"""Time-frequency (STFT) analysis of the core current, total and per floorplan tile, before any regulator design.
Sampling: 1 per cycle (1 GHz). Bands: B0 <10 MHz, B1 10-50, B2 50-200 (die-C / package-L resonance, peak ~92 MHz), B3 200-500.
Inputs: It_<wl>_A.npy (cycles x 16 tiles, mA) from replay2.run(tiles=True), policy A (unmodified BOOM issue)."""
import json, numpy as np, os
from scipy import signal
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import pbuck_freq as M

FS = 1e9; FIG = 'docs/figs'
BANDS = [('B0 <10 MHz', 0.1e6, 10e6), ('B1 10-50 MHz', 10e6, 50e6), ('B2 50-200 MHz', 50e6, 200e6), ('B3 200-500 MHz', 200e6, 500e6)]
TILE = ['fetch', 'decode', 'ren/ROB/RF', 'CSR', 'ALU0', 'ALU1', 'ALU2', 'ALU3', 'LSU0', 'LSU1', 'D$a', 'D$b', 'FPU0', 'FPU1', 'L2if', 'misc']
def Zmag(ff):
    s0 = M.s; M.s = 2j * np.pi * np.maximum(ff, 1.0); Z = np.abs(M.Z_total()); M.s = s0; return Z

def analyse(wl, nperseg=256, hop=64):
    It = np.load(f'It_{wl}_A.npy').astype(float); I = It.sum(1); T = len(I)
    f, t, S = signal.stft(I - I.mean(), fs=FS, window='hann', nperseg=nperseg, noverlap=nperseg - hop, boundary=None, padded=False)
    P = np.abs(S) ** 2
    out = dict(wl=wl, cycles=T, mean_mA=float(I.mean()), std_mA=float(I.std()))
    # 1. band fractions of current power and of droop power (|Z|^2-weighted)
    Z = Zmag(f); tot = P[f > 0.1e6].sum(); totd = (P * (Z ** 2)[:, None])[f > 0.1e6].sum()
    out['band_frac_current'] = {b: float(P[(f >= lo) & (f < hi)].sum() / tot) for b, lo, hi in BANDS}
    out['band_frac_droop'] = {b: float((P * (Z ** 2)[:, None])[(f >= lo) & (f < hi)].sum() / totd) for b, lo, hi in BANDS}
    # 2. burstiness of the resonance band over time
    m2 = (f >= 50e6) & (f < 200e6); P2 = P[m2].sum(0)
    med = np.median(P2); thr = 4 * med; above = P2 > thr
    runs = []; r = 0
    for a in above:
        if a: r += 1
        elif r: runs.append(r); r = 0
    if r: runs.append(r)
    runs = np.array(runs) if runs else np.array([0])
    out['B2_burst'] = dict(p99_over_median=float(np.percentile(P2, 99) / med), frac_windows_above_4x_median=float(above.mean()),
                           burst_len_median_ns=float(np.median(runs) * hop), burst_len_p99_ns=float(np.percentile(runs, 99) * hop),
                           n_bursts=int(len(runs)), corr_with_mean_current=float(np.corrcoef(P2, np.abs(S[f < 2e6]).sum(0))[0, 1]))
    # 3. per tile: share of resonance-band energy and of droop energy (each tile alone through the same |Z|)
    share, dshare, Pt2 = [], [], []
    for k in range(16):
        fk, tk, Sk = signal.stft(It[:, k] - It[:, k].mean(), fs=FS, window='hann', nperseg=nperseg, noverlap=nperseg - hop, boundary=None, padded=False)
        Pk = np.abs(Sk) ** 2; share.append(Pk[m2].sum()); dshare.append((Pk * (Z ** 2)[:, None])[f > 0.1e6].sum()); Pt2.append(Pk[m2].sum(0))
    share = np.array(share); Pt2 = np.array(Pt2)
    out['tile_B2_share'] = {TILE[k]: float(share[k] / share.sum()) for k in range(16)}
    out['tile_droop_share'] = {TILE[k]: float(dshare[k] / sum(dshare)) for k in range(16)}
    # 4. coherence: energy of the sum vs sum of tile energies in each band (1 = independent tiles, N = fully coherent)
    G = {}
    for b, lo, hi in BANDS:
        mb = (f >= lo) & (f < hi)
        sum_tiles = 0.0
        for k in range(16):
            fk, tk, Sk = signal.stft(It[:, k] - It[:, k].mean(), fs=FS, window='hann', nperseg=nperseg, noverlap=nperseg - hop, boundary=None, padded=False)
            sum_tiles += (np.abs(Sk[mb]) ** 2).sum()
        G[b] = float(P[mb].sum() / sum_tiles)
    out['coherence_gain'] = G
    # pairwise magnitude-squared coherence in B2 for the busy tiles
    busy = [2, 4, 5, 6, 8, 9, 10, 12]
    C = np.zeros((len(busy), len(busy)))
    for i, a in enumerate(busy):
        for j, b in enumerate(busy):
            if j <= i: continue
            fc, coh = signal.coherence(It[:, a], It[:, b], fs=FS, nperseg=nperseg)
            C[i, j] = C[j, i] = coh[(fc >= 50e6) & (fc < 200e6)].mean()
    out['B2_pair_coherence'] = {f'{TILE[a]}-{TILE[b]}': float(C[i, j]) for i, a in enumerate(busy) for j, b in enumerate(busy) if j > i}
    # ---- figures
    fig, ax = plt.subplots(2, 1, figsize=(8, 5.6), gridspec_kw=dict(height_ratios=[2.2, 1]))
    span = slice(0, min(P.shape[1], 1500))
    im = ax[0].pcolormesh(t[span] * 1e6, f / 1e6, 10 * np.log10(P[:, span] + 1e-3), shading='auto', cmap='magma')
    ax[0].axhspan(50, 200, color='w', alpha=0.08); ax[0].axhline(92, color='c', lw=0.6, ls='--')
    ax[0].set_ylabel('MHz'); ax[0].set_title(f'{wl}: core current STFT (256-cycle Hann, hop 64), dashed = supply resonance 92 MHz', fontsize=9)
    plt.colorbar(im, ax=ax[0], label='dB (mA²)')
    ax[1].plot(t[span] * 1e6, 10 * np.log10(P2[span] + 1e-3), lw=0.6); ax[1].axhline(10 * np.log10(thr), color='r', lw=0.6, ls='--')
    ax[1].set_ylabel('50-200 MHz\npower (dB)'); ax[1].set_xlabel('time (µs)')
    plt.tight_layout(); plt.savefig(f'{FIG}/stft_{wl}.png', dpi=140); plt.close()
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.4))
    for a, (title, v) in zip(ax, [('share of 50-200 MHz current energy', share / share.sum()), ('share of droop energy (|Z|-weighted)', np.array(dshare) / sum(dshare))]):
        im = a.imshow(v.reshape(4, 4), cmap='viridis'); a.set_title(title, fontsize=8)
        for k in range(16): a.text(k % 4, k // 4, f'{TILE[k]}\n{100*v[k]:.0f}%', ha='center', va='center', fontsize=6, color='w' if v[k] < v.max() * 0.6 else 'k')
        a.set_xticks([]); a.set_yticks([])
    plt.suptitle(f'{wl}: where the resonance-band current comes from', fontsize=9); plt.tight_layout(); plt.savefig(f'{FIG}/tiles_{wl}.png', dpi=140); plt.close()
    plt.figure(figsize=(8, 3.6))
    sp = slice(0, min(Pt2.shape[1], 600))
    plt.pcolormesh(t[sp] * 1e6, np.arange(16), 10 * np.log10(Pt2[:, sp] + 1e-3), shading='auto', cmap='magma')
    plt.yticks(range(16), TILE, fontsize=6); plt.xlabel('time (µs)'); plt.colorbar(label='50-200 MHz power (dB)')
    plt.title(f'{wl}: resonance-band power per tile over time', fontsize=9); plt.tight_layout(); plt.savefig(f'{FIG}/tile_time_{wl}.png', dpi=140); plt.close()
    return out

if __name__ == '__main__':
    R = [analyse(wl) for wl in ('qsort', 'dhrystone', 'mm')]
    json.dump(R, open('stft_analysis.json', 'w'), indent=1)
    for r in R:
        print(f"\n== {r['wl']}: {r['cycles']} cycles, mean {r['mean_mA']:.0f} mA, std {r['std_mA']:.0f} mA")
        print('  current energy by band :', {k: f'{100*v:.0f}%' for k, v in r['band_frac_current'].items()})
        print('  droop energy by band   :', {k: f'{100*v:.0f}%' for k, v in r['band_frac_droop'].items()})
        b = r['B2_burst']
        print(f"  50-200 MHz burstiness: p99/median {b['p99_over_median']:.1f}x, {100*b['frac_windows_above_4x_median']:.1f}% of time above 4x median, "
              f"bursts median {b['burst_len_median_ns']:.0f} ns / p99 {b['burst_len_p99_ns']:.0f} ns ({b['n_bursts']} bursts); corr with low-freq current {b['corr_with_mean_current']:.2f}")
        top = sorted(r['tile_B2_share'].items(), key=lambda x: -x[1])[:6]
        print('  top tiles, 50-200 MHz energy :', ', '.join(f'{k} {100*v:.0f}%' for k, v in top))
        topd = sorted(r['tile_droop_share'].items(), key=lambda x: -x[1])[:6]
        print('  top tiles, droop energy      :', ', '.join(f'{k} {100*v:.0f}%' for k, v in topd))
        print('  coherence gain (energy of sum / sum of tile energies):', {k: f'{v:.2f}' for k, v in r['coherence_gain'].items()})
        pc = sorted(r['B2_pair_coherence'].items(), key=lambda x: -x[1])
        print('  50-200 MHz pair coherence: highest', ', '.join(f'{k} {v:.2f}' for k, v in pc[:4]), '| lowest', ', '.join(f'{k} {v:.2f}' for k, v in pc[-3:]))
