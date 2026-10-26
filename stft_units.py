"""Second pass: (a) coherence between real issue units (undo the fixed floorplan splits), (b) short-window burst statistics."""
import json, numpy as np
from scipy import signal
FS = 1e9
UNITS = {'ALU0': (4, 0.85), 'ALU1': (5, 0.85), 'ALU2': (6, 0.85), 'ALU3': (7, 0.85), 'MEM0': (8, 0.5), 'MEM1': (9, 0.5), 'FPU0': (12, 0.85), 'FPU1': (13, 0.85), 'CSR': (3, 0.7)}
def run(wl):
    It = np.load(f'It_{wl}_A.npy').astype(float); I = It.sum(1)
    U = {k: (It[:, t] - It[:, t].mean()) / fr for k, (t, fr) in UNITS.items()}
    U = {k: v for k, v in U.items() if v.std() > 0.5}
    out = {'wl': wl}
    # (a) coherence gain over units, per band, with 256-cycle STFT
    def band_E(x, lo, hi, nps=256):
        f, t, S = signal.stft(x, fs=FS, window='hann', nperseg=nps, noverlap=nps - 64, boundary=None, padded=False)
        m = (f >= lo) & (f < hi); return (np.abs(S[m]) ** 2).sum()
    tot = sum(U.values())
    out['unit_coherence_gain'] = {}
    for b, lo, hi in [('<10M', 1e5, 1e7), ('10-50M', 1e7, 5e7), ('50-200M', 5e7, 2e8), ('200-500M', 2e8, 5e8)]:
        out['unit_coherence_gain'][b] = float(band_E(tot, lo, hi) / sum(band_E(v, lo, hi) for v in U.values()))
    out['n_active_units'] = len(U)
    pc = {}
    ks = list(U)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            fc, coh = signal.coherence(U[ks[i]], U[ks[j]], fs=FS, nperseg=256)
            pc[f'{ks[i]}-{ks[j]}'] = float(coh[(fc >= 50e6) & (fc < 2e8)].mean())
    out['unit_pair_coherence_50_200'] = dict(sorted(pc.items(), key=lambda x: -x[1]))
    # (b) burst statistics with a 64-cycle window (15.6 MHz resolution), hop 16
    f, t, S = signal.stft(I - I.mean(), fs=FS, window='hann', nperseg=64, noverlap=48, boundary=None, padded=False)
    m2 = (f >= 50e6) & (f < 2e8); P2 = (np.abs(S[m2]) ** 2).sum(0)
    q = {p: float(np.percentile(P2, p)) for p in (50, 90, 99, 99.9)}
    out['B2_short'] = {'p90/p50': q[90] / q[50], 'p99/p50': q[99] / q[50], 'p99.9/p50': q[99.9] / q[50]}
    for pct in (99, 99.9):
        thr = np.percentile(P2, pct); above = P2 > thr; runs = []; r = 0
        for a in above:
            if a: r += 1
            elif r: runs.append(r); r = 0
        runs = np.array(runs) * 16
        gaps = np.diff(np.where(np.diff(above.astype(int)) == 1)[0]) * 16
        out[f'B2_events_above_p{pct}'] = dict(n=int(len(runs)), dur_median_ns=float(np.median(runs)), dur_max_ns=float(runs.max()),
                                              per_us=float(len(runs) / (len(I) / 1e3)), gap_median_ns=float(np.median(gaps)) if len(gaps) else None)
    # correlation of resonance-band power with the windowed MEAN current (program phase)
    win = 64; hop = 16
    mean_cur = np.array([I[k * hop:k * hop + win].mean() for k in range(len(P2))])
    out['corr_B2power_vs_mean_current'] = float(np.corrcoef(P2, mean_cur)[0, 1])
    return out
if __name__ == '__main__':
    R = [run(wl) for wl in ('qsort', 'dhrystone', 'mm')]
    json.dump(R, open('stft_units.json', 'w'), indent=1)
    for r in R:
        print(f"\n== {r['wl']} ({r['n_active_units']} active units)")
        print('  unit coherence gain (1 = independent, >1 = units move together):', {k: round(v, 2) for k, v in r['unit_coherence_gain'].items()})
        top = list(r['unit_pair_coherence_50_200'].items())
        print('  50-200 MHz unit-pair coherence: top', ', '.join(f'{k} {v:.2f}' for k, v in top[:4]), '| bottom', ', '.join(f'{k} {v:.2f}' for k, v in top[-3:]))
        print('  64-cycle window, 50-200 MHz power quantiles vs median:', {k: round(v, 1) for k, v in r['B2_short'].items()})
        for pct in (99, 99.9):
            e = r[f'B2_events_above_p{pct}']
            print(f"  events above p{pct}: {e['n']} ({e['per_us']:.2f}/us), duration median {e['dur_median_ns']:.0f} ns max {e['dur_max_ns']:.0f} ns, gap median {e['gap_median_ns']} ns")
        print(f"  corr(resonance-band power, windowed mean current) = {r['corr_B2power_vs_mean_current']:.2f}")
