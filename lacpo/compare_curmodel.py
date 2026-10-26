"""LACPo (exact training core, MediumBoom, 3 ns) vs the old curmodel traces (PaperMatchBoomV3 replay, 1 ns assumed)
on the kernel phase of qsort / dhrystone / mm. Different core + clock, so compare level, swing, per-cycle steps,
spectrum vs clock-normalised frequency, and per-unit split -- not cycle-by-cycle."""
from __future__ import print_function
import json, numpy as np, pandas as pd
PD = '/home/ahilnandan/Documents/PMIC/pdn_work/'
LB = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
OLD_KERNEL = {'qsort': (40, 235), 'dhrystone': (2, 67), 'mm': (7, 77)}      # kilo-cycles, from current_plots.py
OLD_GROUPS = {'frontend+decode': [0, 1], 'rename/ROB/RF/issue': [2], 'CSR': [3], 'int exec': [4, 5, 6, 7],
              'LSU+D$': [8, 9, 10, 11], 'FP': [12, 13], 'other/L2if': [14, 15]}
def lacpo_group(b):
    if b.startswith('frontend') or 'decode' in b: return 'frontend+decode'
    if 'fp_' in b: return 'FP'
    if 'csr' in b: return 'CSR'
    if 'ALUExeUnit' in b: return 'int exec'
    if 'lsu' in b: return 'LSU+D$'
    if 'glue' in b: return 'other/L2if'
    return 'rename/ROB/RF/issue'
BANDS = [(0, .002), (.002, .02), (.02, .05), (.05, .15), (.15, .3), (.3, .5)]   # cycles/cycle

def stats(I, fclk_mhz):
    d = np.abs(np.diff(I)); x = I - I.mean()
    F = np.abs(np.fft.rfft(x)) ** 2; f = np.fft.rfftfreq(len(x)); F[0] = 0; tot = F.sum()
    band = [100 * F[(f >= a) & (f < b)].sum() / tot for a, b in BANDS]
    absb = 100 * F[(f * fclk_mhz >= 50) & (f * fclk_mhz < 200)].sum() / tot
    # strongest tones above 1/200 cycles
    m = f > 0.005; idx = np.argsort(F[m])[::-1]; ff = f[m][idx]; tones = []
    for v in ff:
        if all(abs(1 / v - 1 / t) > 0.3 for t in tones): tones.append(v)
        if len(tones) == 4: break
    return dict(mean=I.mean(), std=I.std(), min=I.min(), p01=np.percentile(I, .1), p999=np.percentile(I, 99.9), max=I.max(),
                dyn_frac=100 * (I.mean() - np.percentile(I, .1)) / I.mean(),
                dI_p999=np.percentile(d, 99.9), dI_max=d.max(), dI_rel=100 * np.percentile(d, 99.9) / I.mean(),
                bands=band, band_50_200MHz=absb, tone_periods=[1 / t for t in tones])

out = {}
for w in ('qsort', 'dhrystone', 'mm'):
    It = np.load(PD + 'It_%s_A.npy' % w).astype(float); k0, k1 = OLD_KERNEL[w]; It = It[k0 * 1000:k1 * 1000]
    old = stats(It.sum(1), 1000.0)
    old['units'] = {g: It[:, ix].sum(1).mean() for g, ix in OLD_GROUPS.items()}
    old['units_std'] = {g: It[:, ix].sum(1).std() for g, ix in OLD_GROUPS.items()}
    d = pd.read_csv(LB + w + '_out/composed_P_t.csv'); I = d.I_A.values * 1e3
    s = int(np.argmax(pd.Series(I).rolling(2000).mean().values > 90)) ; e = len(I) - 2000    # skip HTIF load + tail
    blocks = [c for c in d.columns if c not in ('Cycle', 'time_ns', 'activity', 'P_mW', 'I_A')]
    new = stats(I[s:e], 333.3)
    G = {}
    for b in blocks: G.setdefault(lacpo_group(b), []).append(b)
    new['units'] = {g: d[bs].values[s:e].sum(1).mean() / 1.1 for g, bs in G.items()}
    new['units_std'] = {g: d[bs].values[s:e].sum(1).std() / 1.1 for g, bs in G.items()}
    new['kernel'] = (s, e)
    out[w] = dict(old=old, lacpo=new)

def row(name, k, fmt='%7.1f'):
    print('%-26s' % name + ''.join((fmt % out[w][src][k]) for w in out for src in ('old', 'lacpo')))
print('%-26s' % '' + ''.join('%7s' % ('old' if s == 'old' else 'LACPo') for w in out for s in ('old', 'lacpo')))
print('%-26s' % 'workload' + ''.join('%7s' % w[:6] for w in out for s in (0, 1)))
for n, k in [('mean mA', 'mean'), ('std mA', 'std'), ('min mA', 'min'), ('p0.1 mA', 'p01'), ('p99.9 mA', 'p999'), ('max mA', 'max'),
             ('dynamic % of mean', 'dyn_frac'), ('|dI|/cycle p99.9 mA', 'dI_p999'), ('|dI|/cycle max mA', 'dI_max'),
             ('|dI| p99.9 % of mean', 'dI_rel'), ('% var in 50-200 MHz', 'band_50_200MHz')]:
    row(n, k)
print('\n% of variance by clock-normalised band (cycles/cycle):')
for i, (a, b) in enumerate(BANDS):
    print('%-26s' % ('  %.3f-%.3f' % (a, b)) + ''.join('%7.1f' % out[w][s]['bands'][i] for w in out for s in ('old', 'lacpo')))
print('\nstrongest tone periods (cycles):')
for w in out:
    for s in ('old', 'lacpo'): print('  %-10s %-6s' % (w, s), ' '.join('%.1f' % p for p in out[w][s]['tone_periods']))
print('\nper-unit mean mA (std mA):')
for g in OLD_GROUPS:
    print('%-26s' % g + ''.join('%7.1f' % out[w][s]['units'].get(g, 0) for w in out for s in ('old', 'lacpo')),
          ' | std' + ''.join('%6.1f' % out[w][s]['units_std'].get(g, 0) for w in out for s in ('old', 'lacpo')))
json.dump(out, open(PD + 'lacpo/compare_curmodel.json', 'w'), indent=1, default=float)
