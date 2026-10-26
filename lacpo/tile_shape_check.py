"""Per-cycle shape of the full-tile model (LACPo + cache model) vs PTPX ground truth, median execution phase."""
import numpy as np, pandas as pd
z = np.load('ptpx_median.npz'); P = z['P'].astype(float)[1:] * 1e3; N = list(z['names'])
g = lambda n: P[12534:, N.index(n)]
d = pd.read_csv('/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/median_out/tile_P_t.csv')[19640:]
def st(x):
    dd = np.abs(np.diff(x)); F = np.abs(np.fft.rfft(x - x.mean())) ** 2; f = np.fft.rfftfreq(len(x)); F[0] = 0
    b = lambda a, c: 100 * F[(f >= a) & (f < c)].sum() / F.sum()
    return [x.mean(), x.std(), np.percentile(x, 99), np.percentile(x, 99.9), x.max(), np.percentile(dd, 99), np.percentile(dd, 99.9), b(0, .02), b(.02, .15), b(.15, .51)]
h = ['mean', 'std', 'p99', 'p99.9', 'max', '|dP|p99', '|dP|p99.9', 'var<.02', '.02-.15', '.15-.5']
print('%-22s' % 'median exec, mW' + ''.join('%9s' % k for k in h))
for nm, x in (('PTPX tile (truth)', g('BoomTile')), ('ours: LACPo+caches', d.P_tile_mW.values), ('LACPo only', d.P_lacpo_mW.values)):
    print('%-22s' % nm + ''.join('%9.1f' % v for v in st(x)))
for nm, a, b in (('PTPX', g('BoomTile/dcache'), g('BoomTile/frontend/icache')), ('ours', d.P_dcache_mW.values, d.P_icache_mW.values)):
    A = a > np.median(a) + 25; B = b > np.median(b) + 25
    print('%-5s D$ spike %.3f  I$ spike %.3f  both %.4f (independent: %.4f)' % (nm, A.mean(), B.mean(), (A & B).mean(), A.mean() * B.mean()))
