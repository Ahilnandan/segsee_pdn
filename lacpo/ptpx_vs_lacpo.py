"""Ground truth (PTPX time-based, median, 3 ns) vs LACPo on our median run: what LACPo covers, what the caches add."""
from __future__ import print_function
import numpy as np, pandas as pd
z = np.load('ptpx_median.npz'); P = z['P'].astype(float) * 1e3; N = list(z['names'])      # mW
g = lambda n: P[:, N.index(n)]
tile = g('BoomTile'); core = g('BoomTile/core'); fe = g('BoomTile/frontend'); ic = g('BoomTile/frontend/icache')
dc = g('BoomTile/dcache'); ptw = g('BoomTile/ptw'); fc = g('BoomTile/frontend/fetch_controller'); bpd = g('BoomTile/frontend/bpdpipeline')
cov = core + fc + bpd; miss = tile - cov
L = pd.read_csv('/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/median_out/composed_P_t.csv').P_mW.values
def act(x):   # drop HTIF load (sustained low) like the comparison script
    s = int(np.argmax(pd.Series(x).rolling(500).mean().values > np.percentile(x, 5) + 15)); return x[s:]
def st(x, nm):
    d = np.abs(np.diff(x)); F = np.abs(np.fft.rfft(x - x.mean()))**2; f = np.fft.rfftfreq(len(x)); F[0] = 0
    b = lambda a, c: 100 * F[(f >= a) & (f < c)].sum() / F.sum()
    print('%-28s n%6d mean %6.1f std %5.1f min %6.1f p99.9 %6.1f max %6.1f |dP|p99.9 %5.1f | var%% <.02 %4.0f .02-.15 %4.0f .15-.5 %4.0f'
          % (nm, len(x), x.mean(), x.std(), x.min(), np.percentile(x, 99.9), x.max(), np.percentile(d, 99.9), b(0, .02), b(.02, .15), b(.15, .51)))
print('mW, per 3 ns cycle')
st(tile, 'PTPX tile (whole run)'); st(cov, 'PTPX LACPo-covered part'); st(miss, 'PTPX NOT covered (caches..)')
st(dc, '  PTPX dcache'); st(ic, '  PTPX icache'); st(ptw, '  PTPX ptw')
st(act(L), 'LACPo our median (active)'); st(L, 'LACPo our median (whole)')
print('corr(covered, not covered) per cycle: %.2f' % np.corrcoef(cov, miss)[0, 1])
print('var(tile) = var(cov) %.0f + var(miss) %.0f + 2cov %.0f' % (cov.var(), miss.var(), 2 * np.cov(cov, miss)[0, 1]))
