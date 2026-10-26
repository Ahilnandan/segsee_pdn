"""Calibrated cache current model (D$, I$) + constant PTW/misc, added to LACPo -> full BoomTile power per cycle.

No aligned training data exists (only PTPX time-based power for THEIR median binary; no VCD for it), so each cache is
  P[n] = floor + sum_k e_k * event_k[n]          (e_k >= 0, events from our VCD, sampled like LACPo features)
fitted on median ONLY by matching the PTPX execution-phase per-cycle distribution (quantiles, tail-weighted) and mean.
Held-out check: whole-run PTPX averages for the other benchmarks (dcache, tile) from power_anal/top_csvs.merged.csv."""
from __future__ import print_function
import json, numpy as np, pandas as pd
from scipy.optimize import minimize
R = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/'
B = R + 'traces/lacpo_cy100/bench/'
HERE = '/home/ahilnandan/Documents/PMIC/pdn_work/lacpo/'
PTPX_EXEC0 = 12534                       # PTPX median: HTIF load before this cycle
WL = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
PTW_MW, MISC_MW = 7.24, 1.5              # PTW flat in PTPX (std 0.2 mW, 7.2-7.3 on all 10 benches); TL buffers/arbs/frontend leftovers

def lag(x, k=1):
    y = np.zeros_like(x); y[k:] = x[:-k]; return y
def dbus(E):   # toggles on the TileLink D-channel data wires shared by both caches (fan-out into both flop arrays)
    return np.maximum(E.get('dc_dbus_hd', 0) * 1.0, E.get('ic_dbus_hd', 0) * 1.0) / 64.0
def feats_dc(E):
    rd, wr = (E['dc_rd'] > 0).astype(float), (E['dc_wr'] > 0).astype(float)
    return {'rd': rd, 'wr': wr, 'rd_wr': rd * wr, 'refill': (E['dc_refill'] > 0) * 1.0, 'wb': (E['dc_wb'] > 0) * 1.0,
            'meta_rd': (E['dc_mrd'] > 0) * 1.0, 'meta_wr': (E['dc_mwr'] > 0) * 1.0, 'rd_data_hd': E['dc_s2d_hd'] / 64.0,
            'wr_next': lag(wr), 'refill_next': lag((E['dc_refill'] > 0) * 1.0), 'dbus_hd': dbus(E), 'dbus_hd_next': lag(dbus(E))}
def feats_ic(E):
    ways = (E['ic_w0'] + E['ic_w1'] + E['ic_w2'] + E['ic_w3']) / 4.0
    hd = (E['ic_d0_hd'] + E['ic_d1_hd'] + E['ic_d2_hd'] + E['ic_d3_hd']) / 256.0
    return {'ways_rd': ways, 'tag_rd': E['ic_tag'] * 1.0, 'dout_hd': hd, 'refill': E['ic_refill'] * 1.0, 'wen': E['ic_wen'] * 1.0,
            'refill_next': lag(E['ic_refill'] * 1.0), 'miss': E['ic_miss'] * 1.0, 'dbus_hd': dbus(E), 'dbus_hd_next': lag(dbus(E))}

Q = np.r_[np.linspace(1, 90, 30), np.linspace(90.5, 99, 18), np.linspace(99.1, 99.9, 9)]
WQ = np.r_[np.ones(30), 3 * np.ones(18), 6 * np.ones(9)]
def fit(F, target, floor):
    names = sorted(F); X = np.stack([F[k] for k in names], 1)
    qt = np.percentile(target, Q); mt = target.mean()
    def model(p): return floor + X.dot(p[1:])
    QD = np.r_[np.linspace(90, 99, 10), np.linspace(99.1, 99.9, 9)]; qd = np.percentile(np.abs(np.diff(target)), QD)
    def loss(p):   # level quantiles + per-cycle step quantiles (droop is driven by steps) + mean
        y = model(p); dq = np.percentile(np.abs(np.diff(y)), QD)
        return (WQ * (np.percentile(y, Q) - qt) ** 2).sum() / WQ.sum() + ((dq - qd) ** 2).mean() + 4 * (y.mean() - mt) ** 2
    best = None
    for seed in range(12):                 # quantile loss is non-convex: multistart
        rng = np.random.RandomState(seed); p0 = np.r_[np.percentile(target, 5), rng.uniform(0, 60, len(names))]; p0[0] = floor
        r = minimize(loss, p0, method='L-BFGS-B', bounds=[(floor, floor)] + [(0, None)] * len(names))
        if best is None or r.fun < best.fun: best = r
    p = best.x; return dict(floor=p[0], **{k: v for k, v in zip(names, p[1:])}), best.fun ** 0.5

def apply(coef, F):
    y = np.full(len(next(iter(F.values()))), coef['floor'])
    for k, v in coef.items():
        if k != 'floor': y = y + v * F[k]
    return y

def describe(x):
    d = np.abs(np.diff(x)); return dict(mean=x.mean(), std=x.std(), p50=np.median(x), p99=np.percentile(x, 99), p999=np.percentile(x, 99.9),
                                        max=x.max(), dP_p999=np.percentile(d, 99.9), frac_gt25=(x > np.median(x) + 25).mean())

if __name__ == '__main__':
    z = np.load(HERE + 'ptpx_median.npz'); P = z['P'].astype(float)[1:] * 1e3; N = list(z['names'])
    tgt = {'dc': P[PTPX_EXEC0:, N.index('BoomTile/dcache')], 'ic': P[PTPX_EXEC0:, N.index('BoomTile/frontend/icache')]}
    Lm = pd.read_csv(B + 'median_out/composed_P_t.csv').P_mW.values
    m = pd.Series(Lm).rolling(100).mean().values; m = np.where(np.isnan(m), 0, m); m[:3000] = 0; s0 = int(np.argmax(m > 90)) - 100
    E = dict(np.load(B + 'median_cache.npz')); E = {k: v[s0:].astype(float) for k, v in E.items()}
    coef = {}
    idle = {'dc': np.median(P[1000:PTPX_EXEC0 - 500, N.index('BoomTile/dcache')]), 'ic': np.median(P[1000:PTPX_EXEC0 - 500, N.index('BoomTile/frontend/icache')])}
    print('idle floors from PTPX load phase (leakage + clock): D$ %.1f  I$ %.1f mW' % (idle['dc'], idle['ic']))
    for c, ff in (('dc', feats_dc), ('ic', feats_ic)):
        coef[c], rms = fit(ff(E), tgt[c], idle[c]); y = apply(coef[c], ff(E))
        print('\n%s fit on median exec (rms quantile err %.2f mW); coefficients mW (pJ = x3):' % (c, rms))
        print('  ' + '  '.join('%s %.1f' % (k, v) for k, v in sorted(coef[c].items(), key=lambda kv: -kv[1]) if v > 0.05))
        a, b = describe(tgt[c]), describe(y)
        print('  %-8s' % '' + ''.join('%9s' % k for k in a)); print('  %-8s' % 'PTPX' + ''.join('%9.2f' % v for v in a.values()))
        print('  %-8s' % 'model' + ''.join('%9.2f' % v for v in b.values()))
    json.dump(coef, open(HERE + 'cache_model_coef.json', 'w'), indent=1)
    # full tile for every benchmark + held-out check against PTPX whole-run averages
    names = [l.strip().split('.')[0] for l in open(R + 'BOOM/power_anal/list_of_top_csvs')]
    avg = {}
    for l in open(R + 'BOOM/power_anal/top_csvs.merged.csv'):
        p = [x.strip() for x in l.split(',')]
        avg[p[0]] = {n: float(p[2 + 4 * i + 3]) * 1e3 for i, n in enumerate(names)}
    print('\nwhole-run mW: ours (LACPo + cache model) vs PTPX averages  [median = fit set, others held out]')
    print('%-10s %8s %8s %8s | %8s %8s %8s | %8s %8s' % ('bench', 'D$ ours', 'D$ PTPX', 'err%', 'tile ours', 'tile PTX', 'err%', 'I$ ours', 'LACPo'))
    rows = []
    for w in WL:
        d = pd.read_csv(B + w + '_out/composed_P_t.csv'); Ew = {k: v.astype(float) for k, v in np.load(B + w + '_cache.npz').items()}
        n = min(len(d), len(Ew['dc_rd'])); Ew = {k: v[:n] for k, v in Ew.items()}
        dc = apply(coef['dc'], feats_dc(Ew)); ic = apply(coef['ic'], feats_ic(Ew)); lac = d.P_mW.values[:n]
        tile = lac + dc + ic + PTW_MW + MISC_MW
        out = pd.DataFrame({'Cycle': np.arange(n), 'P_lacpo_mW': lac, 'P_dcache_mW': dc, 'P_icache_mW': ic,
                            'P_ptw_misc_mW': PTW_MW + MISC_MW, 'P_tile_mW': tile, 'I_tile_A': tile / 1.1e3})
        out.to_csv(B + w + '_out/tile_P_t.csv', index=False)
        if w in avg['BoomTile']:
            ed = 100 * (dc.mean() / avg['dcache'][w] - 1); et = 100 * (tile.mean() / avg['BoomTile'][w] - 1)
            rows.append((w, ed, et))
            print('%-10s %8.1f %8.1f %+8.1f | %8.1f %8.1f %+8.1f | %8.1f %8.1f' % (w, dc.mean(), avg['dcache'][w], ed, tile.mean(), avg['BoomTile'][w], et, ic.mean(), lac.mean()))
        else:
            print('%-10s %8.1f %8s %8s | %8.1f %8s %8s | %8.1f %8.1f' % (w, dc.mean(), '-', '', tile.mean(), '-', '', ic.mean(), lac.mean()))
    ho = [r for r in rows if r[0] != 'median']
    print('held-out mean |err|: D$ %.1f%%  tile %.1f%%  (n=%d)' % (np.mean([abs(r[1]) for r in ho]), np.mean([abs(r[2]) for r in ho]), len(ho)))
