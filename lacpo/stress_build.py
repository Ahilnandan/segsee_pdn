"""Build replay inputs for a new workload (e.g. the di/dt stressmark) WITHOUT re-fitting the uop current kernels, so
every existing result stays unchanged: <w>_tileI.npy (tilecurrents), <w>_attr.npz (workload.build with hpos.npy),
<w>_attr5.npz (workload5.build with hpos_blocks.npy). Usage: python stress_build.py <w>"""
import sys, numpy as np
import json, pandas as pd
import tilecurrents as TC, workload as W1, workload5 as W5, cache_model as CM
B = TC.B; HERE = TC.HERE; w = sys.argv[1]
# cache model with the EXISTING coefficients (cache_model.py's main would re-fit on all benches)
coef = json.load(open(HERE + 'cache_model_coef.json'))
d = pd.read_csv(B + w + '_out/composed_P_t.csv'); Ew = {k: v.astype(float) for k, v in np.load(B + w + '_cache.npz').items()}
n = min(len(d), len(Ew['dc_rd'])); Ew = {k: v[:n] for k, v in Ew.items()}
dc = CM.apply(coef['dc'], CM.feats_dc(Ew)); ic = CM.apply(coef['ic'], CM.feats_ic(Ew)); lac = d.P_mW.values[:n]
tile = lac + dc + ic + CM.PTW_MW + CM.MISC_MW
pd.DataFrame({'Cycle': np.arange(n), 'P_lacpo_mW': lac, 'P_dcache_mW': dc, 'P_icache_mW': ic, 'P_ptw_misc_mW': CM.PTW_MW + CM.MISC_MW,
              'P_tile_mW': tile, 'I_tile_A': tile / 1.1e3}).to_csv(B + w + '_out/tile_P_t.csv', index=False)
print('cache model: D$ %.1f mW  I$ %.1f mW  LACPo %.1f mW  tile %.1f mW mean' % (dc.mean(), ic.mean(), lac.mean(), tile.mean()))
I = TC.tile_currents(w); np.save(B + w + '_tileI.npy', I); print('tileI', I.shape, 'mean total %.1f mA' % I.sum(1).mean())
h1 = np.load(HERE + 'hpos.npy'); d = W1.build(w, h1)
np.savez_compressed(B + w + '_attr.npz', A=d['A'], BG=d['BG'], cls=d['cls'], **d['inst']); print('attr  uops', len(d['cls']), 'identity err %.2e' % d['err'])
h5 = np.load(HERE + 'hpos_blocks.npy'); d = W5.build(w, h5)
np.savez_compressed(B + w + '_attr5.npz', A=d['A'], BG=d['BG'], cls=d['cls'], **d['inst']); print('attr5 uops', len(d['cls']), 'identity err %.2e' % d['err'])
