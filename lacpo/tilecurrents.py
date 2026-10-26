"""Phase 0.7: per-tile current traces (mA, per 3 ns cycle) from per-block LACPo currents + cache model, mapped through
floorplan.json (block current spread over its tiles by area). Also per-tile leakage (mA) for the PDN intra-cycle shape.
Baseline schedule only; the closed-loop simulator re-maps issue-coupled current per uop (Phase 1)."""
import json, sys, numpy as np, pandas as pd
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
HERE = '/home/ahilnandan/Documents/PMIC/pdn_work/lacpo/'
FP = json.load(open(HERE + 'floorplan.json')); N = FP['n']; NT = N * N
COLMAP = {'frontend__bpdpipeline__bpd__counter_table_1': 'bpd', 'frontend__bpdpipeline__btb_handpicked_signals_1': 'btb',
          'frontend__fetch_controller_1': 'fetch', 'core__decode_units_0_1': 'decode', 'core__decode_units_1_1': 'decode',
          'core__rename_stage__freelist_history_1_1': 'rename', 'core__rename_stage___maptable_history_1_1': 'rename',
          'core__fp_rename_stage__freelist_1': 'fp_rename', 'core__fp_rename_stage__maptable_1': 'fp_rename',
          'core__rob_handpicked_signals_1': 'rob', 'core__int_issue_unit_handpicked_signals_3_1': 'int_iq', 'core__iregfile_1': 'iregfile',
          'core__iregister_read_handpicked_signals_1': 'iregread', 'core__ALUExeUnit_handpicked_signals_1': 'alu0',
          'core__csr_exe_unit_1': 'alu1_csr', 'core__csr_1': 'alu1_csr', 'core__fp_pipeline_handpicked_signals_2_1': 'fp',
          'core__mem_issue_unit_1': 'mem_iq', 'core__lsu_1': 'lsu', 'core_glue_1': 'core_glue'}
DSPLIT = {'dcache_data': 0.88, 'dcache_meta': 0.08, 'dcache_ctrl': 0.04}     # PTPX median total-power split of the D$
def weights():
    W = {}
    for b, d in FP['blocks'].items():
        w = np.zeros(NT)
        for k, f in d['tiles'].items():
            c, r = map(int, k.split('_')); w[r * N + c] = f
        W[b] = w
    return W
W = weights()
LEAK_T = sum(FP['blocks'][b]['leak_mW'] * W[b] for b in W) / 1.1           # mA per tile
def block_currents(w):
    d = pd.read_csv(B + w + '_out/composed_P_t.csv'); t = pd.read_csv(B + w + '_out/tile_P_t.csv'); n = min(len(d), len(t))
    blk = {}
    for c, b in COLMAP.items(): blk[b] = blk.get(b, 0) + d[c].values[:n]
    blk['icache'] = t.P_icache_mW.values[:n]
    for b, f in DSPLIT.items(): blk[b] = f * t.P_dcache_mW.values[:n]
    blk['ptw'] = np.full(n, 7.24); blk['tile_misc'] = np.full(n, 1.5)
    return {b: v / 1.1 for b, v in blk.items()}, n                          # mA
def tile_currents(w):
    blk, n = block_currents(w); I = np.zeros((n, NT), np.float32)
    for b, v in blk.items(): I += np.outer(v, W[b]).astype(np.float32)
    return I
if __name__ == '__main__':
    WL = sys.argv[1:] or ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
    for w in WL:
        I = tile_currents(w); np.save(B + w + '_tileI.npy', I)
        tot = I.sum(1); s = int(np.argmax(pd.Series(tot).rolling(100).mean().fillna(0).values * (np.arange(len(tot)) > 3000) > 150)); X = I[s:]
        dI = np.abs(np.diff(X, axis=0))
        print('%-10s cycles %6d  tile mean max %.1f mA (tile %d)  tile p99.9 max %.1f mA (tile %d)  tile step p99.9 max %.1f mA (tile %d)  total mean %.1f'
              % (w, len(X), X.mean(0).max(), X.mean(0).argmax(), np.percentile(X, 99.9, axis=0).max(), np.percentile(X, 99.9, axis=0).argmax(),
                 np.percentile(dI, 99.9, axis=0).max(), np.percentile(dI, 99.9, axis=0).argmax(), tot[s:].mean()))
    np.save(HERE + 'tile_leak_mA.npy', LEAK_T)
