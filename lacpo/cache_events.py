"""Per-cycle cache activity from a BOOM VCD, sampled exactly like the LACPo features (clock rise, values after all
changes at that timestamp): D$ / I$ array enables (with way counts), refill/writeback beats, pipeline valids, and
Hamming activity on the array read-data buses. Output: <out>.npz with one int array per signal name."""
import argparse, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_vcd2features import read_header, stream_sample

T = 'TOP.TestHarness.TestHarness.boom_tile.'
SIG = {  # short name -> (vcd leaf under boom_tile, reduce: 'v' value, 'pop' popcount)
 'dc_rd': ('dcache.data_io_read_valid', 'v'), 'dc_rd_ways': ('dcache.data_io_read_bits_way_en', 'pop'),
 'dc_wr': ('dcache.data_io_write_valid', 'v'), 'dc_wr_ways': ('dcache.data_io_write_bits_way_en', 'pop'),
 'dc_mrd': ('dcache.meta_io_read_valid', 'v'), 'dc_mwr': ('dcache.meta_io_write_valid', 'v'),
 'dc_refill': ('dcache.auto_out_d_valid', 'v'), 'dc_wb': ('dcache.auto_out_c_valid', 'v'), 'dc_acq': ('dcache.auto_out_a_valid', 'v'),
 'dc_req': ('dcache.io_cpu_req_valid', 'v'), 'dc_s1': ('dcache.s1_valid', 'v'), 'dc_s2': ('dcache.s2_valid', 'v'),
 'dc_tlb': ('dcache.dtlb_io_req_valid', 'v'), 'dc_mshr': ('dcache.mshrs_io_req_valid', 'v'),
 'dc_s2d': ('dcache.s2_data_muxed', 'raw'),
 'dc_dbus': ('dcache.auto_out_d_bits_data', 'raw'), 'dc_dbusv': ('dcache.auto_out_d_valid', 'v'),
 'ic_dbus': ('frontend.icache.auto_master_out_d_bits_data', 'raw'), 'ic_dbusv': ('frontend.icache.auto_master_out_d_valid', 'v'),
 'ic_s0': ('frontend.icache.s0_valid', 'v'), 'ic_s1': ('frontend.icache.s1_valid', 'v'), 'ic_s2': ('frontend.icache.s2_valid', 'v'),
 'ic_w0': ('frontend.icache.dataArrayWay_0_RW0_en', 'v'), 'ic_w1': ('frontend.icache.dataArrayWay_1_RW0_en', 'v'),
 'ic_w2': ('frontend.icache.dataArrayWay_2_RW0_en', 'v'), 'ic_w3': ('frontend.icache.dataArrayWay_3_RW0_en', 'v'),
 'ic_tag': ('frontend.icache.tag_array_RW0_en', 'v'), 'ic_refill': ('frontend.icache.refill_one_beat', 'v'),
 'ic_wen': ('frontend.icache.wen', 'v'), 'ic_kill': ('frontend.icache.io_s1_kill', 'v'), 'ic_miss': ('frontend.icache.s2_miss', 'v'),
 'ic_d0': ('frontend.icache.s1_dout_0', 'raw'), 'ic_d1': ('frontend.icache.s1_dout_1', 'raw'),
 'ic_d2': ('frontend.icache.s1_dout_2', 'raw'), 'ic_d3': ('frontend.icache.s1_dout_3', 'raw'),
}
POP8 = np.array([bin(b).count('1') for b in range(256)], dtype=np.int64)
def popcount(x):
    x = x.astype(np.uint64); c = np.zeros(len(x), dtype=np.int64)
    for sh in range(0, 64, 8): c += POP8[((x >> np.uint64(sh)) & np.uint64(255)).astype(np.int64)]
    return c

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--vcd', required=True); ap.add_argument('--out', required=True); a = ap.parse_args()
    n2c = read_header(a.vcd)
    def full(leaf):
        if T + leaf in n2c: return T + leaf
        m = [k for k in n2c if k.endswith('boom_tile.' + leaf)]; return m[0] if m else None
    names = {k: full(v[0]) for k, v in SIG.items()}
    miss = [k for k, v in names.items() if v is None]
    if miss: print('MISSING', miss)
    clk = full('core.lsu.clock'); codes = sorted(set(n2c[v] for v in names.values() if v))
    S = stream_sample(a.vcd, codes, n2c[clk]); col = {c: i for i, c in enumerate(codes)}
    out = {}
    for k, (leaf, red) in SIG.items():
        if names[k] is None: continue
        x = np.array([int(v) & ((1 << 64) - 1) for v in S[:, col[n2c[names[k]]]]], dtype=np.uint64)
        if red == 'pop': out[k] = popcount(x)
        elif red == 'raw':                         # Hamming distance to previous cycle (data-bus activity)
            h = np.zeros(len(x), dtype=np.int64); h[1:] = popcount(np.bitwise_xor(x[1:], x[:-1])); out[k + '_hd'] = h
        else: out[k] = x.astype(np.int64)
    np.savez_compressed(a.out, **out); print('cycles', len(S), 'signals', len(out), '->', a.out)

if __name__ == '__main__':
    main()
