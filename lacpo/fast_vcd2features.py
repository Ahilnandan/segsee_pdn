#!/usr/bin/env python
"""Linear-time drop-in for segseee/OutOfOrder-Optimizations/tools/vcd2features.py (Unnath's LACPo VCD -> feature CSV).

Why: the original samples every net with value_at(), which re-scans the net's whole change list for every clock
edge (O(edges x changes) per net); it did not finish on a 2.5 GB dump in 2h20m. This version reproduces the SAME
semantics in one streaming pass:
  - net naming exactly as Verilog_VCD.parse_vcd (scope path + name; multi-bit names drop the trailing [h:l] token)
  - net -> VCD mapping via the original build_vcd_index / map_lacpo_net (imported, unchanged)
  - samples at every change of the clock net to '1' (sample_at_clock), value = last change with time <= edge time
    (value_at), x/z -> 0, binary strings parsed as int
  - feature engineering via the original engineer_features (imported, unchanged), optional vectorised Hamming
    (--fast-hamming), checked bit-exact by --check against the original on a small dump.
Must run in the LACPo py37 env (imports the original modules)."""
from __future__ import print_function
import argparse, os, sys, time
import numpy as np
import pandas as pd

LACPO = os.environ.get('LACPO_ROOT', '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations')
sys.path.insert(0, os.path.join(LACPO, 'tools')); sys.path.insert(0, os.path.join(LACPO, 'common', 'vcd2csv'))
import vcd2features as V      # noqa: E402  (original, unmodified)

def read_header(path):
    """all full signal names -> code (last definition wins, as the original mapping dict does), timescale line"""
    hier, name2code = [], {}
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if '$enddefinitions' in line: break
            if '$scope' in line: hier.append(line.split()[2])
            elif '$upscope' in line: hier.pop()
            elif '$var' in line:
                ls = line.split(); size = ls[2]; code = ls[3]
                name = ''.join(ls[4:-2]) if int(size) > 1 else ''.join(ls[4:-1])
                name2code['.'.join(hier) + '.' + name] = code
    return name2code

def to_int(s):
    s = s.lower()
    if s.startswith('b'): s = s[1:]
    s = s.replace('x', '0').replace('z', '0')
    try: return int(s, 2) if s else 0
    except ValueError: return 0

def stream_sample(path, codes, clock_code, max_edges=None):
    """one pass: current value per wanted code; at each timestamp boundary, if the clock went to '1' at the
    previous timestamp, record all wanted values (includes every change at that same timestamp)"""
    want = set(codes); cur = {c: 0 for c in codes}
    rows = []; clk_rose = False; started = False
    with open(path, 'r') as fh:
        for line in fh:
            if not started:
                if '$enddefinitions' in line: started = True
                continue
            c0 = line[0] if line else ''
            if c0 == '#':
                if clk_rose:
                    rows.append([cur[c] for c in codes]); clk_rose = False
                    if max_edges and len(rows) >= max_edges: break
                continue
            if c0 in 'bBrR':
                parts = line[1:].split()
                if len(parts) != 2: continue
                val, code = parts
                if code in want: cur[code] = to_int(val)
            elif c0 in '01xXzZ':
                code = line[1:].strip()
                if code == clock_code and c0 == '1': clk_rose = True
                if code in want: cur[code] = 1 if c0 == '1' else 0
        else:
            if clk_rose: rows.append([cur[c] for c in codes])
    return np.array(rows, dtype=object) if rows else np.zeros((0, len(codes)), dtype=object)

def engineer_fast(traces, features):
    """same output as V.engineer_features; Hamming via numpy popcount when all values fit in 63 bits"""
    n = len(next(iter(traces.values()))) if traces else 0
    cols = {'Cycle': np.arange(n, dtype=np.int32)}; cache = {}
    def series(net, slc):
        key = (net, slc)
        if key not in cache:
            raw = traces.get(net)
            cache[key] = np.zeros(n, dtype=np.int64) if raw is None else np.array([V.slice_value(v, slc) for v in raw], dtype=object)
        return cache[key]
    def popcount(x):
        x = x.astype(np.uint64); c = np.zeros(len(x), dtype=np.int64)
        for sh in range(0, 64, 8): c += np.array([bin(b).count('1') for b in range(256)], dtype=np.int64)[((x >> np.uint64(sh)) & np.uint64(255)).astype(np.int64)]
        return c
    def delay(arr, k):
        out = np.zeros_like(arr)
        if k < len(arr): out[k:] = arr[:-k]
        return out
    for feat in features:
        net, slc = V.split_slice(V.strip_engineered(feat))
        vals = series(net, slc)
        if vals.dtype == object and len(vals) and max(int(v).bit_length() for v in vals) <= 63: vals_i = vals.astype(np.int64)
        elif vals.dtype != object: vals_i = vals
        else: vals_i = None
        hd = np.zeros(n, dtype=np.int64)
        if n:
            if vals_i is not None: hd[1:] = popcount(np.bitwise_xor(vals_i[1:], vals_i[:-1]))
            else: hd[1:] = np.array([V.hamming(vals[i - 1], vals[i]) for i in range(1, n)], dtype=np.int64)
        v_out = vals_i if vals_i is not None else vals
        if feat.startswith('v2_'): cols[feat] = delay(v_out, 2)
        elif feat.startswith('v1_'): cols[feat] = delay(v_out, 1)
        elif feat.startswith('d2_'): cols[feat] = delay(hd, 1)
        elif feat.startswith('d1_') or feat.startswith('d_'): cols[feat] = hd
        else: cols[feat] = v_out
    return pd.DataFrame(cols)

def vcd_to_features_fast(vcd_path, features, clock_net, fast_hamming=True, max_edges=None):
    nets = V.needed_nets(features)
    name2code = read_header(vcd_path)
    rel_to_full, by_leaf, by2 = V.build_vcd_index(list(name2code.keys()))
    rename = {}
    for nn in [clock_net] + list(nets):
        m = V.map_lacpo_net(nn, rel_to_full, by_leaf, by2)
        if m: rename[nn] = m
    clock_vcd = rename.get(clock_net) or rel_to_full.get('lsu.clock')
    if not clock_vcd: V.die('Clock %s not in VCD after remap. Have %d boom_tile nets.' % (clock_net, len(rel_to_full)))
    rename[clock_net] = clock_vcd
    siglist = sorted(set(rename.values()))
    print('vcd map: %d/%d nets (+clock) matched' % (len(rename) - 1, len(nets)))
    print('clock', clock_net, '->', clock_vcd)
    codes = sorted(set(name2code[s] for s in siglist))
    S = stream_sample(vcd_path, codes, name2code[clock_vcd], max_edges)
    col = {c: i for i, c in enumerate(codes)}
    traces = {}
    for nn in nets:
        src = rename.get(nn)
        traces[nn] = list(S[:, col[name2code[src]]]) if src else [0] * len(S)
    df = engineer_fast(traces, features) if fast_hamming else V.engineer_features(traces, features)
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcd', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--clock', default='TestDriver.testHarness.TestHarness.boom_tile.core.lsu.clock')
    ap.add_argument('--benchmark', default=''); ap.add_argument('--max-edges', type=int, default=None)
    ap.add_argument('--check', action='store_true', help='also run the ORIGINAL vcd2features path and assert identical output')
    a = ap.parse_args()
    feats = V.collect_invoke_features()
    t0 = time.time(); df = vcd_to_features_fast(a.vcd, feats, a.clock, True, a.max_edges)
    print('fast: cycles %d columns %d in %.1fs' % (len(df), len(df.columns), time.time() - t0))
    if a.check:
        t0 = time.time(); ref, _ = V.vcd_to_features(a.vcd, feats, a.clock, 3)
        print('original: cycles %d columns %d in %.1fs' % (len(ref), len(ref.columns), time.time() - t0))
        assert list(ref.columns) == list(df.columns), 'column mismatch'
        assert len(ref) == len(df), 'row count mismatch %d vs %d' % (len(ref), len(df))
        bad = [c for c in ref.columns if not np.array_equal(np.asarray(ref[c], dtype=object).astype(str), np.asarray(df[c], dtype=object).astype(str))]
        print('IDENTICAL' if not bad else 'MISMATCH in %d columns, e.g. %s' % (len(bad), bad[:5]))
    if a.benchmark: df['Benchmark'] = a.benchmark
    df.to_csv(a.out, index=False); print('wrote', a.out)

if __name__ == '__main__':
    main()
