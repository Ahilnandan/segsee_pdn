"""Phase 0: block list, area shares (PTPX leakage shares, same library), 8x8 tile floorplan, block->tile weights.

Floorplan (columns 0..7, rows 0..7), contiguous column-major fill in this order:
  I$ (cols 0-1+) | BPD/BTB/fetch | decode/rename/ROB/IQs | int RF/read + ALU0 + ALU1/CSR | FP pipeline (bottom middle)
  | LSU + mem IQ | D$ meta/MSHR | D$ data (cols 6-7) | PTW
Output: floorplan.json with blocks (area share, mm^2, tiles {tile: fraction}) and the tile map."""
import json, re, numpy as np
RPT = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/BOOM/power_anal/median.power_hier.avg.rpt'
DIE_MM2 = 4.0; N = 8
def leaks():
    out = {}
    for line in open(RPT).read().splitlines()[14:]:
        m = re.match(r'^(\s*)(\S+)\s+\((\S+)\)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if not m: continue
        try: out.setdefault((len(m.group(1)) // 2, m.group(2)), float(m.group(6)) * 1e3)
        except ValueError: pass
    return out
L = leaks(); g = lambda d, n: L.get((d, n), 0.0)
core = g(1, 'core'); fe = g(1, 'frontend')
named_core = {'rename': g(2, 'rename_stage'), 'fp_rename': g(2, 'fp_rename_stage'), 'rob': g(2, 'rob'), 'int_iq': g(2, 'int_issue_unit'),
              'mem_iq': g(2, 'mem_issue_unit'), 'iregfile': g(2, 'iregfile'), 'iregread': g(2, 'iregister_read'), 'alu0': g(2, 'ALUExeUnit'),
              'alu1_csr': g(2, 'csr_exe_unit') + g(2, 'csr'), 'lsu': g(2, 'lsu'), 'fp': g(2, 'fp_pipeline')}
blocks = [  # (name, leakage mW) in placement order
    ('icache', g(2, 'icache')), ('bpd', g(3, 'bpd')), ('btb', g(3, 'btb')), ('fetch', g(2, 'fetch_controller') + max(0, fe - g(2, 'icache') - g(2, 'bpdpipeline') - g(2, 'fetch_controller'))),
    ('decode', 0.0), ('rename', named_core['rename']), ('fp_rename', named_core['fp_rename']), ('rob', named_core['rob']),
    ('int_iq', named_core['int_iq']), ('iregfile', named_core['iregfile']), ('iregread', named_core['iregread']),
    ('alu0', named_core['alu0']), ('alu1_csr', named_core['alu1_csr']), ('fp', named_core['fp']),
    ('mem_iq', named_core['mem_iq']), ('lsu', named_core['lsu']), ('core_glue', core - sum(named_core.values())),
    ('dcache_meta', g(2, 'meta') + g(2, 'mshrs')), ('dcache_data', g(2, 'data')), ('dcache_ctrl', g(1, 'dcache') - g(2, 'meta') - g(2, 'mshrs') - g(2, 'data')),
    ('ptw', g(1, 'ptw')), ('tile_misc', 0.4)]
blocks = [(n, max(v, 0.15)) for n, v in blocks]          # decode etc. get a minimum footprint
tot = sum(v for _, v in blocks); share = {n: v / tot for n, v in blocks}
# column-major fill over 64 tiles; fractional tiles shared between consecutive blocks
order = [n for n, _ in blocks]; tiles = [(c, r) for c in range(N) for r in range(N)]   # column 0 top->bottom, then column 1...
need = [share[n] * N * N for n in order]; ti = 0; rem = 1.0; assign = {n: {} for n in order}
for n, q in zip(order, need):
    while q > 1e-9 and ti < len(tiles):
        take = min(q, rem); key = '%d_%d' % tiles[ti]; assign[n][key] = assign[n].get(key, 0) + take
        q -= take; rem -= take
        if rem <= 1e-9: ti += 1; rem = 1.0
for n in order:
    s = sum(assign[n].values()); assign[n] = {k: v / s for k, v in assign[n].items()}
# core glue (routing / leftover logic of the core) is physically spread over all core blocks, not packed in one spot
CORE = ['decode', 'rename', 'fp_rename', 'rob', 'int_iq', 'iregfile', 'iregread', 'alu0', 'alu1_csr', 'fp', 'mem_iq', 'lsu']
glue = {}
for n in CORE:
    for k, f in assign[n].items(): glue[k] = glue.get(k, 0) + f * share[n]
sg = sum(glue.values()); glue_tiles_only = {k: v / sg for k, v in glue.items()}
assign['core_glue'] = glue_tiles_only
out = dict(die_mm2=DIE_MM2, n=N, tile_mm2=DIE_MM2 / N / N, tile_um=1e3 * (DIE_MM2 ** 0.5) / N,
           blocks={n: dict(leak_mW=v, share=share[n], mm2=share[n] * DIE_MM2, tiles=assign[n]) for n, v in blocks})
json.dump(out, open('/home/ahilnandan/Documents/PMIC/pdn_work/lacpo/floorplan.json', 'w'), indent=1)
grid = [['' for _ in range(N)] for _ in range(N)]
for n in order:
    for k, f in assign[n].items():
        c, r = map(int, k.split('_')); grid[r][c] = (grid[r][c] + '/' if grid[r][c] else '') + n[:6]
print('tile %.0f um, %.4f mm^2' % (out['tile_um'], out['tile_mm2']))
for n, v in blocks: print('  %-12s leak %6.2f mW  area %5.1f %%  %.3f mm^2  tiles %4.1f' % (n, v, 100 * share[n], share[n] * DIE_MM2, share[n] * 64))
print('\nfloorplan (rows top->bottom):')
for r in range(N): print('  ' + ' | '.join('%-13s' % grid[r][c][:13] for c in range(N)))
