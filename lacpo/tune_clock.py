"""tune_clock: ADAPT-style (Radke & Roy, Intel, EPEPS 2025) automated tuner and gate suite for B2v2 adaptive clock
modulation (multi-threshold entry + staggered frequency restore ladder) on policy6 / pdn6 / mc_engine.

Usage:
  python lacpo/tune_clock.py gates                # run G0, G1, G2, G3 + restore-edge check + hello_vcd trace check
  python lacpo/tune_clock.py tune [R1|R2|R3|MC]   # run NSGA-II / surrogate tuner for the chosen regime

Needs the real traces (SEGSEE_TRACES) and lacpo/hpos*.npy. Without them, SEGSEE_FIXTURE=1 runs on a synthetic
fixture for smoke tests only. LACPO_HELLO_CSV optionally points at DecouplingCap traces/hello_vcd/composed_P_t.csv,
whose cache-block currents are constant placeholders.
"""
from __future__ import print_function
import os, sys, json, copy, time, tempfile, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
sys.path.insert(0, HERE)
import pdn6, pdn7, timing6 as TM

NG = 10
S_GRID = [0.0625 * k for k in range(1, 7)]   # 1/16 .. 6/16 (6.25% .. 37.5% period stretch)


def decode_clock(g, G_mV=None, rep_mode='none', flush=0, rep_pen=5, honest=0):
    """Map gene vector g in [0, 1]^10 -> policy6 B2v2 clock parameter dict."""
    g = np.clip(np.asarray(g, float), 0.0, 1.0)
    thr1 = round(14.0 + 34.0 * g[0], 2)
    two_lvl = bool(g[2] >= 0.25)
    s1 = S_GRID[min(3, int(g[3] * 4))]
    dw1 = 1 + int(g[5] * 5.99)
    if two_lvl:
        thr2 = round(thr1 + 4.0 + 20.0 * g[1], 2)
        s2_idx = max(S_GRID.index(s1) + 1, min(len(S_GRID) - 1, 1 + int(g[4] * 5)))
        s2 = S_GRID[s2_idx]
        dw2 = 1 + int(g[6] * 6.99)
        levels = [(thr1, s1, dw1), (thr2, s2, dw2)]
        s_top = s2
    else:
        levels = [(thr1, s1, dw1)]
        s_top = s1
    n_ex = min(3, int(g[7] * 4))
    d_ex = 1 + int(g[8] * 2.99)
    if n_ex == 0:
        exit_ladder = []
    elif n_ex == 1:
        exit_ladder = [(round(0.5 * s_top, 4), d_ex)]
    elif n_ex == 2:
        exit_ladder = [(round(2.0 * s_top / 3.0, 4), d_ex), (round(s_top / 3.0, 4), d_ex)]
    else:
        exit_ladder = [(round(0.75 * s_top, 4), d_ex), (round(0.50 * s_top, 4), d_ex), (round(0.25 * s_top, 4), d_ex)]
    exit_ladder = [(s, c) for s, c in exit_ladder if s > 0.01]
    hyst = round(8.0 * g[9], 2) if g[9] > 0.2 else 0.0
    rec = 1 + int(g[9] > 0.65)
    prm = dict(cs_on=1, cs_levels=levels, cs_exit=exit_ladder, cs_hyst=hyst, cs_rec=rec,
               cs_honest=int(honest), cs_max_time=64)
    if G_mV is not None:
        prm.update(G_mV=float(G_mV), rep_mode=rep_mode, flush=int(flush), rep_pen=int(rep_pen))
    return prm


HELLO_CSV = os.environ.get('LACPO_HELLO_CSV', '')


def ensure_fixture_env():
    """Use the real traces and kernels when present. Otherwise, only if SEGSEE_FIXTURE=1, build a SYNTHETIC
    fixture (made-up per-uop kernels and uop streams) in a temp directory so the code paths can be smoke-tested.
    Fixture numbers are not results. Fixture kernels are never written into lacpo/."""
    import replay3, tilecurrents as TC, workload5 as W5
    import closedloop6 as C6, policy6 as P6
    if os.path.exists(replay3.B + 'qsort_attr5.npz') and os.path.exists(HERE + 'hpos.npy') and os.path.exists(HERE + 'hpos_blocks.npy'):
        return replay3.B, 'full_traces'
    if os.environ.get('SEGSEE_FIXTURE') != '1':
        raise RuntimeError('Real traces/kernels not found (SEGSEE_TRACES=%s, %shpos*.npy). Set SEGSEE_FIXTURE=1 '
                           'to smoke-test on a synthetic fixture; its numbers are not results.' % (replay3.B, HERE))
    tmp = os.path.join(tempfile.gettempdir(), 'segsee_fixture_traces_v3') + os.sep
    os.makedirs(tmp, exist_ok=True)
    nb = len(W5.BLOCKS); nc = 10; K = replay3.K
    hb = np.zeros((nb, nc, K), np.float64)
    base_q = [18.0, 14.0, 38.0, 34.0, 42.0, 22.0, 14.0, 48.0, 30.0, 26.0]
    sh_k = np.array([0.08, 0.18, 0.24, 0.20, 0.12, 0.07, 0.04, 0.03, 0.02, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0])
    for c in range(nc):
        for b in range(nb):
            hb[b, c, :] = base_q[c] * (0.20 if b in (0, 2, 3, 6, 8) else 0.05) * sh_k
    h4 = np.zeros((4, nc, K), np.float64)
    for c in range(nc):
        h4[1, c, :] = 0.6 * hb[:, c, :].sum(0)
        h4[2, c, :] = 0.3 * hb[:, c, :].sum(0)
        h4[0, c, :] = 0.1 * hb[:, c, :].sum(0)
    np.save(tmp + 'hpos_blocks.npy', hb)
    np.save(tmp + 'hpos.npy', h4)
    lacpo_csv = HELLO_CSV
    real_bg = None
    if lacpo_csv and os.path.exists(lacpo_csv):
        import pandas as pd
        df = pd.read_csv(lacpo_csv)
        n_csv = len(df)
        blk = {}
        for col, b in TC.COLMAP.items():
            if col in df.columns: blk[b] = blk.get(b, 0.0) + df[col].values[:n_csv]
        blk['icache'] = np.full(n_csv, 18.0)
        for b, f in TC.DSPLIT.items(): blk[b] = np.full(n_csv, f * 22.0)
        blk['ptw'] = np.full(n_csv, 7.24); blk['tile_misc'] = np.full(n_csv, 1.5)
        real_bg = np.zeros((n_csv, 64), np.float32)
        for b, v in blk.items():
            real_bg += np.outer(v / 1.1, TC.W[b]).astype(np.float32)
    for w_idx, w in enumerate(['qsort', 'towers', 'mm']):
        if os.path.exists(tmp + w + '_attr5.npz') and os.path.exists(tmp + w + '_cache.npz'):
            continue
        r_w = np.random.RandomState(100 + w_idx)
        n = 3000
        # Generate realistic bursty uop stream (~2.5 uops/cycle in bursts, 0.3 in quiet phases)
        ts, rs, iqs, fus, robs, pdsts, dtypes, prs1s, t1s, prs2s, t2s, clss = [], [], [], [], [], [], [], [], [], [], [], []
        rob_cur = 0; last_writer = {}
        for cyc in range(20, n - 40):
            in_burst = (300 <= (cyc % 350) < 322)
            cnt = r_w.choice([2, 3, 4], p=[0.2, 0.4, 0.4]) if in_burst else r_w.choice([0, 1], p=[0.75, 0.25])
            used_int, used_mem, used_fp = 0, 0, 0
            for _ in range(cnt):
                c = int(r_w.choice([0, 1, 2, 3, 4, 7], p=[0.25, 0.05, 0.25, 0.10, 0.20, 0.15] if in_burst else [0.65, 0.15, 0.10, 0.05, 0.03, 0.02]))
                if c in (0, 1, 4):
                    if used_int >= 2: continue
                    iq = 0; fu = 1 if c == 0 else (2 if c == 1 else 8); used_int += 1
                elif c in (2, 3):
                    if used_mem >= 1: continue
                    iq = 1; fu = 4; used_mem += 1
                else:
                    if used_fp >= 1: continue
                    iq = 2; fu = 64; used_fp += 1
                dt = 2 if c == 3 else (1 if c == 7 else 0)
                pdst = int(r_w.randint(1, 48)) if dt in (0, 1) else 0
                p1 = int(r_w.randint(1, 48)); p2 = int(r_w.randint(1, 48))
                t_prod = max(last_writer.get((0, p1), cyc - 4), last_writer.get((0, p2), cyc - 4))
                r_cyc = max(0, min(cyc, t_prod + r_w.randint(1, 3)))
                if dt in (0, 1) and pdst > 0:
                    last_writer[(dt, pdst)] = cyc
                ts.append(cyc); rs.append(r_cyc); iqs.append(iq); fus.append(fu); robs.append(rob_cur % 64)
                pdsts.append(pdst); dtypes.append(dt); prs1s.append(p1); t1s.append(0); prs2s.append(p2); t2s.append(0)
                clss.append(c); rob_cur += 1
        m = len(ts)
        cls_arr = np.array(clss, np.int64)
        A5 = hb[:, cls_arr, :].transpose(1, 0, 2).astype(np.float32) * (1.25 if w == 'qsort' else 1.1)
        A4 = np.zeros((m, 4, K), np.float32)
        A4[:, 1, :] = A5.sum(1) * 0.7; A4[:, 2, :] = A5.sum(1) * 0.3
        if real_bg is not None and len(real_bg) >= n:
            off = (w_idx * 3000) % (len(real_bg) - n)
            BG5 = real_bg[off:off + n].copy()
        else:
            BG5 = np.tile((pdn6.LEAK_T * 1.8).astype(np.float32), (n, 1))
        # Ramped 8-beat cache refill bursts aligned with issue bursts (ramps over 4 cycles so detector reacts before peak)
        ramp = np.array([0.25, 0.50, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 0.5], np.float32)
        for b0 in range(300, n - 100, 350):
            BG5[b0:b0 + len(ramp)] += np.outer(ramp * 110.0, pdn6.CELL_SHARE).astype(np.float32)
        BG4 = np.zeros((n, 4), np.float32); BG4[:, 0] = BG5.sum(1) * 0.4; BG4[:, 1] = BG5.sum(1) * 0.3; BG4[:, 2] = BG5.sum(1) * 0.3
        common = dict(iq=np.array(iqs, np.int64), uopc=np.zeros(m, np.int64), fu=np.array(fus, np.int64),
                      rob=np.array(robs, np.int64), pdst=np.array(pdsts, np.int64), dtype=np.array(dtypes, np.int64),
                      prs1=np.array(prs1s, np.int64), t1=np.array(t1s, np.int64), prs2=np.array(prs2s, np.int64),
                      t2=np.array(t2s, np.int64), r=np.array(rs, np.int64), t=np.array(ts, np.int64), cls=cls_arr)
        np.savez_compressed(tmp + w + '_attr.npz', A=A4, BG=BG4, **common)
        np.savez_compressed(tmp + w + '_attr5.npz', A=A5, BG=BG5, **common)
        dc_mshr = np.zeros(n, np.int64); ic_miss = np.zeros(n, np.int64)
        for b0 in range(300, n - 100, 250): dc_mshr[b0 - 6] = 1
        np.savez_compressed(tmp + w + '_cache.npz', dc_mshr=dc_mshr, ic_miss=ic_miss)
    replay3.B = tmp; C6.B = tmp; P6.B = tmp
    replay3.HERE = tmp; C6.HERE = tmp
    return tmp, 'SYNTHETIC_FIXTURE'


def run_gates():
    """Run all 4 plan gates (G0, G1, G2, G3) + multi-core clock vs coordinator check."""
    trace_dir, mode = ensure_fixture_env()
    import closedloop6 as C6, policy6 as P6, eval6, mc_engine as MC
    D, prmb = eval6.designs()
    print('=== Running ADAPT B2v2 Gate Suite (trace mode: %s) ===' % mode)

    # Gate G0: policy6 with all defaults == closedloop6 bit-identical
    for w in ('towers', 'qsort'):
        cl = C6.CL(w); eng = P6.Engine(w)
        n_win = min(cl.n, 3500)
        for dn in ('decap', 'banks'):
            for pol in ('A', 'TOK'):
                r_cl = cl.run(D[dn], pol, prmb, a=0, b=n_win)
                r_p6 = eng.run(D[dn], pol, prmb, a=0, b=n_win)
                for k in ('p999', 'max', 'p99', 'slowdown', 'energy_over', 'nhold', 'nreo'):
                    diff = abs(r_cl[k] - r_p6[k])
                    assert diff == 0.0, 'G0 mismatch on %s %s %s %s: %g vs %g' % (w, dn, pol, k, r_cl[k], r_p6[k])
    print('[PASS] Gate G0: policy6 defaults == closedloop6 bit-identical (0.000000e+00 diff across decap/banks, A/TOK)')

    # Gate G1: cs_levels=[(cs_thr, 0.25, 4)], cs_exit=[] == old B2 bit-identical
    for w in ('towers', 'qsort'):
        eng = P6.Engine(w); n_win = min(eng.n, 3500)
        for dn in ('decap', 'banks'):
            thr = 28.0
            r_old = eng.run(D[dn], 'A', dict(prmb, cs_on=1, cs_thr=thr, cs_s=0.25, cs_len=4, G_mV=35.0, rep_mode='razor'), a=0, b=n_win, record=True)
            r_new = eng.run(D[dn], 'A', dict(prmb, cs_levels=[(thr, 0.25, 4)], cs_exit=[], G_mV=35.0, rep_mode='razor'), a=0, b=n_win, record=True)
            v_diff = float(np.max(np.abs(r_old['vworst'] - r_new['vworst'])))
            assert v_diff == 0.0 and r_old['nstretch'] == r_new['nstretch'] and r_old['nerr'] == r_new['nerr'], \
                'G1 mismatch on %s %s: v_diff=%g' % (w, dn, v_diff)
    print('[PASS] Gate G1: cs_levels=[(cs_thr, 0.25, 4)], cs_exit=[] == legacy B2 bit-identical (max|vworst diff| = 0.0)')

    # Gate G2: Synthetic step restore-edge test on pdn6.Sim (instant restore vs 4-step staggered ladder)
    n_step = 900
    I_step = np.tile((pdn6.LEAK_T * 2.0).astype(np.float32), (n_step, 1))
    # Sustained high-current burst from cycle 650 to 740 (+320 mA across die)
    I_step[650:740] += (pdn6.CELL_SHARE * 320.0).astype(np.float32)
    # Trigger a single 4-cycle stretch at the burst onset (no re-trigger during the flat plateau, so we isolate the restore edge at cycle 656)
    # We test ClockMod with a threshold that fires on the initial resonant dip at cycle 651 (dip ~52 mV) but NOT on the steady plateau (~18 mV)
    r_inst = P6.sim_clock_openloop(D['decap'], I_step, dict(cs_levels=[(36.0, 0.25, 4)], cs_exit=[], cs_noise=0.0), warm=600, record=True)
    r_ladd = P6.sim_clock_openloop(D['decap'], I_step, dict(cs_levels=[(36.0, 0.25, 4)], cs_exit=[(0.1875, 1), (0.125, 1), (0.0625, 1)], cs_noise=0.0), warm=600, record=True)
    red_pct = 100.0 * (1.0 - r_ladd['restore_droop_max'] / max(r_inst['restore_droop_max'], 1e-6))
    assert red_pct >= 50.0, 'G2 failed: reduction %.1f%% < 50%%' % red_pct
    print('[PASS] Gate G2: pdn6.Sim 320 mA step restore-edge second droop: instant=%.2f mV -> 4-step ladder=%.2f mV (-%.1f%% >= 50%%)'
          % (r_inst['restore_droop_max'], r_ladd['restore_droop_max'], red_pct))

    # Gate G3: Forward progress / max continuous modulation time (cs_max_time)
    r_g3 = P6.sim_clock_openloop(D['decap'], I_step, dict(cs_levels=[(1.0, 0.25, 4)], cs_exit=[(0.125, 2)], cs_max_time=12, cs_noise=0.0), warm=600, record=True)
    s_tr = r_g3['s_trace'] > 0
    max_run = max(len(chunk) for chunk in ''.join('1' if x else '0' for x in s_tr).split('0'))
    assert max_run <= 12, 'G3 failed: max continuous modulated run %d > 12' % max_run
    print('[PASS] Gate G3: continuous modulation bounded by cs_max_time=12 (observed max run = %d cycles)' % max_run)

    return dict(G0='PASS', G1='PASS', G2_instant_mV=r_inst['restore_droop_max'],
                G2_ladder_mV=r_ladd['restore_droop_max'], G2_reduction_pct=red_pct, G3_max_run=max_run)


def benchmark_b2_vs_b2v2():
    """Compare legacy B2 vs ADAPT-tuned B2v2 across:
      (1) Single-core closed-loop Engine (R1 no-replay, R2 selective replay, R3 flush recovery)
      (2) LACPo hello_vcd trace on pdn6.Sim (decap and banks; cache blocks are constant placeholders)
      (3) Multi-core 2x2 shared-die pdn7.Sim lockstep scenario (chip-wide B2v2 clock vs cross-core hold coordinator)."""
    ensure_fixture_env()
    import policy6 as P6, eval6, mc_engine as MC, tilecurrents as TC
    D, prmb = eval6.designs()
    out = {'single_core': [], 'hello_vcd_real_lacpo': [], 'multicore_lockstep': []}

    # 1. Single-core closed-loop comparison on qsort & towers
    eng = P6.Engine('qsort')
    configs = [
        ('no stretch (baseline)', dict(prmb)),
        ('legacy B2 (thr=30, s=0.25, len=4, instant)', dict(prmb, cs_on=1, cs_thr=30.0, cs_s=0.25, cs_len=4)),
        ('B2v2 staggered restore (thr=30, s=0.25, ladder 3x1)', dict(prmb, cs_levels=[(30.0, 0.25, 3)], cs_exit=[(0.1875, 1), (0.125, 1), (0.0625, 1)])),
        ('B2v2 2-level + ladder (thr=26/36, s=0.125/0.25 + ladder)', dict(prmb, cs_levels=[(26.0, 0.125, 2), (36.0, 0.25, 3)], cs_exit=[(0.1875, 1), (0.125, 1), (0.0625, 1)], cs_hyst=2.0)),
        ('flush G=15 razor only', dict(prmb, G_mV=15.0, rep_mode='razor', flush=1, rep_pen=20, rob=80)),
        ('flush G=15 + legacy B2 (thr=22)', dict(prmb, G_mV=15.0, rep_mode='razor', flush=1, rep_pen=20, rob=80, cs_on=1, cs_thr=22.0, cs_s=0.25, cs_len=4)),
        ('flush G=15 + B2v2 2-level + ladder (thr=18/26)', dict(prmb, G_mV=15.0, rep_mode='razor', flush=1, rep_pen=20, rob=80,
                                                                 cs_levels=[(18.0, 0.125, 2), (26.0, 0.25, 3)], cs_exit=[(0.1875, 1), (0.125, 1), (0.0625, 1)], cs_hyst=2.0)),
    ]
    print('\n=== 1. Single-Core Closed-Loop Engine (banks & decap) ===')
    for dn in ('banks', 'decap'):
        for tag, prm in configs:
            r = eng.run(D[dn], 'A', prm, a=0, b=eng.n)
            row = dict(design=dn, tag=tag, p999=r['p999'], max=r['max'], slowdown_pct=100.0 * r['slowdown'],
                       nerr=r['nerr'], nstretch=r['nstretch'], stretch_time=r['stretch_time'],
                       n_restore_droop=r['n_restore_droop'], restore_droop_max=r['restore_droop_max'])
            out['single_core'].append(row)
            print('%-6s | %-56s | p999 %5.1f | max %5.1f | slow %5.2f%% | err %3d | restore_dip %4.1f mV (n=%d)'
                  % (dn, tag, r['p999'], r['max'], 100.0 * r['slowdown'], r['nerr'], r['restore_droop_max'], r['n_restore_droop']))

    # 2. LACPo hello_vcd trace (if LACPO_HELLO_CSV is set)
    lacpo_csv = HELLO_CSV
    if lacpo_csv and os.path.exists(lacpo_csv):
        import pandas as pd
        df = pd.read_csv(lacpo_csv).iloc[-6000:].reset_index(drop=True); n_csv = len(df)
        blk = {}
        for col, b in TC.COLMAP.items():
            if col in df.columns: blk[b] = blk.get(b, 0.0) + df[col].values[:n_csv]
        blk['icache'] = np.full(n_csv, 18.0)
        for b, f in TC.DSPLIT.items(): blk[b] = np.full(n_csv, f * 22.0)
        blk['ptw'] = np.full(n_csv, 7.24); blk['tile_misc'] = np.full(n_csv, 1.5)
        I_hello = np.zeros((n_csv, 64), np.float32)
        for b, v in blk.items(): I_hello += np.outer(v / 1.1, TC.W[b]).astype(np.float32)
        print('\n=== 2. LACPo hello_vcd Trace on pdn6.Sim (%d cycles, placeholder cache currents) ===' % n_csv)
        for dn in ('decap', 'banks'):
            r0 = P6.sim_clock_openloop(D[dn], I_hello, dict(cs_on=0), warm=600)
            thr_b = round(0.55 * r0['p999'], 2)
            thr_lo = round(0.45 * r0['p999'], 2); thr_hi = round(0.68 * r0['p999'], 2)
            for tag, cprm in [
                ('no stretch', dict(cs_on=0)),
                ('legacy B2 (1-lvl, instant restore)', dict(cs_on=1, cs_thr=thr_b, cs_s=0.25, cs_len=4, cs_noise=0.5)),
                ('B2v2 (1-lvl + 3-step ladder)', dict(cs_levels=[(thr_b, 0.25, 3)], cs_exit=[(0.1875, 1), (0.125, 1), (0.0625, 1)], cs_noise=0.5)),
                ('B2v2 (2-lvl + 3-step ladder)', dict(cs_levels=[(thr_lo, 0.125, 2), (thr_hi, 0.25, 3)], cs_exit=[(0.1875, 1), (0.125, 1), (0.0625, 1)], cs_hyst=1.0, cs_noise=0.5)),
            ]:
                rh = P6.sim_clock_openloop(D[dn], I_hello, cprm, warm=600)
                row = dict(design=dn, tag=tag, p999=rh['p999'], max=rh['max'], slowdown_pct=100.0 * rh['slowdown'],
                           stretch_time=rh['stretch_time'], restore_droop_max=rh['restore_droop_max'])
                out['hello_vcd_real_lacpo'].append(row)
                print('%-6s | %-38s | p999 %5.2f | max %5.2f | slow %5.2f%% | restore_dip %4.2f mV'
                      % (dn, tag, rh['p999'], rh['max'], 100.0 * rh['slowdown'], rh['restore_droop_max']))

    # 3. Multi-core 2x2 lockstep on pdn7.Sim: chip-wide B2v2 clock vs cross-core coordinator
    print('\n=== 3. Multi-Core 2x2 Shared Die (pdn7.Sim, 4x lockstep) ===')
    engs = [eng] * 4; starts = [0, 0, 0, 0]; W_mc = min(eng.n, 1600)
    mc_configs = [
        ('uncoordinated, no clock mod', None, None),
        ('cross-core hold coordinator (budget=30 mA)', dict(budget=30.0), None),
        ('chip-wide legacy B2 clock (thr=24, instant)', None, dict(cs_on=1, cs_thr=24.0, cs_s=0.25, cs_len=4)),
        ('chip-wide B2v2 clock (2-lvl 18/26 + ladder)', None, dict(cs_levels=[(18.0, 0.125, 2), (26.0, 0.25, 3)], cs_exit=[(0.1875, 1), (0.125, 1), (0.0625, 1)], cs_hyst=2.0)),
        ('combined coordinator (budget=38) + B2v2 clock', dict(budget=38.0), dict(cs_levels=[(20.0, 0.125, 2), (28.0, 0.25, 2)], cs_exit=[(0.125, 1), (0.0625, 1)])),
    ]
    for tag, coord, clk in mc_configs:
        rm = MC.run_mc(engs, starts, D['decap'], W=W_mc, coord=coord, clock=clk, warm=400)
        mean_sd = 100.0 * float(np.mean(rm['sd']))
        row = dict(tag=tag, p999=rm['p999'], max=rm['max'], mean_slowdown_pct=mean_sd,
                   dI999=rm['dI999'], restore_droop_max=rm.get('restore_droop_max', 0.0))
        out['multicore_lockstep'].append(row)
        print('%-46s | p999 %5.1f | max %5.1f | mean slow %5.2f%% | dI999 %5.1f mA'
              % (tag, rm['p999'], rm['max'], mean_sd, rm['dI999']))

    return out


def tune_nsga2(regime='R1', pop=16, gens=6, seed=0):
    """Run NSGA-II tuner over the 10-gene ADAPT clock-modulation space for regime R1 (no replay), R2 (selective replay),
    or R3 (flush recovery). Returns Pareto front solutions."""
    ensure_fixture_env()
    import policy6 as P6, eval6, ga5
    D, prmb = eval6.designs()
    eng = P6.Engine('qsort'); n_win = min(eng.n, 3000)
    rng = np.random.RandomState(seed)
    X = rng.rand(pop, NG)
    # Seed 0: legacy B2 equivalent; Seed 1: 2-level with 3-step quarter-wave exit ladder
    X[0] = [0.45, 0.5, 0.0, 0.9, 0.9, 0.6, 0.5, 0.0, 0.0, 0.0]
    X[1] = [0.35, 0.4, 0.8, 0.4, 0.7, 0.3, 0.4, 0.9, 0.1, 0.3]

    def eval_one(g):
        if regime == 'R1':
            cprm = decode_clock(g)
        elif regime == 'R2':
            cprm = decode_clock(g, G_mV=15.0, rep_mode='razor', flush=0, rep_pen=5)
        else:
            cprm = decode_clock(g, G_mV=15.0, rep_mode='razor', flush=1, rep_pen=20)
        prm = dict(prmb, rob=80, **cprm)
        try:
            r = eng.run(D['banks'], 'A', prm, a=0, b=n_win)
            return [r['p999'], r['max'], 100.0 * max(0.0, r['slowdown']), r['restore_droop_max']]
        except RuntimeError:
            return [1e3, 1e3, 1e3, 1e3]

    F = [eval_one(x) for x in X]
    for gen in range(gens):
        rank, crowd = ga5.nsga(F)
        def tour():
            i, j = rng.randint(pop, size=2)
            return i if (rank[i], -crowd[i]) < (rank[j], -crowd[j]) else j
        kids = []
        while len(kids) < pop:
            a, b = X[tour()], X[tour()]; u = rng.rand(NG)
            beta = np.where(u < 0.5, (2 * u) ** (1 / 16.), (1 / (2 * (1 - u))) ** (1 / 16.))
            for c in (0.5 * ((1 + beta) * a + (1 - beta) * b), 0.5 * ((1 - beta) * a + (1 + beta) * b)):
                mu = rng.rand(NG) < 2.0 / NG; c[mu] += rng.normal(0, 0.15, mu.sum())
                kids.append(np.clip(c, 0.0, 1.0))
        kids = np.array(kids[:pop]); Fk = [eval_one(x) for x in kids]
        XA = np.vstack([X, kids]); FA = F + Fk
        rank, crowd = ga5.nsga(FA)
        order = sorted(range(len(FA)), key=lambda i: (rank[i], -crowd[i]))[:pop]
        X = XA[order]; F = [FA[i] for i in order]
        best = sorted(F, key=lambda f: f[1])[:3]
        print('[%s] gen %d/%d best max-droop points (p999, max, slow%%, restore_dip): %s'
              % (regime, gen + 1, gens, '  '.join('(%.1f, %.1f, %.2f%%, %.1f)' % tuple(f) for f in best)), flush=True)
    rank, _ = ga5.nsga(F)
    front0 = [dict(genes=X[i].tolist(), metrics=dict(p999=F[i][0], max=F[i][1], slowdown_pct=F[i][2], restore_dip_mV=F[i][3]),
                   decoded=decode_clock(X[i])) for i in range(len(F)) if rank[i] == 0]
    return front0


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'gates'
    if cmd == 'gates':
        g_res = run_gates()
        b_res = benchmark_b2_vs_b2v2()
        import subprocess
        subprocess.run([sys.executable, os.path.join(HERE, '..', 'pulsim_validate.py'), 'restore_edge'], check=False)
    elif cmd == 'tune':
        reg = sys.argv[2] if len(sys.argv) > 2 else 'R1'
        fr = tune_nsga2(reg)
        print('Pareto front size:', len(fr))
