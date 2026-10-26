"""Phase 5 evaluation on all 10 FULL traces (policy6 engine). Sets:
  A6  circuit-state-aware replay vs baselines at guardband G (timing errors on): razor, razor + error-rate voltage
      scaling, clock stretch + razor, state replay alone, our scheduler + razor, ours (scheduler + state replay + error
      feedback). Designs: GA5 0.35 lowest-max (banks) and 0.35 all-decap.
  B8  scheduler-v2 ablation, droop only (no timing errors): each term alone, all terms, separate-rail banks, stretch.
Weights for the new terms are first guesses (not tuned); GA6 tunes them. Usage: python eval6.py <A6|A6b|A6c|B8|FULL> [workers]."""
import os, sys, json, copy, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
BENCHES = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
def designs():
    import ga5, val5
    fr, _ = val5.front('0.35'); gm = [g for nm, g, _ in fr if nm == 'min_max'][0]
    Db, polb, prmb, _ = ga5.decode(gm, 0.35); Dd, _, _, _ = ga5.decode([0.0] * 16, 0.35)
    Ds = copy.deepcopy(Db)
    for t in Ds['bank']: Ds['bank'][t]['src'] = 1
    return dict(banks=Db, decap=Dd, banks_sep=Ds), prmb
OURS = dict(rep_mode='state', fb_err=2.0, err_arm=1, cls_w=1.0, look_arm=1)
def configs(which):
    J = []
    if which == 'A6':
        for dn in ('banks', 'decap'):
            for G in (50, 40, 30):
                J += [(dn, 'razor', 'A', dict(G_mV=G, rep_mode='razor')),
                      (dn, 'razor+EVS', 'A', dict(G_mV=G, rep_mode='razor', evs_gain=0.01, evs_target=0.1)),
                      (dn, 'stretch+razor', 'A', dict(G_mV=G, rep_mode='razor', cs_on=1, cs_thr=G - 10.0)),
                      (dn, 'state replay', 'A', dict(G_mV=G, rep_mode='state')),
                      (dn, 'sched+razor', 'TOK', dict(G_mV=G, rep_mode='razor')),
                      (dn, 'OURS', 'TOK', dict(OURS, G_mV=G))]
    elif which == 'A6b':                                   # where Razor should break: low guardband, full-flush penalty
        OURS2 = dict(rep_mode='state', fb_err=2.0, err_arm=1, cls_w=1.0)    # no refill-lookahead arming (starved banks in A6/B8)
        for dn in ('banks', 'decap'):
            for G in (20, 10):
                for pen in (5, 20):
                    J += [(dn, 'razor p%d' % pen, 'A', dict(G_mV=G, rep_mode='razor', rep_pen=pen)),
                          (dn, 'state p%d' % pen, 'A', dict(G_mV=G, rep_mode='state', rep_pen=pen)),
                          (dn, 'OURS2 p%d' % pen, 'TOK', dict(OURS2, G_mV=G, rep_pen=pen)),
                          (dn, 'razor+EVS p%d' % pen, 'A', dict(G_mV=G, rep_mode='razor', rep_pen=pen, evs_gain=0.01, evs_target=0.1))]
    elif which == 'A6c':                                   # fixed engine, ROB on, learned forecaster (held-out weights)
        for dn in ('banks', 'decap'):
            F7 = dict(fc_mode='fc7', fc7_design=dn, fc7_lag0=0, fc7_tag='train7', rob=80)
            for G in (30, 20, 10):
                J += [(dn, 'sel razor', 'A', dict(F7, G_mV=G, rep_mode='razor', rep_pen=5)),
                      (dn, 'sel state', 'A', dict(F7, G_mV=G, rep_mode='state', rep_pen=5)),
                      (dn, 'flush razor', 'A', dict(F7, G_mV=G, rep_mode='razor', rep_pen=20, flush=1)),
                      (dn, 'flush state', 'A', dict(F7, G_mV=G, rep_mode='state', rep_pen=20, flush=1)),
                      (dn, 'flush avoid5', 'A', dict(F7, G_mV=G, rep_mode='razor', rep_pen=20, flush=1, avoid_vm=5.0)),
                      (dn, 'flush avoid5+state', 'A', dict(F7, G_mV=G, rep_mode='state', rep_pen=20, flush=1, avoid_vm=5.0)),
                      (dn, 'flush razor+EVS', 'A', dict(F7, G_mV=G, rep_mode='razor', rep_pen=20, flush=1, evs_gain=0.01, evs_target=0.1))]
    elif which == 'FULL':                                  # final single-core stack (power-aware scheduling + replay + banks)
        TH = dict(banks=dict(th=19.0, th_arm=28.0), decap=dict(th=24.0, th_arm=30.0))       # steer7 calibration (fc7 scale)
        for dn in ('decap', 'banks', 'banks_sep'):
            kd = 'decap' if dn == 'decap' else 'banks'
            F7 = dict(fc_mode='fc7', fc7_design=kd, fc7_lag0=0, fc7_tag='train7', rob=80)
            S1 = dict(F7, kv=0.0, hold=1, maxd=2, **TH[kd])
            J += [(dn, 'A', 'A', dict(F7)), (dn, 'S1 sched', 'TOK', dict(S1)), (dn, 'S2 sched+lock', 'TOK', dict(S1, lock_k=2))]
            for G in (30, 20, 10):
                J += [(dn, 'A+replay G%d' % G, 'A', dict(F7, G_mV=G, rep_mode='razor', rep_pen=5)),
                      (dn, 'S1+replay G%d' % G, 'TOK', dict(S1, G_mV=G, rep_mode='razor', rep_pen=5))]
    else:
        T = 'TOK'
        J = [('banks', 'A (no policy)', 'A', {}), ('banks', 'TOK base', T, {}),
             ('banks', '+edge 15 mA', T, dict(edge_dI=15.0)), ('banks', '+edge 25 mA', T, dict(edge_dI=25.0)),
             ('banks', '+window 3.6 MHz', T, dict(win_dI=10.0, win_k=0.3)), ('banks', '+charge debt', T, dict(debt_k=0.02)),
             ('banks', '+refill lookahead', T, dict(look_k=4.0, look_arm=1)), ('banks', '+77 MHz tone', T, dict(tone_thr=3.0, tone_k=1.0)),
             ('banks', '+bank readiness', T, dict(ready_min=0.5, ready_k=8.0)),
             ('banks', '+all terms', T, dict(edge_dI=25.0, win_dI=10.0, win_k=0.3, debt_k=0.02, look_k=4.0, look_arm=1, tone_thr=3.0, tone_k=1.0, ready_min=0.5, ready_k=8.0)),
             ('banks_sep', 'TOK, separate-rail banks', T, {}), ('banks', 'A + stretch thr 40', 'A', dict(cs_on=1, cs_thr=40.0))]
    return J
_eng = {}
def job(a):
    which, dn, tag, pol, extra, bench = a
    import policy6 as P6, analysis6 as AN
    if bench not in _eng: _eng.clear(); _eng[bench] = P6.Engine(bench)
    D, prmb = designs(); prm = dict(prmb, **extra)
    try: r = _eng[bench].run(D[dn], pol, prm, record=True)
    except RuntimeError as ex:                            # engine fell > 20k cycles behind: recovery collapse
        return dict(bench=bench, design=dn, tag=tag, G=extra.get('G_mV'), collapse=True, err=str(ex))
    bd = AN.bands(r['vworst_pdn']); jt = AN.jitter(r['vworst_pdn'])
    keep = {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}
    keep.update(bench=bench, design=dn, tag=tag, G=extra.get('G_mV'), bands=bd, jitter=jt)
    return keep
if __name__ == '__main__':
    import multiprocessing as mp
    which = sys.argv[1]; workers = min(6, int(sys.argv[2]) if len(sys.argv) > 2 else 6)
    jobs = [(which, dn, tag, pol, extra, b) for b in BENCHES for (dn, tag, pol, extra) in configs(which)]
    print('%s: %d jobs' % (which, len(jobs)), flush=True)
    res = []; part = open(HERE + 'eval6_%s.partial.jsonl' % which, 'w')
    with mp.get_context('spawn').Pool(workers) as P:
        for r in P.imap_unordered(job, jobs, chunksize=1):
            res.append(r); part.write(json.dumps(r, default=float) + '\n'); part.flush()
    json.dump(res, open(HERE + 'eval6_%s.json' % which, 'w'), indent=1, default=float)
    print('DONE', which, flush=True)
