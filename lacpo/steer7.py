"""steer7: temporal issue steering (TOK reorder / holds) with the learned fc7 forecaster vs the old forecast.
policy6.py line 83 has code swallowed by a comment (as of 2026-10-24); it is patched IN MEMORY here (file untouched).
usage: steer7.py cal | run"""
import sys, os, json, types, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'; sys.path.insert(0, HERE)
def load_p6():
    if 'policy6' in sys.modules: return sys.modules['policy6']
    s = open(HERE + 'policy6.py').read()
    bad = "# room for up to ~200 % slowdown before 'overflow' ring ="
    if bad in s: s = s.replace(bad, "# room for up to ~200 % slowdown before 'overflow'\n        ring =")
    M = types.ModuleType('policy6'); M.__file__ = HERE + 'policy6.py'; sys.modules['policy6'] = M
    exec(compile(s, M.__file__, 'exec'), M.__dict__); return M
BENCHES = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
HELD = {'qsort', 'mm', 'towers'}
_eng = {}
def job(a):
    dn, tag, pol, extra, bench = a
    P6 = load_p6(); import eval6
    if bench not in _eng: _eng.clear(); _eng[bench] = P6.Engine(bench)
    D, prmb = eval6.designs(); prm = dict(prmb, rob=80); prm.update(extra)
    if prm.get('fc_mode') == 'fc7': prm.update(fc7_design=dn, fc7_lag0=0, fc7_tag='train7')
    t0 = time.time()
    try: r = _eng[bench].run(D[dn], pol, prm)
    except RuntimeError as ex: return dict(bench=bench, design=dn, tag=tag, collapse=True, err=str(ex))
    out = {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}
    if 'flog' in r:
        m = -r['flog'].min(1); out['fq'] = {str(q): float(np.percentile(m, q)) for q in (50, 90, 95, 97, 98, 99, 99.9)}
    out.update(bench=bench, design=dn, tag=tag, secs=time.time() - t0); return out
def configs(TH):
    J = []
    for dn in ('banks', 'decap'):
        J.append((dn, 'A', 'A', {}))
        for fc, th in (('fc7', TH[dn]), ('old', {})):
            F = dict(th, fc_mode='fc7') if fc == 'fc7' else {}
            J += [(dn, '%s reo' % fc, 'TOK', dict(F, reo=1, hold=0)),
                  (dn, '%s reo+hold d2' % fc, 'TOK', dict(F, reo=1, hold=1, maxd=2)),
                  (dn, '%s reo+hold d4' % fc, 'TOK', dict(F, reo=1, hold=1, maxd=4))]
    return J
if __name__ == '__main__':
    import multiprocessing as mp
    mode = sys.argv[1]
    if mode == 'cal':
        jobs = [(dn, 'cal ' + fc, 'A', dict(fc_log=1, **({'fc_mode': 'fc7'} if fc == 'fc7' else {})), b)
                for dn in ('banks', 'decap') for fc in ('fc7', 'old') for b in ('dhrystone', 'spmv')]
        out = HERE + 'steer7_cal.json'
    else:
        TH = json.load(open(HERE + 'steer7_th.json'))
        jobs = [c + (b,) for b in BENCHES for c in configs(TH)]; out = HERE + 'steer7.json'
    print(mode, len(jobs), 'jobs', flush=True); res = []; part = open(out + 'l', 'w')
    with mp.get_context('spawn').Pool(4) as P:
        for r in P.imap_unordered(job, jobs, chunksize=1):
            res.append(r); part.write(json.dumps(r, default=float) + '\n'); part.flush()
            print(r.get('design'), r.get('tag'), r.get('bench'), 'p999=%.2f' % r.get('p999', -1), 'sd=%.4f' % r.get('slowdown', -1), r.get('fq', ''), '%.0fs' % r.get('secs', 0), flush=True)
    json.dump(res, open(out, 'w'), indent=1, default=float); print('DONE', flush=True)
