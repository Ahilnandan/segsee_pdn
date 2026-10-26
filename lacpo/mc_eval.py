"""mc_eval: gate + scenario sweep for the multi-core cross-core issue coordinator (mc_engine).
  python mc_eval.py gate            single core through mc_engine.Core + pdn6.Sim vs policy6.Engine.run(.., 'A', ..) on a
                                    20k window (rob 0 and rob 80). policy6.py line 83 has its array setup swallowed by a
                                    comment (NameError on BG); the gate execs an in-memory copy with that line split.
  python mc_eval.py sweep [design] [workers<=4]   scenarios x coordinator settings, 60k-cycle windows -> mc_eval_<design>.json
Scenarios: lock = 4x qsort lockstep; jit = 4x qsort start offsets U(0..20); diff = qsort rsort mm dhrystone;
           rand = 4x qsort random offsets U(0..20000) (reference)."""
from __future__ import print_function
import os, sys, json, types, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'; sys.path.insert(0, HERE)
A0 = 70000; W = 60000          # all four benches are past the HTIF load phase from ~60k


def policy6_fixed():
    src = open(HERE + 'policy6.py').read()
    bad = "# room for up to ~200 % slowdown before 'overflow' ring ="
    assert bad in src
    src = src.replace(bad, "# room for up to ~200 % slowdown before 'overflow'\n        ring =")
    mod = types.ModuleType('policy6fix'); mod.__file__ = HERE + 'policy6.py'
    exec(compile(src, HERE + 'policy6.py', 'exec'), mod.__dict__); return mod


def gate():
    import eval6, mc_engine as M
    P6 = policy6_fixed(); D, prm = eval6.designs(); E = P6.Engine('qsort')
    for dn in ('decap', 'banks'):
        for rob in (0, 80):
            r0 = E.run(D[dn], 'A', dict(prm, rob=rob), A0, A0 + 20000, record=True)
            r1 = M.run_single(E, A0, A0 + 20000, D[dn], rob=rob)
            print('GATE %-6s rob %2d  policy6 p999 %.4f max %.4f sd %.4f | stepper p999 %.4f max %.4f sd %.4f | tnew identical %s'
                  % (dn, rob, r0['p999'], r0['max'], r0['slowdown'], r1['p999'], r1['max'], r1['slowdown'], bool((r0['tnew'] == r1['tnew']).all())), flush=True)


SCEN = {'lock': (['qsort'] * 4, [0, 0, 0, 0]),
        'jit': (['qsort'] * 4, list(np.random.RandomState(5).randint(0, 21, 4))),
        'diff': (['qsort', 'rsort', 'mm', 'dhrystone'], [0, 0, 0, 0]),
        'rand': (['qsort'] * 4, list(np.random.RandomState(3).randint(0, 20000, 4)))}
CONF = [('none', None), ('B80', dict(budget=80.0)), ('B50', dict(budget=50.0)), ('B30', dict(budget=30.0)), ('B20', dict(budget=20.0)),
        ('B50+sens', dict(budget=50.0, ks=2.0, v0=25.0))]
_eng = {}


def job(a):
    dn, sc, tag, coord = a
    import policy6 as P6, eval6, mc_engine as M
    ws, offs = SCEN[sc]
    for w in ws:
        if w not in _eng: _eng[w] = P6.Engine(w)
    D, _ = eval6.designs()
    r = M.run_mc([_eng[w] for w in ws], [A0 + int(o) for o in offs], D[dn], W=W, coord=coord)
    r.update(design=dn, scen=sc, tag=tag); print(json.dumps(r), flush=True); return r


if __name__ == '__main__':
    if sys.argv[1] == 'gate': gate(); sys.exit()
    import multiprocessing as mp
    dn = sys.argv[2] if len(sys.argv) > 2 else 'decap'; nw = min(4, int(sys.argv[3]) if len(sys.argv) > 3 else 4)
    jobs = [(dn, sc, tag, c) for sc in ('lock', 'jit', 'diff', 'rand') for tag, c in CONF]
    with mp.get_context('spawn').Pool(nw) as P: res = P.map(job, jobs, chunksize=1)
    json.dump(res, open(HERE + 'mc_eval_%s.json' % dn, 'w'), indent=1)
    base = {r['scen']: r for r in res if r['tag'] == 'none'}
    print('%-5s %-9s %6s %6s %7s %22s %7s' % ('scen', 'conf', 'p99.9', 'max', 'dI999', 'sd vs none mean/worst %', 'held'))
    for r in res:
        b = base[r['scen']]; rel = [(e / max(1, e0) - 1) * 100 for e, e0 in zip(r['exec'], b['exec'])]
        print('%-5s %-9s %6.1f %6.1f %7.0f %10.2f / %6.2f %9d' % (r['scen'], r['tag'], r['p999'], r['max'], r['dI999'], np.mean(rel), max(rel), r['nheld']))
