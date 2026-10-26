import time, numpy as np, sys
from smt_replay import make_thread, run
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
# sanity 1: a single thread on a huge-IQ machine with BOOM-like caps should track its solo timeline
th = make_thread('dhrystone', 50000, N)
t0 = time.time(); R = run([th], 'A', iq=10**6); print('solo dhry on SMT machine: makespan', R['makespan'], 'solo span', int(th['tr'].max()), f'{time.time()-t0:.1f}s')
thr = [make_thread('qsort', 50000, N), make_thread('dhrystone', 50000, N), make_thread('mm', 50000, N), make_thread('mm', 200000, N)]
for pol in ('A', 'T', 'FLAT'):
    t0 = time.time(); R = run(thr, pol)
    nr = R['nready']; d = droop(R['I']) if False else None
    print(pol, 'makespan', R['makespan'], 'fin', R['fin'], f"ready mean {nr.mean():.1f} p50 {np.median(nr):.0f} max {nr.max()} <9: {100*(nr<9).mean():.1f}%",
          'nchg', R['nchg'], f'{time.time()-t0:.1f}s')
