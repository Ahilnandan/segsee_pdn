"""Validate the GA front on all 10 FULL execution traces + passive-only reference + jitter/guardband per design.
Usable-area gain = (passive reference area - design area) as % of the 4 mm^2 tile."""
from __future__ import print_function
import json, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn3, ga3, jitter3
WL = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
B = ga3.B

def full_eval(g):
    a, o, cv, per = ga3.evaluate(g, benches=WL, full=True)
    d = ga3.decode(g); spec = ga3.spec_of(d)
    z = np.load(B + 'qsort_rep_%s.npz' % d['var']); s = ga3.STARTS['qsort']; I = z['I'][s + 60000:s + 160000]
    pred = pdn3.predict_triggers(z['Ihat'][s + 60000:s + 160000], d['th_p'], d['lead']) if (d['inj'] and d['pred']) else None
    m = pdn3.simulate(I, dict(A_x=d['A_x'], R_d=d['R_d']), spec, pred=pred, want_v=True)
    j = jitter3.jitter(m['vavg'], m['vmin'])
    return dict(area=a, overhead=o, cv=cv, per=per, jitter=j, design={k: (v if not isinstance(v, np.ndarray) else v.tolist()) for k, v in d.items()})

def passive(Ax):
    g = np.zeros(ga3.NG)
    g[0:4] = np.clip(Ax * pdn3.SHARE / 0.5, 0, 1)                  # A_x[r] = Ax * share  (0.5 mm^2 per gene unit)
    return g

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == '--job':
    kind, v, out = sys.argv[2], sys.argv[3], sys.argv[4]
    run = json.load(open(ga3.HERE + 'ga3_run.json'))
    g = passive(float(v)) if kind == 'p' else np.array(run['X'])[int(v)]
    json.dump(full_eval(g), open(out, 'w'), default=float); sys.exit(0)

if __name__ == '__main__':
    run = json.load(open(ga3.HERE + 'ga3_run.json'))
    X = np.array(run['X']); F = run['F']
    feas = [i for i, f in enumerate(F) if f[2] == 0]
    # nondominated feasible set, thinned to <= 12 designs spread along the front
    nd = [i for i in feas if not any((F[j][0] <= F[i][0] and F[j][1] <= F[i][1]) and (F[j][0] < F[i][0] or F[j][1] < F[i][1]) for j in feas)]
    nd = sorted(set(nd), key=lambda i: F[i][0]); uniq = []
    for i in nd:
        if not uniq or abs(F[i][0] - F[uniq[-1]][0]) > 1e-3 or abs(F[i][1] - F[uniq[-1]][1]) > 1e-3: uniq.append(i)
    sel = uniq if len(uniq) <= 12 else [uniq[int(k)] for k in np.linspace(0, len(uniq) - 1, 12)]
    grid = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
    jobs = [('p', a) for a in grid] + [('d', i) for i in sel]
    import subprocess, tempfile
    tmp = ga3.HERE + 'val_jobs/'; os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PDN3_LIB='libpdn3b.so')
    procs = []; running = []
    for k, (kind, v) in enumerate(jobs):
        while len(running) >= 10:
            running = [p for p in running if p.poll() is None]; __import__('time').sleep(0.5)
        p = subprocess.Popen([sys.executable, __file__, '--job', kind, str(v), tmp + 'job%02d.json' % k], env=env,
                             stdout=subprocess.DEVNULL, stderr=open(tmp + 'job%02d.err' % k, 'w'))
        running.append(p); procs.append(p)
    for k, p in enumerate(procs): print("job", k, "exit", p.wait(), flush=True)
    res = [json.load(open(tmp + 'job%02d.json' % k)) for k in range(len(jobs))]
    pres, dres = res[:len(grid)], res[len(grid):]
    ref = next((a for a, r in zip(grid, pres) if r['cv'] == 0), None)
    print('passive-only reference (uniform decap by area share, policy A), all 10 full traces:')
    for a, r in zip(grid, pres):
        print('  decap %.2f mm2: worst p99.9 %.1f mV, worst max %.1f mV, feasible %s, jitter p99.9 %.1f ps, guardband %.2f %%' % (
            a, max(v['p999'] for v in r['per'].values()), max(v['max'] for v in r['per'].values()), r['cv'] == 0, r['jitter']['Jp_p999'], r['jitter']['guardband_p9999']))
    print('reference area meeting p99.9<=55 & max<=80 on all 10: %s mm2' % ref)
    print('\nGA front (fit on 5 windows) re-validated on all 10 full traces:')
    for i, r in zip(sel, dres):
        d = r['design']; wp = max(v['p999'] for v in r['per'].values()); wm = max(v['max'] for v in r['per'].values())
        parts = ['decap %.2f' % sum(d['A_x']), 'Rd %s' % '/'.join('%.2f' % x for x in d['R_d'])]
        if sum(d['Ca']) > 0: parts.append('activeC %.1fnF G%d amp%.1f' % (sum(d['Ca']) * 1e9, d['G'], d['amp']))
        if d['inj']: parts.append('inj %s L%.0fnH Vin%.1f Ron%.1f Ipk%.0fmA th%.0fmV%s' % ('pkg' if d['ind_pkg'] else 'die', d['L'] * 1e9, d['Vin'], d['Ron'], d['Ipk'] * 1e3, d['th_b'] * 1e3, ' pred(th%.0f,lead%d)' % (d['th_p'], d['lead']) if d['pred'] else ''))
        parts.append('policy %s' % d['var'])
        gain = (100 * (ref - r['area']) / pdn3.DIE_MM2) if ref else float('nan')
        print('  area %.3f mm2 (usable +%.1f %% of die) overhead %.2f %% | worst p99.9 %.1f max %.1f mV feasible %s | jitter p99.9 %.1f ps gb %.2f %% | %s' % (
            r['area'], gain, 100 * r['overhead'], wp, wm, r['cv'] == 0, r['jitter']['Jp_p999'], r['jitter']['guardband_p9999'], ', '.join(parts)))
    json.dump(dict(ref=ref, passive=list(zip(grid, pres)), front=dres), open(ga3.HERE + 'ga3_validate.json', 'w'), indent=1, default=float)
