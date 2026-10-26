"""Closed-loop full system: scheduler (replay2 as generator, sensed tokens + supply signal) <->
switching supply chain (pdnchain) with passive die decap or an ACTIVE bank (added, or replacing
passive decap at equal total capacitance)."""
import sys, json, os, numpy as np
from multiprocessing import Pool
CFG = {  # name: (die passive nF, bank nF or 0)
    'passive5': (5.0, 0.0), 'add0.5': (5.0, 0.5), 'repl4.5+0.5': (4.5, 0.5), 'repl4+1': (4.0, 1.0)}
POL = {'A': ('A', None), 'T': ('T', None), 'SF': ('SF', None), 'TD': ('TD', None), 'SFD': ('SFD', None), 'SFD+fill': ('SFD', 'fill')}
def job(a):
    wl, case, cfg, pol, sg = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    cd, cb = CFG[cfg]; p, dm = POL[pol]
    ir = 25.0 if cb else 0.0
    g = run_gen2(inst, prods, cons, age, pin, p, 100, 0, ir, dummy=dm, THF=200)
    kw = dict(case=case, sched=g, nmax=n, cdie=cd * 1e-9, sig_gain=sg)
    kw.update(dict(fsw_fivr=140e6, L_fivr=2e-9) if case == 'B' else dict(dt=100e-12))
    if cb: kw['bank'] = dict(cb=cb, ir=ir)
    o = simulate_chain(None, **kw)
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']
    return dict(wl=wl, case=case, cfg=cfg, pol=pol, sg=sg, droop=float(v.max()), p999=float(np.percentile(v, 99.9)),
                top20=float(np.sort(v)[-20:].mean()), over=float(1e3 * (o['vmax'].max() - 0.9)),
                slow=100 * (R['makespan'] / R['treal'].max() - 1), nhold=R['nhold'], nsf=R.get('nsf', 0),
                dummyq=R['dummyq'], Isum=float(R['I'].sum()))
if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'fullsys.jsonl'
    jobs = [(wl, 'B', cfg, pol, 1.0) for wl in ('qsort', 'dhrystone', 'mm') for cfg in CFG for pol in POL]
    jobs += [(wl, 'A', cfg, pol, 1.0) for wl in ('qsort', 'dhrystone', 'mm') for cfg in ('passive5', 'repl4.5+0.5') for pol in POL]
    done = set()
    if os.path.exists(out):
        for l in open(out): r = json.loads(l); done.add((r['wl'], r['case'], r['cfg'], r['pol'], r['sg']))
    jobs = [j for j in jobs if j not in done]; print(len(jobs), 'jobs', flush=True)
    with Pool(12) as p, open(out, 'a') as f:
        for r in p.imap_unordered(job, jobs):
            f.write(json.dumps(r) + '\n'); f.flush(); print(r['wl'], r['case'], r['cfg'], r['pol'], '%.1f/%.1f %.2f%%' % (r['droop'], r['p999'], r['slow']), flush=True)
