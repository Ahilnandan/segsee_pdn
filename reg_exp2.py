import json, numpy as np
from multiprocessing import Pool
CTL={'loopON':dict(ctl='vm'),'loopOFF':dict(ctl='open'),'slowint':dict(ctl='vm',kp_vm=0.0,fint=0.5e6)}
def job(a):
    wl, fsw, ctl, pol, sth = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    g = run_gen2(inst, prods, cons, age, pin, pol, 100, 0, 0, reg_sth=sth)
    o = simulate_chain(None, case='B', sched=g, nmax=n, fsw_fivr=fsw * 1e6, L_fivr=2e-9 * 140 / fsw, **CTL[ctl])
    v = 1e3 * (0.9 - o['vmin']); R = o['sched']
    return dict(wl=wl, fsw=fsw, ctl=ctl, pol=pol, sth=sth, droop=float(v.max()), p999=float(np.percentile(v, 99.9)),
                over=float(1e3 * (o['vmax'].max() - 0.9)), vavg=float(o['vavg'].mean()),
                slow=100 * (R['makespan'] / R['treal'].max() - 1), dummyE=100 * R['dummyq'] / (R['I'].sum() - R['dummyq']))
if __name__ == '__main__':
    jobs=[]
    for wl in ('qsort','dhrystone','mm'):
        for fsw in (140,70,35):
            jobs.append((wl,fsw,'loopON','A',5.0))
            for ctl in ('loopOFF','slowint'):
                jobs += [(wl,fsw,ctl,'REG',5.0),(wl,fsw,ctl,'REG',40.0),(wl,fsw,ctl,'REGH',5.0),(wl,fsw,ctl,'A',5.0)]
    with Pool(12) as p, open('reg_exp2.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
