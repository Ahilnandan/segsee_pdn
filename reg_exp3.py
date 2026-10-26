import json, numpy as np
from multiprocessing import Pool
V={ 'REG':        dict(pol='REG', sm='both', sth=5.0, bank=None),
    'REG-vonly':  dict(pol='REG', sm='v', sth=33.0, bank=None),
    'REGH+bank':  dict(pol='REGH', sm='both', sth=5.0, bank=dict(cb=0.5, ir=25.0)),
    'REGH+bidir': dict(pol='REGH', sm='both', sth=5.0, bank=dict(cb=0.5, ir=25.0, harvest_v=True)),
    'REG+bidir':  dict(pol='REG', sm='both', sth=60.0, bank=dict(cb=0.5, ir=25.0, harvest_v=True)),
    'A+bidir':    dict(pol='A', sm='both', sth=5.0, bank=dict(cb=0.5, ir=25.0, harvest_v=True)),
    'A+bank':     dict(pol='A', sm='both', sth=5.0, bank=dict(cb=0.5, ir=25.0)) }
def job(a):
    wl, fsw, var = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    v = V[var]
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    g = run_gen2(inst, prods, cons, age, pin, v['pol'], 100, 0, 0, reg_sth=v['sth'])
    kw = dict(case='B', sched=g, nmax=n, fsw_fivr=fsw * 1e6, L_fivr=2e-9 * 140 / fsw, ctl='vm', kp_vm=0.0, fint=0.5e6, sig_mode=v['sm'])
    if v['bank']: kw.update(bank=v['bank'], cdie=4.5e-9)      # iso-capacitance: bank REPLACES 0.5 nF of passive
    o = simulate_chain(None, **kw)
    w = 1e3 * (0.9 - o['vmin']); R = o['sched']
    return dict(wl=wl, fsw=fsw, var=var, droop=float(w.max()), p999=float(np.percentile(w, 99.9)), over=float(1e3 * (o['vmax'].max() - 0.9)),
                vavg=float(o['vavg'].mean()), slow=100 * (R['makespan'] / R['treal'].max() - 1), dummyE=100 * R['dummyq'] / (R['I'].sum() - R['dummyq']))
if __name__ == '__main__':
    jobs = [(wl, fsw, var) for wl in ('qsort', 'dhrystone', 'mm') for fsw in (140, 70, 35) for var in V]
    with Pool(12) as p, open('reg_exp3.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
