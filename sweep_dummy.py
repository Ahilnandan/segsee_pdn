"""Dummy-instruction sweep on the validated replay (regression-checked: dummy=None
reproduces sweep.jsonl exactly). Adds overshoot (voltage above nominal) because
falling edges -- what 'fill' targets -- show up there, not in droop."""
import sys, json, numpy as np
from multiprocessing import Pool
from sweep import setup
from replay2 import run, droop
from curmodel import FLOOR
G={}
def job(a):
    wl,pol,TH,Qb,Ir,dm=a
    if wl not in G: G[wl]=setup(wl)
    inst,prods,cons,age,pin=G[wl]
    R=run(inst,prods,cons,age,pin,pol,TH,Qb,Ir,dummy=dm)
    d=droop(R['In'])
    I=R['I']; dyn=(I-FLOOR).sum()-R['dummyq']
    return dict(wl=wl,pol=pol,TH=TH,Qb=Qb,Ir=Ir,dummy=dm or 'none',makespan=R['makespan'],base=int(R['treal'].max()),
        droop_max=float(d.max()),over_max=float((-d).max()),events=R['events'],evfull=R['evfull'],
        dummy_pct_dyn=100*R['dummyq']/dyn, dummy_pct_tot=100*R['dummyq']/(I.sum()-R['dummyq']),ndum=R['ndum'],nhold=R['nhold'])
if __name__=='__main__':
    wls=('mm','qsort','dhrystone'); jobs=[]
    for wl in wls:
        for TH in (60,100):
            for Qb,Ir in ((0,0),(250,25)):
                for dm in (None,'fill'): jobs.append((wl,'A',TH,Qb,Ir,dm))
                for dm in (None,'fill','ramp','fill+ramp'): jobs.append((wl,'TD',TH,Qb,Ir,dm))
    with Pool(12) as p, open('sweep_dummy.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs):
            f.write(json.dumps(r)+'\n'); f.flush(); print(r['wl'],r['pol'],r['TH'],r['Qb'],r['dummy'],'done',flush=True)
