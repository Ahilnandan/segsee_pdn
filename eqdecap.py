"""Droop of the UNSCHEDULED, UNBANKED baseline current vs passive on-die decap C.
Used to convert a policy's droop into 'equivalent passive decap'."""
import numpy as np, json, sys
from multiprocessing import Pool
from sweep import setup
from replay2 import run, droop
def one(a):
    wl,C=a
    inst,prods,cons,age,pin=setup(wl)
    R=run(inst,prods,cons,age,pin,'A',1e9,0,0)
    d=droop(R['I'],C=C*1e-9); return dict(wl=wl,C=C,max=float(d.max()),p999=float(np.percentile(d,99.9)))
if __name__=='__main__':
    jobs=[(wl,C) for wl in ('mm','qsort','dhrystone') for C in (5,7.5,10,15,20,30,40,60,80,120)]
    with Pool(10) as p: res=p.map(one,jobs)
    json.dump(res,open('eqdecap.json','w'))
    for r in res: print(r)
