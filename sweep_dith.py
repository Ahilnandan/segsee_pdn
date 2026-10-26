import json
from multiprocessing import Pool
from sweep import setup
from replay2 import run, droop
import numpy as np
G={}
def job(a):
    wl,pol,Qb,Ir=a
    if wl not in G: G[wl]=setup(wl)
    inst,prods,cons,age,pin=G[wl]
    R=run(inst,prods,cons,age,pin,pol,100,Qb,Ir)
    d=droop(R['In'])
    # resonance energy: spectrum of net current around the PDN resonance (100 MHz = period 10 cycles)
    I=R['In'][:R['makespan']]-R['In'][:R['makespan']].mean(); F=np.abs(np.fft.rfft(I)); f=np.fft.rfftfreq(len(I),1e-9)
    band=(f>80e6)&(f<125e6)
    return dict(wl=wl,pol=pol,Qb=Qb,slow=100*(R['makespan']/R['treal'].max()-1),droop=float(d.max()),p999=float(np.percentile(d,99.9)),
                res_energy=float((F[band]**2).sum()/(F[1:]**2).sum()),nsf=R['nsf'])
if __name__=='__main__':
    jobs=[(wl,pol,Qb,Ir) for wl in ('mm','qsort','dhrystone') for pol in ('A','DITH','TD','DITHD') for Qb,Ir in ((0,0),(250,25))]
    with Pool(8) as p, open('sweep_dith.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs): f.write(json.dumps(r)+'\n'); f.flush()
