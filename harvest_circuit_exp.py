import json, numpy as np
from multiprocessing import Pool
def job(a):
    wl,pol,ir,harv=a
    from sweep import setup
    from replay2 import run as replay
    from pulsim_harvest import run
    inst,prods,cons,age,pin=setup(wl)
    if pol=='A': R=replay(inst,prods,cons,age,pin,'A',1e9,0,0)
    else:        R=replay(inst,prods,cons,age,pin,'TD',100,250,ir,THF=150,recyc=(0.8 if harv else 0.0))
    I=R['I']   # pre-bank load current: the circuit applies its OWN bank
    d,hq=run(I,ir_mA=ir,harvest=harv)
    return dict(wl=wl,pol=pol,ir=ir,harv=harv,droop=float(d.max()),p999=float(np.percentile(d,99.9)),hq_nC=hq/1e3,
                slow=100*(R['makespan']/R['treal'].max()-1))
if __name__=='__main__':
    jobs=[(wl,pol,ir,h) for wl in ('qsort','dhrystone','mm') for pol in ('A','TD') for ir in (25,10) for h in (False,True)]
    with Pool(8) as p, open('harvest_circuit.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs): f.write(json.dumps(r)+'\n'); f.flush(); print(r,flush=True)
