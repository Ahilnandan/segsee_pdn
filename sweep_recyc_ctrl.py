import json
from multiprocessing import Pool
from sweep import setup
from replay2 import run, droop
G={}
def job(a):
    wl,Qb,Ir=a
    if wl not in G: G[wl]=setup(wl)
    inst,prods,cons,age,pin=G[wl]
    R=run(inst,prods,cons,age,pin,'TD',100,Qb,Ir)
    return dict(wl=wl,Qb=Qb,Ir=Ir,slow=100*(R['makespan']/R['treal'].max()-1),droop=float(droop(R['In']).max()),events=R['events'],evfull=R['evfull'])
if __name__=='__main__':
    jobs=[(wl,Qb,Ir) for wl in ('mm','qsort','dhrystone') for Qb in (250,500,1000) for Ir in (0,10,25)]
    with Pool(9) as p, open('sweep_recyc_ctrl.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs): f.write(json.dumps(r)+'\n'); f.flush()
