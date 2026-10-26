import json
from multiprocessing import Pool
from sweep import setup
from replay2 import run, droop
from curmodel import FLOOR
G={}
def job(a):
    wl,pol,Qb,Ir,THF=a
    if wl not in G: G[wl]=setup(wl)
    inst,prods,cons,age,pin=G[wl]
    R=run(inst,prods,cons,age,pin,pol,100,Qb,Ir,dummy='fill',THF=THF)
    I=R['I']; dyn=(I-FLOOR).sum()-R['dummyq']
    return dict(wl=wl,pol=pol,Qb=Qb,THF=THF,slow=100*(R['makespan']/R['treal'].max()-1),droop=float(droop(R['In']).max()),E=100*R['dummyq']/dyn)
if __name__=='__main__':
    jobs=[(wl,pol,250,25,THF) for wl in ('mm','qsort','dhrystone') for pol in ('A','TD') for THF in (150,200,300,400)]
    with Pool(10) as p, open('sweep_thf.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs): f.write(json.dumps(r)+'\n'); f.flush()
