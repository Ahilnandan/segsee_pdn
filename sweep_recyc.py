"""Edge-harvesting bank vs dummy fill, same falling-edge rule (THF), same bank (250 pC).
Energy accounting: dummies burn their charge; harvesting loses only (1-eta) of the sunk charge."""
import json
from multiprocessing import Pool
from sweep import setup
from replay2 import run, droop
from curmodel import FLOOR
G={}
def job(a):
    wl,pol,Qb,Ir,THF,eta,mode=a
    if wl not in G: G[wl]=setup(wl)
    inst,prods,cons,age,pin=G[wl]
    if mode=='harvest': R=run(inst,prods,cons,age,pin,pol,100,Qb,Ir,THF=THF,recyc=eta)
    else:               R=run(inst,prods,cons,age,pin,pol,100,Qb,Ir,dummy='fill',THF=THF)
    I=R['I']; dyn=(I-FLOOR).sum()-R['dummyq']
    waste=R['dummyq'] if mode!='harvest' else (1-eta)*R['rq']
    return dict(wl=wl,pol=pol,Qb=Qb,Ir=Ir,THF=THF,eta=eta,mode=mode,slow=100*(R['makespan']/R['treal'].max()-1),
                droop=float(droop(R['In']).max()),Ewaste=100*waste/dyn,harv_nC=R['rq']/1e3,nhold=R['nhold'],events=R['events'],evfull=R['evfull'])
if __name__=='__main__':
    jobs=[]
    for wl in ('mm','qsort','dhrystone'):
        for THF in (100,150,200):
            for Ir in (25,10,0):          # 0 = NO external recharge at all: bank lives only on harvested charge
                jobs.append((wl,'TD',250,Ir,THF,0.8,'harvest'))
                jobs.append((wl,'A',250,Ir,THF,0.8,'harvest'))
    with Pool(10) as p, open('sweep_recyc.jsonl','w') as f:
        for r in p.imap_unordered(job,jobs): f.write(json.dumps(r)+'\n'); f.flush()
