import sys, json, numpy as np, pickle, os
from multiprocessing import Pool
from replay2 import parse, build, unwrap_age, calibrate, run, droop
G={}
def setup(wl):
    f=f'traces/{wl}.setup.pkl'
    if os.path.exists(f): return pickle.load(open(f,'rb'))
    inst,_=parse(f'traces/pm2_{wl}.log.gz'); prods,cons=build(inst); age=unwrap_age(inst)
    pin,_=calibrate(inst,prods,cons,age); d=(inst,prods,cons,age,pin); pickle.dump(d,open(f,'wb')); return d
def job(a):
    wl,pol,TH,Qb,Ir,vr=a
    if wl not in G: G[wl]=setup(wl)
    inst,prods,cons,age,pin=G[wl]
    R=run(inst,prods,cons,age,pin,pol,TH,Qb,Ir,vdd_recharge=vr)
    dr=droop(R['In']); d=R['tnew']-R['treal']
    return dict(wl=wl,vr=vr,pol=pol,TH=TH,Qb=Qb,Ir=Ir,makespan=R['makespan'],base=int(R['treal'].max()),
        events=R['events'],evfull=R['evfull'],unc_nC=R['unc']/1e3,drawn_nC=R['drawn']/1e3,
        droop_max=float(dr.max()),droop_p999=float(np.percentile(dr,99.9)),
        mean_delay=float(d.mean()),p99_delay=float(np.percentile(d,99)),frac_moved=float((d!=0).mean()),
        nsel=R['nsel'],nchg=R['nchg'],nhold=R['nhold'])
if __name__=='__main__':
    wls=sys.argv[1:]
    for wl in wls: setup(wl)
    jobs=[(wl,pol,TH,Qb,Ir,False) for wl in wls for TH in (60,100)
          for (Qb,Ir) in [(0,0),(250,5),(250,10),(250,25),(500,10),(500,25)] for pol in ('A','T','TD')]
    jobs+=[(wl,pol,TH,Qb,Ir,True) for wl in wls for TH in (60,100) for (Qb,Ir) in [(250,10),(250,25)] for pol in ('A','TD')]
    with Pool(12) as p, open('sweep.jsonl','w') as f:
        for res in p.imap_unordered(job,jobs):
            f.write(json.dumps(res)+'\n'); f.flush()
            print(res['wl'],res['pol'],res['TH'],res['Qb'],res['Ir'],res['vr'],'done',flush=True)
