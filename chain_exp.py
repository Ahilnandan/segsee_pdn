"""Run every scheduling policy's real current trace (FULL trace, no windowing) through
real switching power-delivery chains (pdnchain.py).

Policies (TH=100 mA, bank 250 pC / 25 mA where used):
  raw        : A, no bank            -- hardware as-is
  hold+bank  : TD + bank             -- token-gated hold (earlier result)
  dummy+bank : A + bank + fill       -- dummies only, zero holds
  combined   : TD + bank + fill      -- everything
Topologies:
  offchip    : case A, off-chip high-current 4-phase buck direct
  fivrF      : case B, off-chip 4-phase -> on-die 4-phase FIVR at F MHz (L scaled 2nH*140/F, constant ripple)
  fivrF+ff   : same + causal scheduler->FIVR feed-forward (3-cycle lead, ops already issued)
"""
import sys, json, pickle, os, numpy as np
from multiprocessing import Pool
POL={'raw':('A',1e9,0,0,None,None),'hold+bank':('TD',100,250,25,None,None),
     'dummy+bank':('A',100,250,25,'fill',150),'combined150':('TD',100,250,25,'fill',150),
     'combined200':('TD',100,250,25,'fill',200)}
def trace(wl,pol):
    f=f'traces/{wl}.{pol}.trace.pkl'
    if os.path.exists(f): return pickle.load(open(f,'rb'))
    from sweep import setup
    from replay2 import run
    inst,prods,cons,age,pin=setup(wl)
    p,TH,Qb,Ir,dm,THF=POL[pol]
    R=run(inst,prods,cons,age,pin,p,TH,Qb,Ir,dummy=dm,THF=THF,lead_rec=3)
    d=dict(In=R['In'].astype(np.float32),predL=R['predL'].astype(np.float32),makespan=R['makespan'],base=int(R['treal'].max()),dummyq=R['dummyq'],Isum=float(R['I'].sum()))
    pickle.dump(d,open(f,'wb')); return d
def job(a):
    wl,pol,topo,fsw,ff,cdie=a
    from pdnchain import simulate_chain
    d=trace(wl,pol); I=d['In'].astype(float)
    kw=dict(cdie=cdie*1e-9)
    if topo=='B':
        kw.update(case='B',fsw_fivr=fsw*1e6,L_fivr=2e-9*140/fsw)
        if ff: kw.update(ff_lead=3,ffsig_mA=d['predL'][:len(I)].astype(float))
    else: kw.update(case='A',dt=100e-12)
    o=simulate_chain(I,**kw)
    vmin,vmax=o['vmin'],o['vmax']
    W=5000; cs=np.cumsum(np.concatenate([[0],o['vavg']])); ma=np.empty_like(vmin)
    idx=np.arange(len(vmin)); lo=np.maximum(0,idx-W//2); hi=np.minimum(len(vmin),idx+W//2)
    ma=(cs[hi]-cs[lo])/(hi-lo)   # 5 us centred moving average = slow regulation level
    r=dict(wl=wl,pol=pol,topo=topo,fsw=fsw,ff=ff,cdie=cdie,droop=float(1e3*(0.9-vmin.min())),droop999=float(1e3*(0.9-np.percentile(vmin,0.1))),fast=float(1e3*(ma-vmin).max()),fast999=float(1e3*np.percentile(ma-vmin,99.9)),dc_drift=float(1e3*(0.9-ma.min())),
           over=float(1e3*(vmax.max()-0.9)),slow=100*(d['makespan']/d['base']-1),dummyE=100*d['dummyq']/max(1,d['Isum']-d['dummyq']))
    if 'v18min' in o: r['v18droop']=float(1e3*(1.8-o['v18min'].min()))
    return r
if __name__=='__main__':
    wls=('mm','qsort','dhrystone')
    for wl in wls:                       # build traces first (serial per wl, parallel across)
        pass
    with Pool(12) as p: p.starmap(trace,[(wl,pol) for wl in wls for pol in POL])
    jobs=[]
    for wl in wls:
        for pol in POL:
            jobs.append((wl,pol,'A',0,False,5))
            for fsw in (35,70,100,140,200):
                for ff in (False,True): jobs.append((wl,pol,'B',fsw,ff,5))
            for cdie in (2.5,3.5):     # less on-die decap at the 140 MHz point
                jobs.append((wl,pol,'A',0,False,cdie)); jobs.append((wl,pol,'B',140,False,cdie)); jobs.append((wl,pol,'B',140,True,cdie))
    out=sys.argv[1] if len(sys.argv)>1 else 'chain_exp.jsonl'
    done=set()
    if os.path.exists(out):
        for l in open(out): r=json.loads(l); done.add((r['wl'],r['pol'],r['topo'],r['fsw'],r['ff'],r['cdie']))
    jobs=[j for j in jobs if j not in done]
    print(len(jobs),'jobs',flush=True)
    with Pool(12) as p, open(out,'a') as f:
        for r in p.imap_unordered(job,jobs):
            f.write(json.dumps(r)+'\n'); f.flush(); print(r['wl'],r['pol'],r['topo'],r['fsw'],r['ff'],r['cdie'],'%.1f'%r['droop'],flush=True)
