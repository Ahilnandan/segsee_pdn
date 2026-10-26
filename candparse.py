"""Parse BACT_CAND into uop instances: (iq, uopc, fu, rob, pdst, prs[], r_real, t_real).
r_real = first cycle the instance was a ready candidate, t_real = cycle it was granted."""
import re, gzip, pickle, sys, os, numpy as np
def parse(path):
    cache=path+'.cand2.pkl'
    if os.path.exists(cache): return pickle.load(open(cache,'rb'))
    open_={}; inst=[]
    with gzip.open(path,'rt',errors='ignore') as f:
        for ln in f:
            if not ln.startswith('BACT_CAND'): continue
            a=ln.split()
            if len(a)<12: continue
            tsc,iq,slot,g,u,fu,rob,pd,p1,p2,p3=map(int,a[1:12])
            k=(iq,rob,u)
            o=open_.get(k)
            if o is None or tsc-o[1]>1: o=[tsc,tsc]      # new instance, or it dropped out of the ready set (kill/flush/reuse)
            o[1]=tsc; open_[k]=o
            if g:
                inst.append((iq,u,fu,rob,pd,(p1,p2,p3),open_.pop(k)[0],tsc))
    inst.sort(key=lambda x:(x[7],x[0]))
    t=np.array([x[7] for x in inst]); g=np.diff(t); i=int(np.argmax(g))
    inst=inst[i+1:]                     # HTIF load-phase cut
    pickle.dump((inst,int(g[i])),open(cache,'wb'))
    return inst,int(g[i])
if __name__=='__main__':
    for wl in sys.argv[1:]:
        inst,gap=parse(f'traces/pm2_{wl}.log.gz')
        r=np.array([x[6] for x in inst]); t=np.array([x[7] for x in inst])
        from collections import Counter
        print(wl, len(inst),'grants; cut gap',gap,'; span',t[-1]-t[0],'cyc; iq mix',Counter(x[0] for x in inst),
              '; wait(t-r) mean %.2f p90 %.0f'%((t-r).mean(),np.percentile(t-r,90)))
        # per-cycle per-(iq,fu) max simultaneous grants = empirical port capacity
        cap=Counter(); per=Counter((x[7],x[0],x[2]) for x in inst)
        for (c,iq,fu),n in per.items(): cap[(iq,fu)]=max(cap[(iq,fu)],n)
        pq=Counter((x[7],x[0]) for x in inst); qc=Counter()
        for (c,iq),n in pq.items(): qc[iq]=max(qc[iq],n)
        print('  max grants/cycle per iq',dict(qc),' per (iq,fu_code)',dict(cap))
