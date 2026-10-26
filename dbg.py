import numpy as np
from collections import Counter
from replay2 import *
inst,_=parse('traces/pm2_mm.log.gz')
cap=Counter(); per=Counter((x[7],x[0],x[2]) for x in inst)
for (c,a,b),k in per.items(): cap[(a,b)]=max(cap[(a,b)],k)
prods,cons=build(inst)
R=run(inst,prods,cons,cap,'A',1e9,0,0)
d=R['tnew']-R['treal']; bad=np.nonzero(d)[0]
print(len(bad),'moved; early',(d<0).sum(),'late',(d>0).sum()); 
first=bad[np.argmin(R['treal'][bad])]
t0=min(x[6] for x in inst)
for i in sorted(bad,key=lambda i:R['treal'][i])[:6]:
    print(i,inst[i],'real',R['treal'][i],'new',R['tnew'][i],'prods',[(p,inst[p][7]-t0,R['tnew'][p]) for p in prods[i]])
# same-cycle over-subscription in real trace?
print('real per-(iq,fu) caps',dict(cap))
