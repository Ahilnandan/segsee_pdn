import numpy as np, pdn3, sys, json
from multiprocessing import Pool
B='/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
WL=['qsort','rsort','mm','dhrystone','spmv','multiply','memcpy','median','vvadd','towers']
def exec_start(I):
    import pandas as pd
    m=pd.Series(I.sum(1)).rolling(100).mean().values; m=np.where(np.isnan(m),0,m); m[:3000]=0; return max(0,int(np.argmax(m>150))-100)
def job(a):
    w,Ax=a; I=np.load(B+w+'_rep_A.npz')['I'].astype(float); s=exec_start(I)
    m=pdn3.simulate(I[s:],dict(A_x=[Ax*x for x in pdn3.SHARE]))
    return w,Ax,s,len(I)-s,m['p99'],m['p999'],m['max'],m['over_max']
if __name__=='__main__':
    jobs=[(w,a) for w in WL for a in (0.0,0.2,0.4,0.8)]
    with Pool(20) as P: res=P.map(job,jobs)
    print('%-10s %6s %8s %8s | p99 p99.9 max over (mV)'%('bench','decap','start','len'))
    for r in res: print('%-10s %6.2f %8d %8d | %5.1f %5.1f %5.1f %5.1f'%r)
    json.dump(res,open('passive_scan.json','w'))
