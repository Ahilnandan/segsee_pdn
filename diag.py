import numpy as np
from collections import deque
from curmodel import parse_iss, current_from_times, port, PORTW
evs,_=parse_iss('traces/pm2_mm.log.gz')
# (1) how often does a port's ready pool exceed its real issue width?  (pool built from real issue times)
t0=evs[0][0]; pools={p:deque() for p in PORTW}; ei=0; N=len(evs); t=0; dec={p:0 for p in PORTW}; more={p:0 for p in PORTW}
while ei<N or any(pools.values()):
    while ei<N and evs[ei][0]-t0<=t: c,u=evs[ei]; pools[port(u)].append(u); ei+=1
    for p,W in PORTW.items():
        if pools[p]:
            dec[p]+=1; more[p]+= len(pools[p])>W
            for _ in range(min(W,len(pools[p]))): pools[p].popleft()
    t+=1
for p in PORTW: print(f"(1) port {p}: {dec[p]} issue decisions, pool > width in {100*more[p]/max(1,dec[p]):.2f}%")
print("    note: replay pools are filled from REAL issue times, so they only contain what the real scheduler")
print("    already chose that cycle -- the true ready set (BACT_CAND) is larger. See discussion.")
# (3) band-limited events: bank only needs to cover the envelope the passive decap can't filter
I=current_from_times(evs)
def ma(x,k): return np.convolve(x,np.ones(k)/k,mode='same')
for k in [1,5,10]:
    If=ma(I,k) if k>1 else I
    r=np.array([If[t]-If[t-5:t].min() for t in range(5,len(If))])
    print(f"(3) envelope smoothed over {k:2d} cyc: 5-cyc rise p50 {np.percentile(r,50):5.0f} p90 {np.percentile(r,90):5.0f} p99 {np.percentile(r,99):5.0f} mA")
