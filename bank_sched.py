"""Does a switched-cap bank keep up with real event rates, and can issue reordering
space events out enough for the recharge to keep up?
Bank = token bucket: capacity Qb (pC), refill Ir (mA -> pC/cycle @1GHz), each
'event' (5-cycle current rise > TH) needs q = rise*W/2 pC (taper injection over W)."""
import sys, numpy as np
from collections import deque
from curmodel import parse_iss, current_from_times, SHAPE, COST, port, PORTW, FLOOR
W=5; SLEN=10

def events(I, th):
    out=[]; last=-10**9
    for t in range(W, len(I)):
        r = I[t]-I[t-W:t].min()
        if r>th and t-last>=W: out.append((t,r)); last=t
    return out

def bucket(I, th, Qb, Ir):
    ev=events(I,th); tok=Qb; tlast=0; unc=0.0; full=0; tot=0.0
    for t,r in ev:
        tok=min(Qb, tok+Ir*(t-tlast)); tlast=t
        q=r*W/2; s=min(q,tok); tok-=s; unc+=q-s; tot+=q; full+= (s>=q-1e-9)
    return len(ev), full, unc, tot

def replay(evs, pol, th, Qb, Ir, look=12, max_defer=6):
    t0=evs[0][0]; n=evs[-1][0]-t0+40; H=n+2000
    I=np.full(H,FLOOR); pools={p:deque() for p in PORTW}
    ei=0; N=len(evs); t=0; last=0; tok=Qb
    while ei<N or any(pools.values()):
        while ei<N and evs[ei][0]-t0<=t:
            c,u=evs[ei]; pools[port(u)].append((c-t0,u)); ei+=1
        tok=min(Qb,tok+Ir)
        lo=I[max(0,t-W):t].min() if t>0 else FLOOR
        for pc,Wd in PORTW.items():
            pool=pools[pc]
            if not pool: continue
            cand=list(pool)[:look]
            if pol=='A':
                pick=cand[:Wd]
            else:
                forced=[p for p in cand if t-p[0]>=max_defer]
                def need(sel):   # charge the bank would need if we issue `sel` now
                    tmp=I[t:t+SLEN].copy()
                    for _,u in sel: tmp[:len(SHAPE[u])]+=SHAPE[u]
                    r=tmp.max()-lo
                    return r*W/2 if r>th else 0.0
                oldest=cand[:Wd]
                if need(oldest)<=tok: pick=oldest                     # bank can cover it: no change
                else:
                    cheap=sorted(cand,key=lambda p:COST[p[1]])[:Wd]    # T: pure reorder, still fill width
                    pick=cheap
                    if pol=='TD' and need(cheap)>tok:                  # TD: may also hold ops back
                        pick=list(forced[:Wd])
                        for p in cheap:
                            if len(pick)>=Wd: break
                            if p in pick: continue
                            if need(pick+[p])<=tok: pick.append(p)
                    for p in forced:                                   # fairness backstop always wins
                        if p not in pick and len(pick)<Wd: pick.append(p)
            q=0.0
            for p in pick:
                pool.remove(p); _,u=p; I[t:t+len(SHAPE[u])]+=SHAPE[u]; last=t
            r=I[t:t+SLEN].max()-lo
            if r>th: tok-=min(tok, r*W/2)
        t+=1
    return I[:last+SLEN], last+1

if __name__=='__main__':
    wl=sys.argv[1]; path=f'traces/pm2_{wl}.log.gz'
    evs,gap=parse_iss(path)
    print(f"=== {wl}: {len(evs)} issue events after HTIF cut (load-phase gap {gap} cyc) ===")
    Ihw=current_from_times(evs)
    IA,cA=replay(evs,'A',0,0,0)
    print(f"apparatus check: replay-A vs real HW issue times -> max |dI| {np.abs(IA[:len(Ihw)]-Ihw[:len(IA)]).max():.2f} mA, "
          f"cycles {cA} vs {len(Ihw)-20}")
    r=np.array([Ihw[t]-Ihw[t-W:t].min() for t in range(W,len(Ihw))])
    print("5-cycle rise pctl (mA): p50 %.0f p90 %.0f p99 %.0f p99.9 %.0f max %.0f" % tuple(np.percentile(r,[50,90,99,99.9,100])))
    for th in [60,100]:
        e=events(Ihw,th); ts=np.array([x[0] for x in e]); g=np.diff(ts) if len(ts)>1 else np.array([0])
        q=np.array([x[1] for x in e])*W/2
        print(f"\nTH={th} mA: {len(e)} events = 1 per {len(Ihw)/max(1,len(e)):.0f} cycles; gap p10/p25/p50 "
              f"{np.percentile(g,10):.0f}/{np.percentile(g,25):.0f}/{np.percentile(g,50):.0f} cyc; "
              f"q p50 {np.median(q):.0f} pC p99 {np.percentile(q,99):.0f} pC; mean demand {q.sum()/len(Ihw):.2f} mA")
        print(f"  {'bank':>14} {'Ir':>5} | {'policy':>6} {'events':>7} {'fully served':>13} {'uncovered':>10} {'slowdown':>9}")
        for Qb,Ir in [(250,5),(250,10),(250,25),(500,10)]:
            for pol in ['A','T','TD']:
                I,c = (IA,cA) if pol=='A' else replay(evs,pol,th,Qb,Ir)
                ne,full,unc,tot=bucket(I,th,Qb,Ir)
                print(f"  {Qb:>11} pC {Ir:>4}mA | {pol:>6} {ne:7d} {100*full/max(1,ne):12.1f}% "
                      f"{100*unc/max(1e-9,tot):9.1f}% {100*(c-cA)/cA:8.2f}%")
