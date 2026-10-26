"""Dependency-aware counterfactual issue replay on real BOOM ready sets (BACT_CAND),
with an on-die switched-cap bank (token bucket) and a series-RLC droop check.

Uop instance i: ready r_i and grant t_i from the real trace; producers = latest instance
with pdst==prs issued at t_p <= r_i within LWIN cycles (conservative: may include false
deps across int/FP files, which only makes reordering look MORE expensive).
Counterfactual ready'_i = r_i + max(0, max_p (t'_p - t_p)); never earlier than r_i.

Bank: holds I_net[t] <= min(I_net[t-W:t]) + TH (a slew clamp on the fast edge) while it
has charge; each cycle refills Ir pC (1 GHz core) up to Qb. Charge drawn = injected mA*1ns.
Policies (per cycle, all issue queues jointly, width/FU caps from the trace):
  A  : oldest-first (priority = real grant time -> reproduces hardware exactly)
  T  : if the bank cannot cover A's pick, choose the same NUMBER of uops minimising the
       charge the bank would need (pure reorder, work-conserving)
  TD : T, then hold uops back until the need fits the tokens; a uop held MAXD times is forced
  Bank Qb=0 + TD == pipeline damping (P1-style slew limiting, no active decap)."""
import sys, numpy as np, heapq
from collections import defaultdict, Counter
from scipy.linalg import expm
from candparse import parse
from curmodel import SHAPE, FLOOR
W=5; LWIN=10**9; MAXD=6; H=6   # LWIN: no window -- a preg cannot be reallocated while a consumer is in flight
WIDTH={1:4,2:2,4:2}
PORTW_INT=4
BK=4   # bridge lookahead (cycles): wakeups already scheduled in the wake table

def build(inst):
    n=len(inst); last={}; prods=[[] for _ in range(n)]
    order=sorted(range(n),key=lambda i:inst[i][7])
    # producers: process in real grant order; lookup uses latest issued pdst writer with t_p<=r_i
    hist=defaultdict(list)                       # pdst -> list of (t_p, idx) increasing
    for i in order:
        iq,u,fu,rob,pd,prs,r,t=inst[i]
        for p in prs:
            if p==0: continue
            h=hist.get(p)
            if not h: continue
            k=len(h)-1
            while k>=0 and h[k][0]>=r: k-=1
            if k>=0 and r-h[k][0]<=LWIN: prods[i].append(h[k][1])
        if pd!=0: hist[pd].append((t,i))
    cons=[[] for _ in range(n)]
    for i,ps in enumerate(prods):
        for p in set(ps): cons[p].append(i)
    return prods,cons

# Real PaperMatch port map (execution-units.scala is_nth with int_width=4; fp-pipeline w==0)
FU_ALU,FU_JMP,FU_MEM,FU_MUL,FU_DIV,FU_CSR,FU_FPU,FU_FDV,FU_I2F,FU_F2I=1,2,4,8,16,32,64,128,256,512
PORTS={1:[FU_ALU|FU_JMP|FU_I2F, FU_ALU|FU_CSR, FU_ALU|FU_MUL, FU_ALU|FU_DIV],
       2:[FU_MEM, FU_MEM],
       4:[FU_FPU|FU_FDV|FU_F2I, FU_FPU]}

def unwrap_age(inst, R=384):
    """program-order age from wrapped rob_idx (in-flight span < ROB size)."""
    order=sorted(range(len(inst)),key=lambda i:inst[i][7]); age=[0]*len(inst); ref=None
    for i in order:
        rb=inst[i][3]
        if ref is None: v=rb
        else:
            k=(ref-rb)//R; v=min((rb+(k+j)*R for j in (-1,0,1,2)),key=lambda x:abs(x-ref))
        age[i]=v; ref=v if ref is None else max(ref,v)
    return age

def port_fill(order, fu, iq, busy_ok):
    """BOOM greedy: each uop (age order) takes the lowest free port whose fu mask matches."""
    used={q:[False]*len(p) for q,p in PORTS.items()}; out=[]
    for i in order:
        q=iq[i]; P=PORTS[q]
        for w in range(len(P)):
            if not used[q][w] and (fu[i] & busy_ok(q,w,P[w])):
                used[q][w]=True; out.append(i); break
    return out

def feasible(sel, fu, iq, busy_ok):
    """exact bipartite check: can every uop in sel get a distinct compatible port?"""
    by={}
    for i in sel: by.setdefault(iq[i],[]).append(i)
    for q,us in by.items():
        P=[busy_ok(q,w,m) for w,m in enumerate(PORTS[q])]
        if len(us)>len(P): return False
        match=[-1]*len(P)
        def aug(u,seen):
            for w in range(len(P)):
                if fu[u]&P[w] and not seen[w]:
                    seen[w]=True
                    if match[w]<0 or aug(match[w],seen): match[w]=u; return True
            return False
        for u in us:
            if not aug(u,[False]*len(P)): return False
    return True

DIVLAT=None
DSH=np.array([10.,20.,20.,10.])   # dummy op = ALU-class shape (curmodel 'alu')
def run(inst, prods, cons, age, pin, pol, TH, Qb, Ir, vdd_recharge=False, rob=384, dummy=None, THF=None, lead_rec=0, recyc=0.0, sf_db=5.0, reg_k=2.0, reg_sth=5.0, mp_h=None, mp_vth=None, mp_Hp=13, mp_maxd=None, wb=None):
    n=len(inst); t0=min(x[6] for x in inst)
    r=np.array([x[6] for x in inst])-t0; treal=np.array([x[7] for x in inst])-t0
    r=np.where(pin, treal, r)                       # pinned: unobservable structural stall persists
    iq=[x[0] for x in inst]; fu=[x[2] for x in inst]; sh=[SHAPE.get(x[1],SHAPE[0]) for x in inst]
    pend=np.array([len(set(p)) for p in prods]); dmax=np.zeros(n,np.int64)
    tnew=np.full(n,-1,np.int64); held=np.zeros(n,np.int64)
    L=int(treal.max())*3+1000; I=np.full(L,FLOOR); In=np.full(L,FLOOR); inj=np.zeros(L)
    _res=[0]   # int ports reserved this cycle by overdue DMR checks (backpressure: late checks block commit)
    busy_ok=lambda q,w,m: (0 if (q==1 and w>=PORTW_INT-_res[0]) else m)   # busy units handled via pinning
    wake=defaultdict(list)
    for i in range(n):
        if pend[i]==0: wake[int(r[i])].append(i)
    # ROB-capacity gate: a uop may not become ready while any uop >= `rob` slots older is unissued
    # (real ROB frees at commit, after issue, so this is an optimistic/lower-bound constraint).
    uages=sorted(set(age)); cnt=Counter(age); ptr=0; gated=[]; ngate=0
    ready=[]; done=0; t=0; tok=float(Qb); drawn=0.0; ev=0; evfull=0; unc=0.0; inev=False; short=False
    nsel=0; nchg=0; nhold=0
    lastdone=0; lastprog=0
    dprev=0; dummyq=0.0; ndum=0; THF=TH if THF is None else THF
    predL=np.zeros(L) if lead_rec else None
    rq=0.0   # charge the bank harvested from falling edges (pC)
    SIG=0.0  # supply-state signal from the circuit (mA of surplus supply current; 0 = no co-sim)
    csum=[float(s.sum()) for s in sh]; nsf=0
    if mp_h is not None:
        # PREDICTIVE HOLDS (MPCH): identified supply kernel h (mV per mA per 1-cycle pulse); predicted die-voltage
        # deviation over mp_Hp cycles = conv(h, issued-current history + committed tails + candidates - dc)
        _h=np.asarray(mp_h,float); _K=len(_h); _HS=10
        _Hm=np.zeros((mp_Hp,_HS))
        for c in range(mp_Hp):
            for m in range(min(_HS,c+1)): _Hm[c,m]=_h[c-m]
        _dc=FLOOR; _md=MAXD if mp_maxd is None else mp_maxd; nmp=0
        def _shp(i):
            v=np.zeros(_HS); m=min(_HS,len(sh[i])); v[:m]=sh[i][:m]; return v
        def _vbase(persist=False):
            lo_=max(0,t-_K+1)
            fut=I[t:t+mp_Hp].copy()
            if persist:   # future = committed tails, floored at the recent REAL current level (persistence forecast)
                fut=np.maximum(fut,float((I[max(0,t-8):t]-Ifl[max(0,t-8):t]).mean()) if t>0 else FLOOR)
            x=np.concatenate([np.zeros(_K-1-(t-lo_)), I[lo_:t]-_dc, fut-_dc])
            return np.array([_h[::-1]@x[c:c+_K] for c in range(mp_Hp)])
    # ELASTIC WORK BANK: background work that must be done anyway (e.g. AES-CTR memory-encryption pads), consumed at
    # wb['rate'] ops/cycle from a buffer of wb['cap'] ops; issued on int ports the real picks leave free.
    # mode const: steady rate | greedy: any free slot while buffer not full | mpc: predictor places it where it flattens
    # the predicted supply voltage (fills overshoot/current dips, defers during predicted droop).
    # mode 'redund': opportunistic redundant re-execution of recently issued int ops (soft-error checking), no demand.
    wbs=None
    if wb is not None:
        wbs=dict(buf=wb.get('cap',64)/2, credit=0.0, starve=0.0, n=0, q=0.0, ring=[], cov=0, dq=[], late=0, maxq=0)
    Ifl=np.zeros(L)   # filler current (excluded from the predictor's DC reference -> no self-reinforcing fill)
    _rng=np.random.default_rng(12345)
    while done<n:
        if done!=lastdone: lastdone=done; lastprog=t
        if t-lastprog>5000 or t>=L-50: raise RuntimeError(f'stall at t={t} done={done}/{n} ready={len(ready)} gated={len(gated)}')
        if wb is not None and wb.get('mode','').startswith('dmr'):
            _res[0]=min(PORTW_INT,sum(1 for x in wbs['dq'] if x[0]<=t))
        if t in wake:
            for i in wake.pop(t):
                if ptr<len(uages) and age[i]-rob>=uages[ptr]: heapq.heappush(gated,(age[i],i)); ngate+=1
                else: ready.append(i)
        pick=[]; held_now=False
        if ready:
            ready.sort(key=lambda i:age[i])
            lo=In[max(0,t-W):t].min() if t>0 else FLOOR
            def need(sel):
                proj=I[t:t+H].copy()
                for i in sel: s=sh[i]; m=min(len(s),H); proj[:m]+=s[:m]
                return float(np.maximum(0,proj[:W]-(lo+TH)).sum())
            base=port_fill(ready,fu,iq,busy_ok); pick=base
            if pol in ('DITH','DITHD'):   # SPREAD-SPECTRUM issue: random heavy-first / light-first each cycle
                SIG=(sf_db+1.0)*(1.0 if _rng.random()<0.5 else -1.0)
            if pol in ('SF','SFD','DITH','DITHD') and base and abs(SIG)>sf_db:
                # SUPPLY-FOLLOWING reorder (work-conserving: same issue count as age order).
                # SIG>0: the supply is delivering a surplus (inductor current > load) -> issue the
                # heaviest ready ops now; SIG<0: deficit -> issue the lightest. Ops passed over MAXD
                # times are forced first (no starvation).
                forced=[i for i in ready if held[i]>=MAXD]
                order=forced+sorted([i for i in ready if held[i]<MAXD],key=lambda i:((-csum[i]) if SIG>0 else csum[i],age[i]))
                sel=[]
                for i in order:
                    if len(sel)>=len(base): break
                    if feasible(sel+[i],fu,iq,busy_ok): sel.append(i)
                if len(sel)==len(base):
                    if set(sel)!=set(base): nsf+=1
                    for i in base:
                        if i not in sel: held[i]+=1
                    pick=sel; base=sel
            if pol in ('SF','SFD','DITH','DITHD'):
                if pol in ('SFD','DITHD') and base and need(base)>tok+Ir:
                    forced=[i for i in pick if held[i]>=MAXD]
                    rest=sorted([i for i in pick if held[i]<MAXD],key=lambda i:-sh[i].sum())
                    kept=list(forced)
                    for i in reversed(rest):
                        if need(kept+[i])<=tok+Ir: kept.append(i)
                    if len(kept)<len(pick): nhold+=1; held_now=True
                    for i in pick:
                        if i not in kept: held[i]+=1
                    pick=kept
            elif pol in ('REG','REGH') and base and SIG<-sf_db:
                # INSTRUCTION-DOMAIN SHUNT REGULATOR: supply deficit -> current budget shrinks with |SIG|;
                # hold the costliest picks that exceed it (forced after MAXD as always)
                budget=max(0.0,reg_k*(SIG+60.0))
                forced=[i for i in pick if held[i]>=MAXD]
                rest=sorted([i for i in pick if held[i]<MAXD],key=lambda i:-sh[i].sum())
                kept=list(forced)
                for i in reversed(rest):
                    if need(kept+[i])<=budget: kept.append(i)
                if len(kept)<len(pick): nhold+=1; held_now=True
                for i in pick:
                    if i not in kept: held[i]+=1
                pick=kept
            elif pol not in ('A','REG','REGH') and base and need(base)>tok+Ir:
                nsel+=1
                sel=[]; pool=list(ready)
                while len(sel)<len(base):
                    best=None; bv=None
                    for i in pool:
                        if not feasible(sel+[i],fu,iq,busy_ok): continue
                        v=need(sel+[i])
                        if bv is None or v<bv-1e-9: best,bv=i,v
                    if best is None: break
                    sel.append(best); pool.remove(best)
                if len(sel)==len(base) and need(sel)<need(base)-1e-9: pick=sel
                if pol=='TD':
                    forced=[i for i in pick if held[i]>=MAXD]
                    rest=sorted([i for i in pick if held[i]<MAXD],key=lambda i:-sh[i].sum())
                    kept=list(forced)
                    for i in reversed(rest):                  # cheapest first
                        if need(kept+[i])<=tok+Ir: kept.append(i)
                    if len(kept)<len(pick): nhold+=1; held_now=True
                    pick=kept
                    for i in base:
                        if i not in pick: held[i]+=1
                if set(pick)!=set(base): nchg+=1
        if mp_h is not None and mp_vth is not None and pick:
            lo_=max(0,t-_K+1)
            x=np.concatenate([np.zeros(_K-1-(t-lo_)), I[lo_:t+mp_Hp]-_dc])
            vb=np.array([_h[::-1]@x[c:c+_K] for c in range(mp_Hp)])
            forced=[i for i in pick if held[i]>=_md]
            kept=forced+sorted([i for i in pick if held[i]<_md],key=lambda i:csum[i])    # light first
            def _J(sel):
                v=np.zeros(_HS)
                for i in sel: v+=_shp(i)
                return -float((vb+_Hm@v).min())
            while len(kept)>len(forced) and _J(kept)>mp_vth: kept.pop()                # drop heaviest
            if len(kept)<len(pick):
                nmp+=1; nhold+=1
                for i in pick:
                    if i not in kept: held[i]+=1
            pick=kept
        if mp_h is not None: _dc+=(I[t]-Ifl[t]-_dc)/2000.0
        for i in pick:
            ready.remove(i); tnew[i]=t; s=sh[i]; I[t:t+len(s)]+=s; done+=1; cnt[age[i]]-=1
            d=t-treal[i]
            for c in cons[i]:
                if d>dmax[c]: dmax[c]=d
                pend[c]-=1
                if pend[c]==0: wake[int(max(r[c]+max(0,dmax[c]),t+1))].append(c)
        if wbs is not None:
            free=PORTW_INT-sum(1 for i in pick if iq[i]==1); md=wb['mode']
            if md in ('dmr_eager','dmr_mpc'):
                # FULL DMR: every issued int op must be re-executed (checked) within wb['D'] cycles
                dq=wbs['dq']
                for i in pick:
                    if iq[i]==1: dq.append((t+wb['D'],i))
                nf=free; take=[]          # free already excludes the reserved ports' real use: reserved ports serve overdue checks
                if md=='dmr_eager':
                    take=dq[:nf]
                else:
                    due=[x for x in dq if x[0]<=t+1][:nf]
                    rest=[x for x in dq if x[0]>t+1]
                    v=_vbase(True); add=np.zeros(_HS)
                    for x in due: add+=_shp(x[1])
                    best=(float(((v+_Hm@add)**2).sum()),0); a2=add.copy()
                    for kk in range(1,min(nf-len(due),len(rest))+1):
                        a2=a2+_shp(rest[kk-1][1]); c_=float(((v+_Hm@a2)**2).sum())
                        if c_<best[0]-1e-9: best=(c_,kk)
                    take=due+rest[:best[1]]
                for x in take:
                    s_=sh[x[1]]; I[t:t+len(s_)]+=s_; wbs['q']+=float(s_.sum()); wbs['cov']+=1   # steady load: stays in the DC reference
                ids=set(id(x) for x in take); wbs['dq']=dq=[x for x in dq if id(x) not in ids]
                wbs['late']+=sum(1 for x in dq if x[0]<=t); wbs['maxq']=max(wbs['maxq'],len(dq))
            elif md=='redund':
                for i in pick:
                    if iq[i]==1: wbs['ring'].append(i)
                wbs['ring']=wbs['ring'][-64:]
            if md not in ('redund','dmr_eager','dmr_mpc'):
                wbs['buf']-=wb['rate']
                if wbs['buf']<0: wbs['starve']+=-wbs['buf']; wbs['buf']=0.0
            room=wb.get('cap',64)-wbs['buf'] if md not in ('redund','dmr_eager','dmr_mpc') else 1e9
            nmax=int(max(0,min(free,room)))
            nwb=0
            if md=='const':
                wbs['credit']+=wb['rate']; nwb=int(min(nmax,np.floor(wbs['credit']))); wbs['credit']-=nwb
            elif md=='greedy':
                nwb=nmax
            elif md=='mpc' and nmax>0:
                v=_vbase(); WS=np.zeros(_HS); WS[:4]=DSH
                need=min(nmax,2) if wbs['buf']<wb.get('low',8) else 0
                best=None
                for nn in range(need,nmax+1):
                    vv=v+_Hm@(nn*WS)                  # h<0: added current lowers the predicted voltage
                    cst=float((vv**2).sum())
                    if best is None or cst<best[0]-1e-9: best=(cst,nn)
                nwb=best[1]
            elif md=='redund' and nmax>0 and wbs['ring']:
                # OVERSHOOT-CANCELLING redundant re-execution: act only if the predicted overshoot exceeds vo; re-execute
                # the FEWEST recent int ops that bring it back under vo, never pushing the predicted dip past the droop target
                v=_vbase(True); vo=wb['vo']; vd=wb.get('vd', mp_vth if mp_vth is not None else 1e9)
                if v.max()>vo:
                    add=np.zeros(_HS); k=0
                    ring=wbs['ring']
                    while k<min(nmax,len(ring)):
                        a2=add+_shp(ring[-1-k]); vv=v+_Hm@a2
                        if -vv.min()>vd: break
                        add=a2; k+=1
                        if vv.max()<=vo: break
                    if k:
                        sel=ring[-k:]; del ring[-k:]
                        for i in sel: s_=sh[i]; I[t:t+len(s_)]+=s_; Ifl[t:t+len(s_)]+=s_; wbs['q']+=float(s_.sum())
                        wbs['cov']+=k
            if nwb:
                I[t:t+4]+=nwb*DSH; wbs['buf']+=nwb; wbs['n']+=nwb; wbs['q']+=nwb*float(DSH.sum())
        if pol=='REG' and SIG>reg_sth:    # supply surplus -> shunt it through dummy ALU ops on free int ports
            free=PORTW_INT-sum(1 for i in pick if iq[i]==1)
            nd=max(0,min(free,int(SIG/DSH[1])))
            if nd: I[t:t+4]+=nd*DSH; dummyq+=nd*DSH.sum(); ndum+=nd
        if dummy:
            # dummy ALU ops on int ports the real picks left free this cycle
            free=PORTW_INT-sum(1 for i in pick if iq[i]==1); nd=0
            if 'ramp' in dummy:          # ramp up <=1/cycle while real ops are held, down <=1/cycle after
                want=free if held_now else 0
                nd=min(dprev+1,want) if want>dprev else max(dprev-1,want,0)
            imminent=True
            if 'bridge' in dummy:        # only bridge GAPS: real work pending now or waking within BK cycles
                imminent=bool(ready) or bool(gated) or any((t+j) in wake for j in range(1,BK+1))
            if ('fill' in dummy or 'bridge' in dummy) and t>0 and imminent:  # falling-edge fill: I[t+1] may not fall more than THF below recent max
                tgt=In[max(0,t-W):t].max()-THF; gap=tgt-I[t+1]
                if gap>0: nd=max(nd,int(np.ceil(gap/DSH[1])))
            nd=max(0,min(nd,free))
            if nd: I[t:t+4]+=nd*DSH; dummyq+=nd*DSH.sum(); ndum+=nd
            dprev=nd
        if lead_rec: predL[t]=I[t+lead_rec]-I[t]   # CAUSAL: change the issue stage already knows is coming (ops issued <=t)
        while ptr<len(uages) and cnt[uages[ptr]]==0: ptr+=1
        while gated and (ptr>=len(uages) or gated[0][0]-rob<uages[ptr]):
            wake[t+1].append(heapq.heappop(gated)[1])
        tok=min(float(Qb),tok+Ir)
        lo=In[max(0,t-W):t].min() if t>0 else FLOOR
        ex=max(0.0,I[t]-(lo+TH)); b=0.0
        if ex>0:
            b=min(ex,tok); tok-=b; drawn+=b; inj[t]=b; unc+=ex-b
            if not inev: ev+=1; inev=True; short=False
            if b<ex-1e-9: short=True
        elif inev:
            evfull+= (not short); inev=False
        In[t]=I[t]-b
        if recyc and t>0 and tok<Qb:
            # EDGE-HARVESTING bank: on a falling edge the bank SINKS current from the grid (like a dummy
            # load: same smoothing) but pumps it into its reservoir at efficiency recyc instead of burning it.
            gap=(In[max(0,t-W):t].max()-THF)-In[t]
            if gap>0:
                sk=min(gap,(Qb-tok)/recyc); tok+=recyc*sk; In[t]+=sk; rq+=sk
        if vdd_recharge and tok<Qb: In[t]+=min(Ir,Qb-tok)
        t+=1
    end=t+40
    return dict(tnew=tnew,treal=treal,I=I[:end],In=In[:end],makespan=int(tnew.max()),
                events=ev,evfull=evfull,unc=unc,drawn=drawn,nsel=nsel,nchg=nchg,nhold=nhold,ngate=ngate,dummyq=dummyq,ndum=ndum,predL=(predL[:end] if lead_rec else None),rq=rq,nsf=nsf,wb=(None if wbs is None else {k:v for k,v in wbs.items() if k not in ('ring','dq')}))

def calibrate(inst, prods, cons, age, maxit=30):
    """pin uops the port model would issue EARLIER than hardware (unobservable busy/pause
    stalls) to their real time; iterate to a fixed point. Returns pin mask."""
    pin=np.zeros(len(inst),bool)
    for it in range(maxit):
        R=run(inst,prods,cons,age,pin,'A',1e9,0,0)
        d=R['tnew']-R['treal']; early=(d<0)&~pin
        if not early.any() and not (d!=0).any(): return pin,it
        pin|=early
        if not early.any(): break
    return pin,it

def droop(Iamp_mA, L=507e-12, C=5e-9, R=79.6e-3, dt=1e-9, sub=4):
    """series-RLC PDN, exact ZOH discretisation on sub-steps; returns droop trace (mV) vs DC-at-floor."""
    A=np.array([[-R/L,-1/L],[1/C,0.0]]); B=np.array([0.0,-1/C]); h=dt/sub
    M=np.zeros((3,3)); M[:2,:2]=A*h; M[:2,2]=B*h; E=expm(M); Ad=E[:2,:2]; Bd=E[:2,2]
    I=np.asarray(Iamp_mA)*1e-3; i0=I[0]
    x=np.array([i0,-R*i0]); out=np.empty(len(I))    # state: inductor current, die-voltage deviation
    for k in range(len(I)):
        u=I[k]; mn=0.0
        for _ in range(sub):
            x=Ad@x+Bd*u; mn=min(mn,x[1])
        out[k]=mn
    return -out*1e3

if __name__=='__main__':
    wl=sys.argv[1]
    inst,gap=parse(f'traces/pm2_{wl}.log.gz')
    prods,cons=build(inst); age=unwrap_age(inst)
    pin,it=calibrate(inst,prods,cons,age)
    R=run(inst,prods,cons,age,pin,'A',1e9,0,0)
    d=R['tnew']-R['treal']
    print(f"{wl}: {len(inst)} uops; pinned {pin.sum()} ({100*pin.mean():.2f}%) after {it} iters; "
          f"check1 A==HW: moved {int((d!=0).sum())}, max|d| {int(np.abs(d).max())}, makespan {R['makespan']} vs {int(R['treal'].max())}")
    from collections import Counter
    print("  pinned by fu_code:",Counter(inst[i][2] for i in np.nonzero(pin)[0]).most_common(6))
