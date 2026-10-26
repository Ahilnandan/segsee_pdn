"""APPROXIMATE per-uop current model, rebuilt 2026-09-30 after the Tier-A NNLS library
(per_instruction_library.json) was lost with /tmp. Energies/peak timing follow the
Tier-A table (ADD~45pJ 1-3cyc; MUL~115pJ peak@+2; LD~64pJ sharp@+0; CSR/serialising
~260pJ spike@+4; FP~40pJ). sum(mA*cycle) = E(pJ)/0.75V. Absolute values +-30-40%;
use for relative/event-rate questions, not absolute droop claims."""
import re, gzip, numpy as np
FLOOR=27.6
S = {
 'alu':   [10,20,20,10],            'logic':[15,35,30,13],   'shift':[20,40,35,25],
 'mul':   [15,30,55,35,18],         'div':  [20]*10,
 'ld':    [63,15,7],                'sta':  [57,25,10],      'std':  [48,20,8],
 'br':    [22,45,25,15],            'jmp':  [20,30,15,8],
 'ser':   [0,0,0,5,330,12],         # CSR/fence/eret/sfence/wfi: serialising drain spike
 'fp':    [10,20,15,8],             'fmv':  [24,25,12,6],    'fdiv': [15]*10,
 'nop':   [0],
}
def cls(u):
    if u==0: return 'nop'
    if u==1: return 'ld'
    if u==2 or u==67: return 'sta'
    if u==3: return 'std'
    if u in (6,7,8,19,20,21): return 'logic'
    if u in (11,12,13,14,22,23,46,47,48,49,50,51): return 'shift'
    if 24<=u<=29: return 'br'
    if 30<=u<=35 or u in (40,41,42,65,105,106,107): return 'ser'
    if u in (36,37,38): return 'jmp'
    if 52<=u<=56: return 'mul'
    if 57<=u<=64: return 'div'
    if 68<=u<=79: return 'fmv'
    if u in (101,102,103,104): return 'fdiv'
    if 72<=u<=100: return 'fp'
    return 'alu'
SHAPE={u:np.array(S[cls(u)],float) for u in range(0,110)}
COST={u:float(s.sum()) for u,s in SHAPE.items()}
def port(u):
    c=cls(u)
    if c in ('ld','sta','std'): return 'mem'
    if c in ('fp','fmv','fdiv'): return 'fp'
    return 'int'
PORTW={'int':4,'mem':2,'fp':2}      # PaperMatchBoomV3Config real issue widths (NOT MediumBoom 2/1/1)

def parse_iss(path):
    rx=re.compile(r'BACT_ISS\s+(\d+)\s+\d+\s+(\d+)\s+\d+'); ev=[]
    with gzip.open(path,'rt',errors='ignore') as f:
        for ln in f:
            if ln.startswith('BACT_ISS'):
                m=rx.match(ln)
                if m: ev.append((int(m.group(1)),int(m.group(2))))
    ev.sort()
    # HTIF load-phase cut: drop everything before the longest issue-free gap (memory note)
    c=np.array([e[0] for e in ev]); g=np.diff(c); i=int(np.argmax(g))
    return ev[i+1:], int(g[i])

def current_from_times(ev, n=None):
    t0=ev[0][0]; n=n or (ev[-1][0]-t0+20)
    I=np.full(n,FLOOR)
    for c,u in ev:
        s=SHAPE[u]; k=c-t0
        if k+len(s)<=n: I[k:k+len(s)]+=s
    return I
