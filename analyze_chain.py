"""Tables + iso-droop deliverables from chain_exp.jsonl.
iso-fsw : lowest FIVR switching frequency at which a policy's worst droop <= the reference
          (raw hardware current, FIVR 140 MHz, 5 nF, no feed-forward). Linear interpolation
          between simulated fsw points; '<35' means even the lowest point beats the reference.
iso-C   : least on-die decap (of 2.5/3.5/5 nF) at which droop <= reference, per topology.
metric  : 'droop' = worst absolute droop below 0.9 V (what the guardband must cover).
"""
import json, sys, numpy as np
R = [json.loads(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else 'chain_exp.jsonl')]
M = sys.argv[2] if len(sys.argv) > 2 else 'droop'
POLS = ['raw', 'hold+bank', 'dummy+bank', 'combined150', 'combined200']
def get(wl, pol, topo, fsw, ff, c=5):
    x = [r for r in R if r['wl'] == wl and r['pol'] == pol and r['topo'] == topo and r['fsw'] == fsw and r['ff'] == ff and r['cdie'] == c]
    return x[0] if x else None
for wl in ('mm', 'qsort', 'dhrystone'):
    print(f'\n=== {wl}  ({M}, mV; slowdown; dummy energy) ===')
    print(f'{"policy":12s} {"slow%":>6s} {"dumE%":>6s} | {"offchip":>8s} {"FIVR140":>8s} {"F140+ff":>8s} | over(off/F140) | iso-fsw (no ff / ff) | iso-C off/F140/F140ff')
    ref = get(wl, 'raw', 'B', 140, False)
    for pol in POLS:
        a = get(wl, pol, 'A', 0, False); b = get(wl, pol, 'B', 140, False); bf = get(wl, pol, 'B', 140, True)
        if not (a and b and bf): print(f'{pol:12s} (incomplete)'); continue
        iso = []
        for ff in (False, True):
            pts = sorted([(r['fsw'], r[M]) for r in R if r['wl'] == wl and r['pol'] == pol and r['topo'] == 'B' and r['ff'] == ff and r['cdie'] == 5])
            ok = [f for f, d in pts if d <= ref[M]]
            if not ok: iso.append('>200')
            else:
                f0 = min(ok); i = [f for f, _ in pts].index(f0)
                if i == 0: iso.append('<=%d' % f0)
                else:
                    (fa, da), (fb, db) = pts[i - 1], pts[i]
                    iso.append('%.0f' % (fa + (ref[M] - da) * (fb - fa) / (db - da)) if db != da else str(fb))
        isoc = []
        for topo, fsw, ff in (('A', 0, False), ('B', 140, False), ('B', 140, True)):
            refc = get(wl, 'raw', topo, fsw, ff) if topo == 'A' else ref
            cs = [c for c in (2.5, 3.5, 5) if (g := get(wl, pol, topo, fsw, ff, c)) and g[M] <= refc[M]]
            isoc.append(('%.1f' % min(cs)) if cs else '>5')
        print(f'{pol:12s} {a["slow"]:6.2f} {a["dummyE"]:6.1f} | {a[M]:8.1f} {b[M]:8.1f} {bf[M]:8.1f} | {a["over"]:5.0f}/{b["over"]:5.0f}    | {iso[0]:>6s} / {iso[1]:<6s}      | {"/".join(isoc)}')
    print(f'  reference = raw, FIVR 140 MHz, 5 nF: {ref[M]:.1f} mV;  raw offchip 5 nF: {get(wl,"raw","A",0,False)[M]:.1f} mV')
