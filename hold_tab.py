import json
R=[json.loads(l) for l in open('hold_exp.jsonl')]
key=lambda r:(r['mix'],r['frac'],str(r['hold']),r['par'])
R.sort(key=key)
base={(r['mix'],r['frac']):r for r in R if r['hold'] is None}
print(f"{'mix':9s} {'Qb':>3s} {'hold':5s} {'par':>5s} | {'slow%':>6s} {'#held':>6s} | {'p99.9':>6s} {'d%':>6s} {'worst':>6s} {'d%':>6s} {'top20':>6s} {'d%':>6s}")
for r in R:
    b=base.get((r['mix'],r['frac']))
    if b is None: continue
    pc=lambda k:100*(r[k]/b[k]-1)
    print(f"{r['mix']:9s} {r['frac']:3.1f} {str(r['hold']):5s} {r['par']:5.2f} | {r['slow']:6.2f} {r['nheld']:6d} | {r['p999']:6.1f} {pc('p999'):6.1f} {r['worst']:6.1f} {pc('worst'):6.1f} {r['top20']:6.1f} {pc('top20'):6.1f}")
