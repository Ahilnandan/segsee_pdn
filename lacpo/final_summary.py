"""Final single-core summary from eval6_FULL.json: performance at equal reliability.
Guardband needed: no replay -> worst max droop over the 10 benches (every droop must be covered); with selective replay
-> the test guardband G, valid only if no error went undetected. Reference = all-decap, no policy, no replay.
Frequency view: f ~ 1 / d(V) at V = VNOM - guardband (timing6 alpha-power law), divided by (1 + mean slowdown).
Voltage view: same frequency, supply lowered by the guardband saved -> dynamic power ~ V^2 (bank losses added back)."""
import json, collections, numpy as np
import timing6 as TM
VNOM = 1.1
R = [r for r in json.load(open('eval6_FULL.json')) if not r.get('collapse')]
g = collections.OrderedDict()
for r in R: g.setdefault((r['design'], r['tag']), []).append(r)
ref = max(r['max'] for r in g[('decap', 'A')])
f = lambda G: TM.d(VNOM - ref * 1e-3) / TM.d(VNOM - G * 1e-3)
print('reference guardband (all-decap, no policy, worst max over 10 benches): %.1f mV' % ref)
print('%-10s %-18s %7s %7s %8s %8s %8s %7s %9s %9s %9s' % ('design', 'config', 'mp999', 'wmax', 'slow%', 'slowmax', 'err/1M', 'undet', 'guard mV', 'perf %', 'power %'))
rows = []
for (dn, tag), rs in g.items():
    n = sum(r['nuops'] for r in rs); undet = sum(r['nundet'] for r in rs); err = 1e6 * sum(r['nerr'] for r in rs) / n
    sl = np.mean([r['slowdown'] for r in rs]); slw = max(r['slowdown'] for r in rs); wmax = max(r['max'] for r in rs); loss = np.mean([r['energy_over'] for r in rs])
    G = rs[0]['G'] if rs[0]['G'] is not None else wmax
    ok = undet == 0
    perf = 100 * (f(G) / (1 + sl) - 1) if ok else float('nan')
    Vop = VNOM - (ref - G) * 1e-3; power = 100 * ((Vop / VNOM) ** 2 * (1 + loss) - 1) if ok else float('nan')
    rows.append(dict(design=dn, config=tag, mp999=float(np.mean([r['p999'] for r in rs])), wmax=wmax, slow=sl, slowmax=slw, err=err, undet=undet, guard=G, perf=perf, power=power, loss=loss))
    print('%-10s %-18s %7.1f %7.1f %8.2f %8.2f %8.0f %7d %9.1f %+9.1f %+9.1f' % (dn, tag, rows[-1]['mp999'], wmax, 100 * sl, 100 * slw, err, undet, G, perf, power))
json.dump(dict(ref_guard=ref, rows=rows), open('final_summary.json', 'w'), indent=1)
