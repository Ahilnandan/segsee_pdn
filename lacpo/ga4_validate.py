"""Validate ga4 fronts: all 10 FULL traces, jitter/guardband (qsort 100k cycles), DVFS step per design.
Passive reference = smallest uniform decap meeting both limits on the full traces (pdn4)."""
import os, sys, json, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn4, ga4, jitter3
HERE = ga4.HERE

def front(arch, kmax=6):
    r = json.load(open(HERE + 'ga4_%s.json' % arch)); X = np.array(r['X']); F = r['F']
    feas = [i for i, f in enumerate(F) if f[2] == 0]
    nd = [i for i in feas if not any(F[j][0] <= F[i][0] and F[j][1] <= F[i][1] and (F[j][0] < F[i][0] or F[j][1] < F[i][1]) for j in feas)]
    nd = sorted(set(nd), key=lambda i: F[i][0]); u = []
    for i in nd:
        if not u or abs(F[i][0] - F[u[-1]][0]) > 2e-3 or abs(F[i][1] - F[u[-1]][1]) > 2e-3: u.append(i)
    sel = u if len(u) <= kmax else [u[int(k)] for k in np.linspace(0, len(u) - 1, kmax)]
    return [(arch, X[i].tolist(), F[i]) for i in sel]

def dvfs(arch, d, p, spec):
    S1, S2, N = 3000, 13000, 23000
    z = np.load(ga4.B + 'qsort_rep_%s.npz' % d['var']); I = z['I'][60000:60000 + N].astype(np.float32)
    vref = np.full(N, 1.0); vref[S1:S2] = 1.1; kw = {}
    if arch == 'B':
        sw = np.ones(N, np.uint8); kw['role_off'] = sw
    if arch == 'D':
        vpre = np.full(N, 1.0); vpre[S1 - 2000:S2] = 1.1; kw['vref_off'] = vpre
    m = pdn4.simulate(I, p, spec, vref=vref, want_v=True, warm=4000, **kw); v = np.convolve(m['vavg'].mean(1), np.ones(100) / 100, 'same')
    def st(k0, k1):
        t = vref[k0 + 5]; ok = np.abs(v[k0:k1] - t) < 0.01 * t; need = int(1e-6 / pdn4.T_CYC)
        for k in range(len(ok) - need):
            if ok[k:k + need].all(): return k * pdn4.T_CYC * 1e9
        return float('nan')
    return st(S1, S2), st(S2, N)

def job(item):
    arch, g, F = item
    if arch == 'P':
        Ax = g; d = dict(arch='P', var='A', th_p=6, lead=3, pred_on=False, boost=0.0); p = dict(A_x=[Ax * s for s in pdn4.SHARE]); spec = None; ro = None
        per = {}
        for w in ga4.WL:
            I, Ih = ga4.trace(w, 'A'); s = ga4.STARTS[w]; m = pdn4.simulate(I[s:], p); per[w] = dict(p999=m['p999'], max=m['max'], o=0.0)
        area = Ax
    else:
        a, o, cv, info = ga4.evaluate(arch, g, full=True); d = ga4.decode(arch, g); p, spec, ro = ga4.build_design(d); per = info['per']; area = a
    I, Ih = ga4.trace('qsort', d['var']); s = ga4.STARTS['qsort']
    m = ga4.sim(d, p, spec, ro, I[s + 60000:s + 160000], Ih[s + 60000:s + 160000]) if arch != 'P' else pdn4.simulate(I[s + 60000:s + 160000], p, want_v=True)
    if 'vmin' not in m:
        pred = pdn4.predict_triggers(Ih[s + 60000:s + 160000], d['th_p'], d['lead']) if d['pred_on'] else None
        m = pdn4.simulate(I[s + 60000:s + 160000], p, spec, pred=pred, role_off=(np.ones(100000, np.uint8) if ro == 'helper' else None), want_v=True)
    j = jitter3.jitter(m['vavg'], m['vmin'])
    up, dn = dvfs(arch, d, p, spec) if arch != 'P' else dvfs('A', d, p, spec)
    return dict(arch=arch, area=area, ga_F=F, worst_p999=max(v['p999'] for v in per.values()), worst_max=max(v['max'] for v in per.values()),
                overhead=float(np.mean([v['o'] for v in per.values()])), per=per, jitter_p999=j['Jp_p999'], guardband=j['guardband_p9999'],
                dvfs_up_ns=up, dvfs_dn_ns=dn, design={k: v for k, v in d.items()})

if __name__ == '__main__':
    import multiprocessing as mp
    items = [('P', a, None) for a in (0.6, 0.7, 0.8, 0.9)]
    for arch in 'ABDEFG':
        if os.path.exists(HERE + 'ga4_%s.json' % arch): items += front(arch)
    with mp.get_context('spawn').Pool(10) as P: res = P.map(job, items)
    json.dump(res, open(HERE + 'ga4_validate.json', 'w'), indent=1, default=float)
    ok = lambda r: r['worst_p999'] <= 55 and r['worst_max'] <= 80
    ref = min([r['area'] for r in res if r['arch'] == 'P' and ok(r)] or [float('nan')])
    print('passive reference (pdn4, full traces): %.2f mm2' % ref)
    for r in res:
        d = r['design']
        print('%s area %.3f (%+5.1f %% die vs ref) overhead %5.2f %% | p99.9 %5.1f max %5.1f %s | jitter %4.1f ps gb %4.2f %% | DVFS up %6.0f dn %6.0f ns | var %s' % (
            r['arch'], r['area'], 100 * (ref - r['area']) / pdn4.DIE_MM2, 100 * r['overhead'], r['worst_p999'], r['worst_max'], 'OK ' if ok(r) else 'MISS',
            r['jitter_p999'], r['guardband'], r['dvfs_up_ns'], r['dvfs_dn_ns'], d.get('var')))
