"""mc_full: multi-core coordinator x per-core PDN design (decap / banks / banks_sep) with droop, slowdown, bank loss and
post-hoc timing-error estimate (timing6.errors per core, that core's worst-tile V per cycle, G = 30/20/10 mV).
Post-hoc errors: selective replay costs ~0 slowdown on one core at rob=80, so schedule/voltage are taken as unperturbed.
  python mc_full.py [workers<=4]  -> mc_full.json + table"""
from __future__ import print_function
import os, sys, json, collections, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'; sys.path.insert(0, HERE)
import mc_engine as M, mc_eval as E, pdn7, timing6 as TM
GS = (30, 20, 10)


def run_full(engs, starts, design1, W=60000, coord=None, rob=80, warm=600, q_mV=5.0, sensor_delay=2):
    nc = len(engs); cores = [M.Core(engs[k], starts[k], starts[k] + W, rob) for k in range(nc)]
    d4 = pdn7.place_design(design1); has_bank = bool(d4.get('bank'))
    sim = pdn7.Sim(d4, Im_A=pdn7.place([c.BG[:3000].mean(0) for c in cores]) * 1e-3)
    prm = dict(M.CDEF, **(coord or {})); on = coord is not None and prm['budget'] is not None
    st = dict(ew=[0.0] * nc, lag=[0] * nc); S_hist = []; vq = collections.deque([pdn7.VNOM] * sensor_delay, maxlen=sensor_delay)
    TI = [pdn7.tiles_of(k) for k in range(nc)]; vc = []; arm = np.zeros((1, pdn7.T), np.uint8); vref = np.array([pdn7.VNOM]); t = 0
    while any(c.active(t) for c in cores):
        cands = [c.cand(t) for c in cores]
        pk = M.coordinate(cores, cands, t, S_hist, st, prm, vq[0]) if (on and any(cands)) else [(cd, []) for cd in cands]
        It = [c.issue(t, p, h) for c, (p, h) in zip(cores, pk)]; S_hist.append(sum(c.Ihat[t] for c in cores))
        vmin, _, _ = sim.run(pdn7.place(It)[None, :].astype(np.float32), vref=vref, arm=arm, want_res=has_bank)
        vk = [float(vmin[0][ti].min()) for ti in TI]; vc.append(vk); v = min(vk)
        vq.append(np.round(v / (q_mV * 1e-3)) * q_mV * 1e-3); t += 1
    tend = min(c.tend for c in cores); VC = np.array(vc); dr = (pdn7.VNOM - VC[warm:tend].min(1)) * 1e3
    s = sim.stats; loss = float((s[2] + s[3] + s[5]) / s[0]) if s[0] else float('nan')
    err = {}
    for G in GS:
        n = det = und = 0.0
        for k, c in enumerate(cores):
            m = c.tnew >= warm
            e = TM.errors(c.tnew[m], c.fu[m], VC[:, k], G, w=TM.W)
            n += sum(x['n'] for x in e.values()); det += sum(x['det'] for x in e.values()); und += sum(x['undet'] for x in e.values())
        err[G] = dict(per1M=(det + und) / n * 1e6, undet=und, det=det, n=n)
    return dict(p999=float(np.percentile(dr, 99.9)), max=float(dr.max()), sd=[c.slowdown() for c in cores],
                exec=[int(c.tnew.max()) for c in cores], loss=loss, stats=list(map(float, s)), err=err)


def job(a):
    dn, sc, tag, coord = a
    import policy6 as P6, eval6
    ws, offs = E.SCEN[sc]; engs = {w: P6.Engine(w) for w in set(ws)}; D, _ = eval6.designs()
    r = run_full([engs[w] for w in ws], [E.A0 + int(o) for o in offs], D[dn], W=E.W, coord=coord)
    r.update(design=dn, scen=sc, tag=tag); print(json.dumps(r), flush=True); return r


if __name__ == '__main__':
    import multiprocessing as mp
    nw = min(4, int(sys.argv[1]) if len(sys.argv) > 1 else 4)
    jobs = [(dn, sc, tag, c) for dn in ('decap', 'banks', 'banks_sep') for sc in ('lock', 'jit', 'diff', 'rand')
            for tag, c in (('none', None), ('B30', dict(budget=30.0)))]
    with mp.get_context('spawn').Pool(nw) as P: res = P.map(job, jobs, chunksize=1)
    json.dump(res, open(HERE + 'mc_full.json', 'w'), indent=1)
    base = {(r['design'], r['scen']): r for r in res if r['tag'] == 'none'}
    print('%-9s %-4s %-4s %6s %6s %6s %13s ' % ('design', 'scen', 'conf', 'p99.9', 'max', 'loss', 'sd mean/worst') +
          ' '.join('G%d e/1M undet' % G for G in GS) + '  Gmin')
    for r in res:
        b = base[(r['design'], r['scen'])]; rel = [(e / max(1, e0) - 1) * 100 for e, e0 in zip(r['exec'], b['exec'])]
        sd = rel if r['tag'] != 'none' else [x * 100 for x in r['sd']]
        ok = [G for G in GS if r['err'][G]['undet'] == 0 and r['err'][G]['per1M'] < 2000]
        print('%-9s %-4s %-4s %6.1f %6.1f %6.3f %6.2f/%6.2f ' % (r['design'], r['scen'], r['tag'], r['p999'], r['max'], r['loss'], np.mean(sd), max(sd)) +
              ' '.join('%8.0f %5.1f' % (r['err'][G]['per1M'], r['err'][G]['undet']) for G in GS) + '  %s' % (min(ok) if ok else '>30'))
