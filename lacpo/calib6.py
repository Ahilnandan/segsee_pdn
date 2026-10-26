"""calib6: why the B4 scheduler-v2 terms (window, charge debt, 77 MHz tone, bank readiness) never acted in eval6 B8, and a
calibration so each one does. policy6.py is NOT modified: IEngine re-compiles policy6.Engine.run from its own source with
logging lines spliced in (same arithmetic, so a logged run is bit-identical to the plain one).
  inst  : full traces, banks design, TOK base + rob=80. Per cycle (when the ready set is non-empty and TOK decides):
          s_short - s_long (mA), debt (nC), tone magnitude (mA), bank_frac, forecast min fb (mV), the base threshold
          th0 (no extra terms), len(ready) and len(base pick). Saved to calib6_inst.npz.
  calib : thresholds/weights from the pooled distributions (calib6_calib.json); also prints, for the B8 guesses, how
          often each term made extra > 0 and how often it flipped viol().
  eval  : each calibrated term alone, all four together (no refill arming), TOK base; hold=0 and hold=1; 10 benches.
Usage: python calib6.py <inst|calib|eval|report> [workers<=4]"""
from __future__ import print_function
import os, sys, json, inspect, textwrap, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
BENCHES = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
B8 = dict(win=dict(win_dI=10.0, win_k=0.3), debt=dict(debt_k=0.02), tone=dict(tone_thr=3.0, tone_k=1.0),
          ready=dict(ready_min=0.5, ready_k=8.0))

def make_engine():
    import policy6 as P6
    src = textwrap.dedent(inspect.getsource(P6.Engine.run))
    hook_sig = "            if ready and t >= freeze:\n"
    hook_dec = "                    fb = forecast([(i, 0) for i in base])\n"
    hook_sig = hook_sig[4:]; hook_dec = hook_dec[4:]      # textwrap.dedent removed 4 spaces of class indentation
    assert src.count(hook_sig) == 1 and src.count(hook_dec) == 1
    sig = ("        _tn = abs(TONE.dot(Ihat[max(0, t - 24):t + 4][-28:]) if t >= 24 else 0.0) / 28.0\n"
           "        _sig = (t, s_short - s_long, sq[0], _tn, bank_frac, margin, len(ready))\n")
    dec = ("                _th0 = float(np.clip(prm['th'] + prm['kv'] * np.clip(margin, -20, 20) + prm['kb'] * (bank_frac if has_bank else 0.0)"
           " - (prm['fb_err'] * err_ewma * 1000.0 if prm['fb_err'] else 0.0), 2.0, 40.0))\n"
           "                self._LOG.append(_sig + (float(fb.min()), _th0, th_eff, len(base)))\n")
    src = src.replace(hook_sig, sig + hook_sig).replace(hook_dec, hook_dec + dec)
    ns = dict(P6.__dict__); exec(compile(src, '<calib6-run>', 'exec'), ns)
    class IEngine(P6.Engine):
        def run_log(self, *a, **k):
            self._LOG = []; r = ns['run'](self, *a, **k)
            r['log'] = np.array(self._LOG, float); return r
    return IEngine

COLS = ['t', 'win', 'debt', 'tone', 'bank', 'margin', 'nready', 'fbmin', 'th0', 'th_eff', 'nbase']
def setup(extra):
    import eval6
    D, prmb = eval6.designs(); prm = dict(prmb, rob=80, **extra); return D, prm

def job_inst(bench):
    E = make_engine()(bench); D, prm = setup({})
    r = E.run_log(D['banks'], 'TOK', prm)
    np.save(HERE + 'calib6_inst_%s.npy' % bench, r['log'])
    return bench, {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}

def contrib(L, name, p):
    """extra (mV) that term `name` with params p adds on each logged decision."""
    c = {k: L[:, i] for i, k in enumerate(COLS)}
    if name == 'win': return p['win_k'] * np.maximum(0, c['win'] - p['win_dI'])
    if name == 'debt': return p['debt_k'] * np.maximum(0, c['debt'])
    if name == 'tone': return p['tone_k'] * np.maximum(0, c['tone'] - p['tone_thr'])
    if name == 'ready': return p['ready_k'] * np.maximum(0, p['ready_min'] - c['bank'])
def flips(L, ex):
    """decisions where viol (forecast part) changes: fbmin >= -th0 but < -clip(th0 - extra) (edge term not used here)."""
    th0 = L[:, COLS.index('th0')]; fb = L[:, COLS.index('fbmin')]
    return (fb >= -th0) & (fb < -np.clip(th0 - ex, 2.0, 40.0))

def load_logs():
    return {b: np.load(HERE + 'calib6_inst_%s.npy' % b) for b in BENCHES}

TARGET = 0.02      # each term active on the top ~2 % of decisions (pooled over benches)
def calibrate(logs, verbose=True):
    P = np.concatenate([logs[b] for b in BENCHES]); c = {k: P[:, i] for i, k in enumerate(COLS)}
    slack = c['th0'] + c['fbmin']               # mV; > 0 = no violation at base threshold
    q = lambda x, p: float(np.percentile(x, p))
    if verbose:
        print('pooled decisions: %d' % len(P))
        for k in ('win', 'debt', 'tone', 'bank', 'margin', 'fbmin', 'th0', 'nready', 'nbase'):
            x = c[k]; print('  %-6s p1 %8.3f p50 %8.3f p95 %8.3f p98 %8.3f p99.5 %8.3f max %8.3f' % (k, q(x, 1), q(x, 50), q(x, 95), q(x, 98), q(x, 99.5), x.max()))
        print('  slack  p1 %.2f p5 %.2f p50 %.2f ; base viol frac %.4f ; choice (nready>nbase) frac %.3f' % (q(slack, 1), q(slack, 5), q(slack, 50), (slack < 0).mean(), (c['nready'] > c['nbase']).mean()))
        print('  frac th0 at its 18.3 mV ceiling %.3f ; margin>=20 (kv term saturated) %.3f ; base viol with choice %.4f' % ((c['th0'] > 18.3).mean(), (c['margin'] >= 20).mean(), ((slack < 0) & (c['nready'] > c['nbase'])).mean()))
        print('B8 guesses: frac extra>0 / frac flips')
        for nm, p in B8.items():
            ex = contrib(P, nm, p); fl = flips(P, ex)
            print('  %-6s %.5f / %.5f  (extra max %.3f mV; flips with a choice nready>nbase %.5f)' % (nm, (ex > 0).mean(), fl.mean(), ex.max(), (fl & (c['nready'] > c['nbase'])).mean()))
    # Calibrated: threshold at the (1-TARGET) quantile of each signal; weight so that the term's extra at the p99.5 of
    # its signal equals the median base slack on those top cycles (i.e. ~half the triggered cycles can flip).
    cal = {}
    top = 100 * (1 - TARGET)
    thr = q(c['win'], top); hi = c['win'] >= thr; k = float(np.median(slack[hi])) / max(1e-9, q(c['win'], 99.5) - thr)
    cal['win'] = dict(win_dI=thr, win_k=k)
    thr = q(c['debt'], top); hi = c['debt'] >= thr
    # debt term has no threshold in policy6 (debt_k * max(0, debt)); emulate by keeping k so it only bites near the top
    k = float(np.median(slack[hi])) / max(1e-9, q(c['debt'], 99.5))
    cal['debt'] = dict(debt_k=k)
    thr = q(c['tone'], top); hi = c['tone'] >= thr; k = float(np.median(slack[hi])) / max(1e-9, q(c['tone'], 99.5) - thr)
    cal['tone'] = dict(tone_thr=thr, tone_k=k)
    thr = q(c['bank'], 100 * TARGET); lo = c['bank'] <= thr; k = float(np.median(slack[lo])) / max(1e-9, thr - q(c['bank'], 0.5))
    cal['ready'] = dict(ready_min=thr, ready_k=k)
    if verbose:
        print('calibrated (target top %.0f%%): frac extra>0 / frac flips' % (100 * TARGET))
        for nm, p in cal.items():
            ex = contrib(P, nm, p); fl = flips(P, ex)
            print('  %-6s %s  %.4f / %.4f (with choice %.4f)' % (nm, {a: round(v, 4) for a, v in p.items()}, (ex > 0).mean(), fl.mean(), (fl & (c['nready'] > c['nbase'])).mean()))
            for b in BENCHES:
                exb = contrib(logs[b], nm, p); print('      %-9s %.4f / %.4f' % (b, (exb > 0).mean(), flips(logs[b], exb).mean()))
    return cal

def configs(cal):
    J = []
    for hold in (0, 1):
        h = dict(hold=hold)
        J.append(('TOK base h%d' % hold, h))
        for nm in ('win', 'debt', 'tone', 'ready'): J.append(('+%s h%d' % (nm, hold), dict(h, **cal[nm])))
        al = dict(h)
        for nm in ('win', 'debt', 'tone', 'ready'): al.update(cal[nm])
        J.append(('+all4 h%d' % hold, al))
        J.append(('+all4+edge15 h%d' % hold, dict(al, edge_dI=15.0)))
    return J

def job_eval(a):
    tag, extra, bench = a
    import policy6 as P6
    D, prm = setup(extra)
    E = P6.Engine(bench)
    try: r = E.run(D['banks'], 'TOK', prm)
    except RuntimeError as e:                 # policy6 ring overflow: run fell > 20000 cycles behind the trace (runaway stall)
        return dict(bench=bench, tag=tag, extra=extra, overflow=1, p999=np.nan, max=np.nan, slowdown=20000.0 / E.n, energy_over=np.nan, nreo=0, nhold=0)
    keep = {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}; keep.update(bench=bench, tag=tag, extra=extra)
    return keep

def report(res):
    tags = []
    for r in res:
        if r['tag'] not in tags: tags.append(r['tag'])
    by = {(r['tag'], r['bench']): r for r in res}
    print('%-18s %9s %9s %8s %8s %9s %9s %9s %9s %8s' % ('config', 'p999 wst', 'p999 mean', 'max wst', 'max mean', 'slow wst', 'slow mean', 'Eover mn', 'nreo', 'nhold'))
    for tg0 in tags:
        tg = tg0; R = [by[(tg, b)] for b in BENCHES if (tg, b) in by]
        g = lambda k: np.array([x[k] for x in R], float)
        nov = sum(x.get('overflow', 0) for x in R)
        if nov: tg = tg + ' [%d OVF]' % nov
        print('%-18s %9.2f %9.2f %8.2f %8.2f %8.3f%% %8.3f%% %8.4f %9d %8d' % (tg, np.nanmax(g('p999')), np.nanmean(g('p999')), np.nanmax(g('max')), np.nanmean(g('max')),
              100 * g('slowdown').max(), 100 * g('slowdown').mean(), np.nanmean(g('energy_over')), int(g('nreo').sum()), int(g('nhold').sum())))
    print('per-bench p999 / max deltas vs base (mV), nonzero only:')
    for tg in tags:
        if tg.startswith('TOK base'): continue
        base = 'TOK base ' + tg.split()[-1]
        d = ['%s OVERFLOW' % b for b in BENCHES if by[(tg, b)].get('overflow') or by[(base, b)].get('overflow')]
        d += ['%s %+.2f/%+.2f' % (b, by[(tg, b)]['p999'] - by[(base, b)]['p999'], by[(tg, b)]['max'] - by[(base, b)]['max'])
             for b in BENCHES if not (by[(tg, b)].get('overflow') or by[(base, b)].get('overflow')) and abs(by[(tg, b)]['p999'] - by[(base, b)]['p999']) + abs(by[(tg, b)]['max'] - by[(base, b)]['max']) > 1e-9]
        print('  %-18s %s' % (tg, ', '.join(d) if d else 'IDENTICAL'))

if __name__ == '__main__':
    import multiprocessing as mp
    mode = sys.argv[1]; workers = min(4, int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    if mode == 'inst':
        with mp.get_context('spawn').Pool(workers) as P: out = P.map(job_inst, BENCHES, chunksize=1)
        json.dump(dict(out), open(HERE + 'calib6_inst.json', 'w'), indent=1, default=float); print('DONE inst', flush=True)
    elif mode == 'calib':
        cal = calibrate(load_logs()); json.dump(cal, open(HERE + 'calib6_calib.json', 'w'), indent=1)
    elif mode == 'eval':
        cal = json.load(open(HERE + 'calib6_calib.json'))
        jobs = [(tg, ex, b) for b in BENCHES for (tg, ex) in configs(cal)]
        print('eval: %d jobs' % len(jobs), flush=True)
        res = []
        with mp.get_context('spawn').Pool(workers) as P:
            for r in P.imap_unordered(job_eval, jobs, chunksize=1):
                res.append(r); print('%-20s %-9s p999 %.2f max %.2f slow %.4f%s' % (r['tag'], r['bench'], r['p999'], r['max'], r['slowdown'], ' OVERFLOW' if r.get('overflow') else ''), flush=True)
                json.dump(res, open(HERE + 'calib6_eval.json', 'w'), indent=1, default=float)
        json.dump(res, open(HERE + 'calib6_eval.json', 'w'), indent=1, default=float); report(res); print('DONE eval', flush=True)
    elif mode == 'report':
        report(json.load(open(HERE + 'calib6_eval.json')))
