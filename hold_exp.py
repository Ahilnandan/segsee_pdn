"""Predictive holds (MPCH) vs reactive bank-gated holds (TDH) vs issue-count ramp limiting (RAMP):
droop vs slowdown trade-off, SMT model, bank 0.5x / 1.0x of A's minimum."""
import json, numpy as np
from multiprocessing import Pool
import smt_exp, smt_exp2
from smt_exp import MIXES, N
from soc_exp import MINB
def job(a):
    mix, frac, hold, par = a
    from smt_replay import make_thread, run, pdn_predictor
    from replay2 import droop
    thr = [make_thread(w, s, N) for w, s in MIXES[mix]]
    mb, k = MINB[mix]; Qb = frac * mb; TH, IR = 100.0 * k, 25.0 * k
    Lp, Cp, Rp = 507e-12 / k, 5e-9 * k, 79.6e-3 / k; P = pdn_predictor(L=Lp, C=Cp, R=Rp)
    kw = dict(TH=TH, Qb=Qb, Ir=IR, record_ready=False, pdn=P)
    RA = run(thr, 'A', **kw); vA = droop(RA['In'], L=Lp, C=Cp, R=Rp); a999 = float(np.percentile(vA, 99.9))
    if hold is None: R = RA
    elif hold == 'MPCH': R = run(thr, 'A', hold='MPCH', vth=par * a999, **kw)
    elif hold == 'TDH': R = run(thr, 'A', hold='TDH', hth=par, **kw)
    else: R = run(thr, 'A', hold='RAMP', ramp=int(par), **kw)
    v = droop(R['In'], L=Lp, C=Cp, R=Rp)
    return dict(mix=mix, frac=frac, hold=hold, par=par, slow=100 * (R['makespan'] / RA['makespan'] - 1), nheld=R['nheld'],
                p99=float(np.percentile(v, 99)), p999=float(np.percentile(v, 99.9)), worst=float(v.max()),
                top20=float(np.sort(v)[-20:].mean()), unc_nC=R['unc'] / 1e3)
if __name__ == '__main__':
    P = [(None, 0)] + [('MPCH', q) for q in (0.95, 0.9, 0.85, 0.8, 0.7)] + [('TDH', h) for h in (1.0, 0.8, 0.6, 0.4, 0.25)] + \
        [('RAMP', d) for d in (4, 3, 2, 1)]
    jobs = [(m, f, h, p) for m in ('smt8-mix', 'smt8-b', 'smt8-c', 'smt4-mix') for f in (0.5, 1.0) for h, p in P]
    with Pool(22) as p, open('hold_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
