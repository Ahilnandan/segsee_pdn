"""Model-predictive (droop-objective) selection vs oldest-first / SOC, undersized and full bank, SMT model."""
import json, numpy as np
from multiprocessing import Pool
import smt_exp, smt_exp2
from smt_exp import MIXES, N
from soc_exp import MINB
def job(a):
    mix, frac, pol, res, q = a
    from smt_replay import make_thread, run, pdn_predictor
    from replay2 import droop
    thr = [make_thread(w, s, N) for w, s in MIXES[mix]]
    mb, k = MINB[mix]; Qb = frac * mb; TH, IR = 100.0 * k, 25.0 * k
    Lp, Cp, Rp = 507e-12 / k, 5e-9 * k, 79.6e-3 / k; P = pdn_predictor(L=Lp, C=Cp, R=Rp)
    RA = run(thr, 'A', TH=TH, Qb=Qb, Ir=IR, record_ready=False)
    vA = droop(RA['In'], L=Lp, C=Cp, R=Rp); vth = float(np.percentile(vA, q)) if q else None
    R = RA if pol == 'A' else run(thr, pol, TH=TH, Qb=Qb, Ir=IR, res=res, pdn=P, vth=vth, record_ready=False)
    v = droop(R['In'], L=Lp, C=Cp, R=Rp)
    return dict(mix=mix, frac=frac, pol=pol, res=res, q=q, vth=vth, slow=100 * (R['makespan'] / RA['makespan'] - 1),
                unc_nC=R['unc'] / 1e3, nshort=R['nshort'], p99=float(np.percentile(v, 99)), p999=float(np.percentile(v, 99.9)),
                worst=float(v.max()), top20=float(np.sort(v)[-20:].mean()), nchg=R['nchg'])
if __name__ == '__main__':
    P = [('A', 0, 0), ('SOC', 0.3, 0), ('MPCmin', 0.3, 0), ('MPCmin', 0.0, 0)] + \
        [('MPC', r, q) for r in (0.0, 0.3) for q in (99, 99.9)]
    jobs = [(m, f, p, r, q) for m in ('smt8-mix', 'smt8-b', 'smt8-c', 'smt4-mix') for f in (0.5, 1.0) for p, r, q in P]
    with Pool(22) as p, open('mpc_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
