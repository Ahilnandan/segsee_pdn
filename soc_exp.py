"""State-of-charge-aware scheduling vs oldest-first / min-need (T), undersized bank, SMT model."""
import json, sys, numpy as np
from multiprocessing import Pool
import smt_exp, smt_exp2
from smt_exp import MIXES, N
MINB = {}
for fn in ('smt_exp.jsonl', 'smt_exp2.jsonl'):
    for r in map(json.loads, open(fn)):
        if r['pol'] == 'A': MINB[r['mix']] = (r['min_bank_pC'], r['k'])
def job(a):
    mix, frac, pol, res = a
    from smt_replay import make_thread, run
    from replay2 import droop
    thr = [make_thread(w, s, N) for w, s in MIXES[mix]]
    mb, k = MINB[mix]; Qb = frac * mb; TH, IR = 100.0 * k, 25.0 * k
    R = run(thr, pol, TH=TH, Qb=Qb, Ir=IR, res=res, record_ready=False)
    RA = run(thr, 'A', TH=TH, record_ready=False)
    v = droop(R['In'], L=507e-12 / k, C=5e-9 * k, R=79.6e-3 / k)
    I = R['In'][:R['makespan']]; dI = np.diff(I)
    return dict(mix=mix, frac=frac, pol=pol, res=res, Qb=Qb, k=k, slow=100 * (R['makespan'] / RA['makespan'] - 1),
                unc_nC=R['unc'] / 1e3, drawn_nC=R['drawn'] / 1e3, nshort=R['nshort'],
                p999=float(np.percentile(v, 99.9)), worst=float(v.max()), top20=float(np.sort(v)[-20:].mean()),
                dIup999=float(np.percentile(dI, 99.9)), dIdn999=float(-np.percentile(dI, 0.1)),
                soc_min_frac=float(R['soc'][200:R['makespan']].min() / Qb), soc_p1=float(np.percentile(R['soc'][200:R['makespan']], 1) / Qb))
if __name__ == '__main__':
    P = [('A', 0), ('T', 0)] + [(p, r) for p in ('RES', 'SOC') for r in (0.0, 0.3, 0.6)]
    jobs = [(m, f, p, r) for m in ('smt8-mix', 'smt8-b', 'smt8-c', 'smt4-mix') for f in (0.5, 1.0) for p, r in P]
    with Pool(22) as p, open('soc_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
