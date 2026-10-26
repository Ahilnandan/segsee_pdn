"""How much does work-conserving reordering buy when there is Fu-et-al.-level issue choice?"""
import json, sys, numpy as np
from multiprocessing import Pool
N = 60000
MIXES = {
    'smt4-mix':  [('qsort', 50000), ('dhrystone', 50000), ('mm', 50000), ('mm', 200000)],
    'smt4-int':  [('qsort', 50000), ('qsort', 250000), ('dhrystone', 20000), ('dhrystone', 150000)],
    'smt4-mm':   [('mm', 30000), ('mm', 110000), ('mm', 190000), ('qsort', 150000)],
    'smt8-mix':  [('qsort', 50000), ('dhrystone', 50000), ('mm', 50000), ('mm', 200000),
                  ('qsort', 250000), ('dhrystone', 150000), ('mm', 110000), ('qsort', 150000)],
}
def metrics(R, k):
    from replay2 import droop
    I = R['I']; v = droop(I, L=507e-12 / k, C=5e-9 * k, R=79.6e-3 / k)
    x = I[200:R['makespan']] - I[200:R['makespan']].mean()
    F = np.abs(np.fft.rfft(x)) ** 2; f = np.fft.rfftfreq(len(x), 1e-9)
    band = F[(f > 50e6) & (f < 200e6)].sum() / F[1:].sum()
    dI = np.abs(np.diff(I[:R['makespan']]))
    return dict(droop=float(v.max()), p999=float(np.percentile(v, 99.9)), p99=float(np.percentile(v, 99)),
                band=float(band), dIp999=float(np.percentile(dI, 99.9)), Imean=float(I[:R['makespan']].mean()))
def job(a):
    mix, pol, maxd = a
    from smt_replay import make_thread, run
    thr = [make_thread(w, s, N) for w, s in MIXES[mix]]
    # per-core-share scaling: the SMT core draws k x a single BOOM core, so its PDN (C*k, L/k, R/k),
    # bank threshold and recharge rate scale by k -> droop numbers comparable to the single-thread work
    solo = np.mean([(sum(x.sum() for x in th['sh']) / th['tr'].max()) + 27.6 for th in thr])
    R0 = run(thr, 'A', MAXD=maxd, record_ready=False)
    k = float(R0['I'][:R0['makespan']].mean() / solo)
    TH, IR = 100.0 * k, 25.0 * k
    R = run(thr, pol, MAXD=maxd, TH=TH)
    out = dict(mix=mix, pol=pol, maxd=maxd, makespan=R['makespan'], fin=R['fin'], nchg=R['nchg'],
               ready_mean=float(R['nready'].mean()), ready_lt9=float((R['nready'] < 9).mean()), ready_max=int(R['nready'].max()))
    # bank demand: unlimited bank -> total charge it must supply; then the smallest bank with no shortfall
    RB = run(thr, pol, MAXD=maxd, TH=TH, Qb=1e12, Ir=1e12, record_ready=False); out['bank_charge_nC'] = RB['drawn'] / 1e3
    lo, hi = 0.0, 4000.0 * k
    for _ in range(11):
        mid = (lo + hi) / 2; Rq = run(thr, pol, MAXD=maxd, TH=TH, Qb=mid, Ir=IR, record_ready=False)
        if Rq['unc'] <= 0.001 * max(RB['drawn'], 1.0): hi = mid
        else: lo = mid
    out['min_bank_pC'] = hi; out['min_bank_pC_per_k'] = hi / k
    out.update(metrics(R, k)); out['k'] = k
    return out
if __name__ == '__main__':
    jobs = [(m, p, md) for m in MIXES for p in ('A', 'T', 'FLAT') for md in ((8,) if p == 'A' else (8, 32))]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 20) as p, open('smt_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush(); print(json.dumps(r), flush=True)
