"""Predictive issue holds on the PARALLEL (direct + FIVR) supply, verified in the switching circuit.
Scheduler: SMT-8/4 model (smt_replay) with 4-bit tokens; predictor = Pulsim-identified impulse kernel of the
supply (kernels.json) or the RLC stand-in (mismatch test). Evaluation: the scheduled current (scaled to the
per-core share) is driven through the full Pulsim circuit (case D, 16-phase FIVR at 35 or 140 MHz)."""
import json, numpy as np
from multiprocessing import Pool
import smt_exp, smt_exp2, par_exp
from smt_exp import MIXES, N
KER = json.load(open('kernels.json'))
POLS = [('A', {}), ('MPCH-fir', dict(q=0.8)), ('MPCH-fir', dict(q=0.68)), ('MPCH-rlc', dict(q=0.68)),
        ('GRO-fir', dict(lo=0.75)), ('GRO-fir', dict(lo=0.85)), ('GRO-fir', dict(lo=1.0)), ('GRO-fir', dict(lo=1.2))]
def fir_pred(I, h):
    x = I - I.mean(); return -np.convolve(x, h)[:len(I)]
def job(a):
    mix, sup, (pol, par) = a
    from smt_replay import make_thread, run, pdn_predictor
    from replay2 import droop
    from pdnchain import simulate_chain
    thr = [make_thread(w, s, N) for w, s in MIXES[mix]]
    solo = np.mean([(sum(x.sum() for x in th['sh']) / th['tr'].max()) + 27.6 for th in thr])
    h = np.array(KER[sup]); Pf = dict(kind='fir', h=h, Hp=13); est = dict(quant=True)
    RA = run(thr, 'A', pdn=Pf, est=est, record_ready=False)
    k = float(RA['I'][:RA['makespan']].mean() / solo)
    if pol == 'A': R = RA
    elif pol == 'MPCH-fir':
        a999 = float(np.percentile(fir_pred(RA['In'][:RA['makespan']], h), 99.9))
        R = run(thr, 'T', pdn=Pf, est=est, hold='MPCH', vth=par['q'] * a999, MAXD=19, record_ready=False)
    elif pol == 'MPCH-rlc':
        Lp, Cp, Rp = 507e-12 / k, 5e-9 * k, 79.6e-3 / k; P = pdn_predictor(L=Lp, C=Cp, R=Rp, Hp=13)
        a999 = float(np.percentile(droop(RA['In'], L=Lp, C=Cp, R=Rp), 99.9))
        R = run(thr, 'T', pdn=P, est=est, hold='MPCH', vth=par['q'] * a999, MAXD=19, record_ready=False)
    else:
        a999 = float(np.percentile(fir_pred(RA['In'][:RA['makespan']], h), 99.9))
        R = run(thr, 'A', pdn=Pf, est=est, gro=dict(lo=par['lo'] * a999, lat=2), MAXD=11, record_ready=False)
    fsw = 35e6 if sup == 'D35' else 140e6
    o = simulate_chain(R['In'][:R['makespan'] + 20] / k, case='D', fsw_fivr=fsw, t3=par_exp.dpar(fsw, False), warmup=3000, **par_exp.CF)
    v = 1e3 * (0.9 - o['vmin'])
    return dict(mix=mix, sup=sup, pol=pol, par=par, k=k, slow=100 * (R['makespan'] / RA['makespan'] - 1), nheld=R['nheld'],
                p999=float(np.percentile(v, 99.9)), worst=float(v.max()), top20=float(np.sort(v)[-20:].mean()),
                over=float(1e3 * (o['vmax'].max() - 0.9)))
if __name__ == '__main__':
    jobs = [(m, s, p) for m in ('smt8-mix', 'smt8-b', 'smt8-c', 'smt4-mix') for s in ('D35', 'D140') for p in POLS]
    with Pool(23) as p, open('hold_par.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
