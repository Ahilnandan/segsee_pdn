"""Active bank inside the parallel supply (FIVR 35 MHz, 4x L) + predictive holds, verified in Pulsim."""
import json, numpy as np
from multiprocessing import Pool
import smt_exp, smt_exp2, par_exp
from smt_exp import MIXES, N
KER = json.load(open('kernels.json'))
KW = dict(case='D', fsw_fivr=35e6, t3=par_exp.dpar(35e6, False, 448e-9), **dict(par_exp.CF, L_fivr=448e-9))
SUP = {'D35L4': dict(cdie=5e-9), 'D35L4+bank': dict(cdie=4.5e-9, bank=dict(cb=0.5, ir=25.0))}
POLS = [('A', None, None), ('MPCH', 'D35L4', 0.68), ('MPCH', 'D35L4', 0.8), ('MPCH', 'D35L4bank', 0.68)]
def fir_pred(I, h): x = I - I.mean(); return -np.convolve(x, h)[:len(I)]
def job(a):
    mix, sup, (pol, ker, q) = a
    from smt_replay import make_thread, run
    from pdnchain import simulate_chain
    thr = [make_thread(w, s, N) for w, s in MIXES[mix]]
    solo = np.mean([(sum(x.sum() for x in th['sh']) / th['tr'].max()) + 27.6 for th in thr])
    est = dict(quant=True)
    h0 = np.array(KER['D35L4'])
    RA = run(thr, 'A', pdn=dict(kind='fir', h=h0, Hp=13), est=est, record_ready=False)
    k = float(RA['I'][:RA['makespan']].mean() / solo)
    if pol == 'A': R = RA
    else:
        h = np.array(KER[ker]); a999 = float(np.percentile(fir_pred(RA['In'][:RA['makespan']], h), 99.9))
        R = run(thr, 'T', pdn=dict(kind='fir', h=h, Hp=13), est=est, hold='MPCH', vth=q * a999, MAXD=19, record_ready=False)
    o = simulate_chain(R['In'][:R['makespan'] + 20] / k, warmup=3000, **KW, **SUP[sup])
    v = 1e3 * (0.9 - o['vmin'])
    return dict(mix=mix, sup=sup, pol=pol, ker=ker, q=q, slow=100 * (R['makespan'] / RA['makespan'] - 1), nheld=R['nheld'],
                p999=float(np.percentile(v, 99.9)), worst=float(v.max()), over=float(1e3 * (o['vmax'].max() - 0.9)))
if __name__ == '__main__':
    jobs = [(m, s, p) for m in ('smt8-mix', 'smt8-b', 'smt8-c', 'smt4-mix') for s in SUP for p in POLS
            if not (p[1] == 'D35L4bank' and s == 'D35L4')]
    with Pool(23) as p, open('hold_par2.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
    R = [json.loads(l) for l in open('hold_par2.jsonl')]
    ref = {r['mix']: r for r in R if r['pol'] == 'A' and r['sup'] == 'D35L4'}
    for r in sorted(R, key=lambda r: (r['mix'], r['sup'], str(r['ker']), str(r['q']))):
        b = ref[r['mix']]
        print(f"{r['mix']:9s} {r['sup']:11s} {r['pol']:4s} {str(r['ker']):10s} {str(r['q']):5s} | slow {r['slow']:+5.2f}% | p99.9 {r['p999']:5.1f} ({100*(r['p999']/b['p999']-1):+5.1f}%) "
              f"worst {r['worst']:5.1f} ({100*(r['worst']/b['worst']-1):+5.1f}%) over {r['over']:5.1f}")
