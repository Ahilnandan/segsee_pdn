import json, numpy as np, sys
from multiprocessing import Pool
import st_exp, wb_exp
_O = {}
def o999(wl):
    if wl in _O: return _O[wl]
    import replay2
    from sweep import setup
    inst, prods, cons, age, pin = setup(wl)
    R = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0); m = R['makespan']; I = R['I'][:m]
    _O[wl] = float(np.percentile(np.convolve(I - I.mean(), st_exp.KER['D35L4'])[:m], 99.9)); return _O[wl]
def job(a):
    wl, qo, mp = a
    name = f'redund vo{qo} ' + ('+MPCH' if mp else 'alone')
    wb_exp.CFG[name] = dict(wb=dict(mode='redund', vo=qo * o999(wl)), **(dict(mp=0.68) if mp else {}))
    return wb_exp.job((wl, name))
if __name__ == '__main__':
    jobs = [(wl, qo, mp) for wl in ('qsort', 'dhrystone', 'mm') for qo, mp in ((0.6, True), (0.8, True), (0.6, False))]
    with Pool(7) as p, open('wb_redund.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
