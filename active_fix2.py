"""Step 2b: holds tuned per tier under a performance budget (slowdown <= 1%), then combined with the feed-forward injector."""
import json, numpy as np
from multiprocessing import Pool
import active_fix as A
def job(a):
    wl, tn = a
    import replay2
    from sweep import setup
    T = A.TIERS[tn]; inst, prods, cons, age, pin = setup(wl)
    RA = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0); mA = RA['makespan']; IA = RA['I'][:mA]
    h = A.kernel(T); a999 = float(np.percentile(-np.convolve(IA - IA.mean(), h)[:mA], 99.9))
    base = A.droop(IA, T); Ii, Pi = A.injector(IA, A.quant_hat(RA, inst)[:mA]); inj = A.droop(Ii, T)
    rows = []
    for q in (0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1):
        RH = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0, mp_h=h, mp_vth=q * a999, mp_maxd=19)
        mH = RH['makespan']; IH = RH['I'][:mH]; slow = 100 * (mH / mA - 1)
        Ib, Pb = A.injector(IH, A.quant_hat(RH, inst)[:mH])
        rows.append(dict(q=q, slow=slow, holds=A.droop(IH, T), both=A.droop(Ib, T), both_loss_mW=Pb * 1e3))
    ok = [r for r in rows if r['slow'] <= 1.0]
    bh = min(ok, key=lambda r: r['holds'][0]) if ok else None
    bb = min(ok, key=lambda r: r['both'][0]) if ok else None
    return dict(wl=wl, tier=tn, base=base, injector=inj, inj_loss_mW=Pi * 1e3, best_holds=bh, best_both=bb, sweep=rows)
if __name__ == '__main__':
    jobs = [(wl, tn) for wl in A.WL for tn in A.TIERS]
    with Pool(9) as p: R = p.map(job, jobs)
    json.dump(R, open('active_fix2.json', 'w'), indent=1)
    for r in R:
        bh, bb = r['best_holds'], r['best_both']
        print(f"{r['wl']:9s} {r['tier']:38s} | passive {r['base'][0]:5.1f} | "
              f"holds<=1% {bh['holds'][0] if bh else float('nan'):5.1f} (q {bh['q'] if bh else '-'}, {bh['slow'] if bh else 0:+.2f}%) | "
              f"injector {r['injector'][0]:5.1f} ({r['inj_loss_mW']:.1f} mW) | both {bb['both'][0] if bb else float('nan'):5.1f} / max {bb['both'][1] if bb else float('nan'):5.1f} "
              f"(q {bb['q'] if bb else '-'}, {bb['slow'] if bb else 0:+.2f}%)")
