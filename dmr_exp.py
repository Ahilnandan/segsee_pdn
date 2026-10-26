"""Full DMR (every int op re-executed within D cycles, backpressure when late) - eager vs predictor-placed checking,
single-thread BOOM closed loop in Pulsim on the final supply."""
import json, numpy as np
from multiprocessing import Pool
import wb_exp
wb_exp.CFG.clear()
wb_exp.CFG.update({
 'A (no DMR)':           dict(),
 'DMR eager':            dict(wb=dict(mode='dmr_eager', D=16)),
 'DMR placed D16':       dict(wb=dict(mode='dmr_mpc', D=16)),
 'DMR placed D64':       dict(wb=dict(mode='dmr_mpc', D=64)),
 'DMR eager +MPCH':      dict(mp=0.68, wb=dict(mode='dmr_eager', D=16)),
 'DMR placed D16 +MPCH': dict(mp=0.68, wb=dict(mode='dmr_mpc', D=16)),
})
if __name__ == '__main__':
    jobs = [(wl, nm) for wl in ('qsort', 'dhrystone', 'mm') for nm in wb_exp.CFG]
    with Pool(18) as p, open('dmr_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(wb_exp.job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
    R = [json.loads(l) for l in open('dmr_exp.jsonl')]; ref = {r['wl']: r for r in R if r['name'] == 'A (no DMR)'}
    order = list(wb_exp.CFG)
    for r in sorted(R, key=lambda r: (r['wl'], order.index(r['name']))):
        b = ref[r['wl']]; pc = lambda k: f"{100*(r[k]/b[k]-1):+5.1f}%"
        print(f"{r['wl']:9s} {r['name']:21s} | slow {r['slow']:+5.2f}% | p99.9 {r['p999']:5.1f} ({pc('p999')}) worst {r['worst']:5.1f} ({pc('worst')}) "
              f"over {r['over']:5.1f} ({pc('over')}) | check charge {r['extra_pct']:4.1f}% cov {100*r['cov']/r['nint']:.0f}%")
