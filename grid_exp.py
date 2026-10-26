"""Centralised vs distributed micro-regulators (global vs local control) on a 4x4 on-die grid, full traces."""
import json, numpy as np
from multiprocessing import Pool
CFG = [('central/global', 'central', 'global'), ('distributed/global', 'distributed', 'global'), ('distributed/local', 'distributed', 'local')]
def job(a):
    wl, Rg, (name, place, ctl) = a
    import replay2
    from sweep import setup
    from gridpdn import simulate_grid
    inst, prods, cons, age, pin = setup(wl)
    R = replay2.run(inst, prods, cons, age, pin, 'A', 1e9, 0, 0, tiles=True); It = R['It'][:R['makespan']]
    o = simulate_grid(It, place=place, ctl=ctl, Rg=Rg, Lg=Rg * 10e-12)   # on-die mesh is resistive: L/R ~ 10 ps
    v = 1e3 * (0.9 - o['vmin']); g = 1e3 * (0.9 - o['vmin'].mean(1))
    per = np.percentile(v, 99.9, axis=0)
    return dict(wl=wl, Rg=Rg, name=name, worst_tile_p999=float(np.percentile(v.max(1), 99.9)), mean_tile_p999=float(np.percentile(g, 99.9)),
                worst=float(v.max()), per_tile_p999=per.round(1).tolist(), hot=int(per.argmax()),
                spread_p999=float(np.percentile(v.max(1) - v.mean(1), 99.9)),
                over=float(1e3 * (o['vmax'].max() - 0.9)))
if __name__ == '__main__':
    jobs = [(wl, Rg, c) for wl in ('qsort', 'dhrystone', 'mm') for Rg in (0.05, 0.2, 0.5) for c in CFG]
    with Pool(24) as p, open('grid_exp.jsonl', 'w') as f:
        for r in p.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush()
    R = [json.loads(l) for l in open('grid_exp.jsonl')]
    for r in sorted(R, key=lambda r: (r['wl'], r['Rg'], [c[0] for c in CFG].index(r['name']))):
        print(f"{r['wl']:9s} Rg {r['Rg']:.2f} {r['name']:19s} | worst-tile p99.9 {r['worst_tile_p999']:6.1f} mean-tile {r['mean_tile_p999']:6.1f} "
              f"local excess p99.9 {r['spread_p999']:5.1f} | max {r['worst']:6.1f} over {r['over']:6.1f} | hottest tile {r['hot']}")
