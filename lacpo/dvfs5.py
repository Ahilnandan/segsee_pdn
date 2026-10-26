"""Phase 4.2: DVFS step study with the role-split on-die buck (libpdn5b.so, current-limited).
Rail steps 1.10 -> 1.00 V at cycle S1 and back at S2, under a real qsort load window. Compares:
  off    : off-chip loop only (slew-limited reference)
  reg2   : on-die role-2 regulator on 8 core tiles, off-chip loop still running
  inj1   : on-die role-1 injector on the same tiles, off-chip in current-share mode (role_off=1)
Metrics: first entry into +-25 mV of the new target (ns), steady error over the last 1000 cycles of each segment, transition min/max."""
import os, sys, json
os.environ['PDN5_LIB'] = 'libpdn5b.so'
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn5
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
T = pdn5.T; N = 12000; S1, S2 = 2000, 7000; V_HI, V_LO = 1.10, 1.00
CORE = [r * 8 + c for c in (3, 4, 5) for r in range(0, 8, 3)][:8]
def buck_dict():
    return {t: dict(role=2, Vin=1.3, L=2e-9, Ron=0.05, Ipk=0.05, th=0.04, dwell=2, fc=30e6, phases=2) for t in CORE}
def design(mode):
    d = dict(A_x=0.3 * pdn5.CELL_SHARE)
    if mode in ('reg2', 'inj1'): d['buck'] = buck_dict() if mode == 'reg2' else {t: dict(role=1, Vin=1.3, L=2e-9, Ron=0.05, Ipk=0.05, th=0.04, dwell=2, fc=30e6, phases=2) for t in CORE}
    return d
def run(mode):
    z = np.load(B + 'qsort_attr5.npz'); BG = z['BG'].astype(np.float64)
    load = BG[5000:5000 + N].copy() if len(BG) > 5000 + N else BG[:N].copy()
    im0 = BG[:3000].mean(0) * 1e-3
    sim = pdn5.Sim(design(mode), Im_A=im0)
    vref = np.where(np.arange(N) < S1, V_HI, np.where(np.arange(N) < S2, V_LO, V_HI))
    role_on = np.zeros((N, T), np.uint8)
    if mode in ('reg2', 'inj1'):
        role_on[:] = 0
        for t in CORE: role_on[:, t] = 2 if mode == 'reg2' else 1
    role_off = np.full(N, 1 if mode == 'inj1' else 0, np.uint8)
    vmin, vavg, _ = sim.run(load.astype(np.float32), vref=vref.astype(np.float64), role_off=role_off, role_on=role_on)
    vdie = vavg.mean(1); return vdie, vref
def entry(v, tgt, t0, t1, band=0.025):
    ok = np.where(np.abs(v[t0:t1] - tgt) <= band)[0]
    return int(ok[0]) if len(ok) else None
if __name__ == '__main__':
    out = {}
    for mode in ('off', 'reg2', 'inj1'):
        v, vref = run(mode)
        ts1 = entry(v, V_LO, S1, S2); ts2 = entry(v, V_HI, S2, N)
        ss1 = float(v[S2 - 1000:S2].mean() - V_LO) * 1e3; ss2 = float(v[N - 1000:].mean() - V_HI) * 1e3
        seg1 = v[S1:S2]; seg2 = v[S2:]
        ns = lambda c: None if c is None else round(c * pdn5.T_CYC * 1e9, 1)
        out[mode] = dict(entry25_down_ns=ns(ts1), entry25_up_ns=ns(ts2), steady_err_down_mV=round(ss1, 2), steady_err_up_mV=round(ss2, 2),
                         min_down_mV=float((seg1.min() - V_LO) * 1e3), max_down_mV=float((seg1.max() - V_LO) * 1e3),
                         max_up_mV=float((seg2.max() - V_HI) * 1e3), min_up_mV=float((seg2.min() - V_HI) * 1e3))
        print(mode, out[mode], flush=True)
    json.dump(out, open(HERE + 'dvfs5.json', 'w'), indent=1)
