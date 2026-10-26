"""Stressmark evaluation, open loop (tile currents -> pdn6): droop per stressmark phase vs the 10 benchmarks, for
all-decap 0.35, pump-fed banks and separate-rail banks (eval6 designs). Phases are found as active regions separated by
long idle gaps (dependent-divide loops) and labelled in program order counted from the end:
  d0 d1 d2 d3 d4 d6 (first-droop sweep: chain length D) | slow (second-droop pattern) | step1..step4 (2k heavy / idle)."""
import os, sys, json, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn6, eval6
B = '/home/ahilnandan/Documents/PMIC/segseee/OutOfOrder-Optimizations/traces/lacpo_cy100/bench/'
BENCHES = ['qsort', 'rsort', 'mm', 'dhrystone', 'spmv', 'multiply', 'memcpy', 'median', 'vvadd', 'towers']
LABELS = ['d0', 'd1', 'd2', 'd3', 'd4', 'd6', 'slow', 'step1', 'step2', 'step3', 'step4']
def droop(d, I):
    sim = pdn6.Sim(d, Im_A=I[:3000].mean(0) * 1e-3)
    vmin, _, _ = sim.run(I.astype(np.float32), vref=np.full(len(I), pdn6.VNOM))
    return (pdn6.VNOM - vmin.min(1)) * 1e3
def segments(tot, gap=600):
    sm = np.convolve(tot, np.ones(64) / 64, 'same'); thr = np.percentile(sm, 20) + 0.35 * (np.percentile(sm, 95) - np.percentile(sm, 20))
    act = sm > thr; seg = []; t = 0; n = len(act)
    while t < n:
        if act[t]:
            s = t
            while t < n and (act[t] or (t + gap < n and act[t:t + gap].any())): t += 1
            seg.append((s, t))
        t += 1
    return [(a, b) for a, b in seg if b - a > 500]
def peak_mhz(x):
    x = x - x.mean(); f = np.fft.rfftfreq(len(x), 3e-9); P = np.abs(np.fft.rfft(x)) ** 2; m = f > 5e6
    return f[m][np.argmax(P[m])] / 1e6
def div_gaps(w, n, win=400):
    """Phases between idle gaps found from the issue trace: gaps = windows dominated by divide uops (fu 16)."""
    z = np.load(B + w + '_issue.npz'); t = z['t']; fu = z['fu']; ok = t < n
    dv = np.bincount(t[ok & (fu == 16)], minlength=n)[:n].astype(float); al = np.bincount(t[ok], minlength=n)[:n].astype(float)
    k = np.ones(win); gap = (np.convolve(dv, k, 'same') > 0.15 * np.maximum(np.convolve(al, k, 'same'), 1))
    seg = []; ins = False
    for i, gp in enumerate(gap):
        if not gp and not ins: s0 = i; ins = True
        elif gp and ins: seg.append((s0, i)); ins = False
    return [(a, b) for a, b in seg if b - a > 1000]
if __name__ == '__main__':
    w = sys.argv[1] if len(sys.argv) > 1 else 'stress'
    labels = LABELS if w == 'stress' else ['ring setup', 'chase1 64KB', 'chase1 2MB', 'chase3 2MB', 'chase2x 2MB']
    D, _ = eval6.designs(); designs = [('all-decap 0.35', D['decap']), ('pump banks', D['banks']), ('separate-rail banks', D['banks_sep'])]
    I = np.asarray(np.load(B + w + '_tileI.npy')).astype(float); tot = I.sum(1)
    seg = segments(tot) if w == 'stress' else div_gaps(w, len(tot)); print('found %d phases; using the last %d' % (len(seg), len(labels)))
    seg = seg[-len(labels):]; out = {'segments': {}, 'bench': {}}; LAB = labels
    dr = {nm: droop(d, I) for nm, d in designs}
    print('%-12s %8s %8s %8s | %s' % ('phase', 'cycles', 'mean mA', 'peak MHz', ' | '.join('%-22s' % nm for nm, _ in designs)))
    for lab, (a, b) in zip(LAB, seg):
        row = {nm: (float(np.percentile(dr[nm][a:b], 99.9)), float(dr[nm][a:b].max())) for nm, _ in designs}
        out['segments'][lab] = dict(a=int(a), b=int(b), mean_mA=float(tot[a:b].mean()), peak_MHz=float(peak_mhz(tot[a:b])), droop=row)
        print('%-12s %8d %8.1f %8.1f | %s' % (lab, b - a, tot[a:b].mean(), peak_mhz(tot[a:b]), ' | '.join('p999 %5.1f max %5.1f   ' % row[nm] for nm, _ in designs)))
    print('benchmarks (full traces, worst over 10):')
    for nm, d in (designs if '--bench' in sys.argv else []):
        r = []
        for wb in BENCHES:
            x = droop(d, np.asarray(np.load(B + wb + "_tileI.npy", mmap_mode="r")).astype(float))[3000:]
            r.append((wb, float(np.percentile(x, 99.9)), float(x.max())))
        out['bench'][nm] = r; wp = max(r, key=lambda z: z[1]); wm = max(r, key=lambda z: z[2])
        print('  %-20s worst p99.9 %.1f (%s)  worst max %.1f (%s)' % (nm, wp[1], wp[0], wm[2], wm[0]))
    out['whole_trace'] = {nm: (float(np.percentile(dr[nm][3000:], 99.9)), float(dr[nm][3000:].max())) for nm, _ in designs}
    print('whole %s trace:' % w, {k: tuple(round(x, 1) for x in v) for k, v in out['whole_trace'].items()})
    json.dump(out, open(os.path.dirname(os.path.abspath(__file__)) + '/stress_eval_%s.json' % w, 'w'), indent=1)
