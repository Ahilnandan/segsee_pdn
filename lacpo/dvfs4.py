"""DVFS reaction time per architecture (pdn4): reference steps 1.0 -> 1.1 V at S1 and back at S2 while qsort runs.
Settle = first time after the step from which the 300 ns moving average of the mean die voltage stays within +-1 % of the
new reference for 1 us. Also reports droop (p99.9 over the run, vs the reference), loss, and transition overshoot."""
import os, sys, json, numpy as np
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdn4, ga3
T = pdn4.T_CYC; S1, S2, N = 3000, 13000, 23000
def settle(v, ref, k0, k1):
    v = np.convolve(v, np.ones(100) / 100, 'same'); tgt = ref[k0 + 5]; ok = np.abs(v[k0:k1] - tgt) < 0.01 * tgt; need = int(1e-6 / T)
    for k in range(len(ok) - need):
        if ok[k:k + need].all(): return k * T * 1e9
    return float('nan')
def case(name, p, spec=None, role_off=None, role_on=None, vref_off=None, Ax=0.4):
    z = np.load(ga3.B + 'qsort_rep_A.npz'); I = z['I'][60000:60000 + N].astype(np.float32)
    vref = np.full(N, 1.0); vref[S1:S2] = 1.1
    pp = dict(dict(A_x=[Ax * s for s in pdn4.SHARE]), **p)
    m = pdn4.simulate(I, pp, spec, vref=vref, role_off=role_off, role_on=role_on, vref_off=vref_off, want_v=True, warm=4000)
    v = m['vavg'].mean(1)
    return dict(name=name, up_ns=settle(v, vref, S1, S2), dn_ns=settle(v, vref, S2, N), p999=m['p999'], loss=m['loss_total'] * 100,
                over=float((v[S1:S2].max() - 1.1) * 1e3), under=float((1.0 - v[S2:].min()) * 1e3), area=pdn4.area(dict(pdn4.DEFAULT, **pp), spec)['total'])
reg = lambda role, L=100e-9, fsw=100e6, fc=20e6, N=4, Ron=1.0, Ax=0.4: [dict(reg=dict(role=role, L=L, on_die=False, Vin=1.8, Ron=Ron, phases=N, fsw=fsw, fc=fc)) for s in pdn4.SHARE]
ldo_s = lambda: [dict(ldo=dict(series=1, R_pass=0.2, gm=10.0, tau=1e-9, off=0.0)) for s in pdn4.SHARE]
ldo_p = lambda: [dict(ldo=dict(series=0, V=1.2, R_pass=0.3, gm=5.0, tau=1e-9, off=0.02)) for s in pdn4.SHARE]
ones = np.ones(N, np.uint8)
# B role swap: off-chip = current helper only during the transition windows (+-3 us), on-die = regulator throughout
sw = np.zeros(N, np.uint8); sw[S1 - 200:S1 + 1000] = 1; sw[S2 - 200:S2 + 1000] = 1
role_on_B = np.repeat(np.where(sw > 0, 2, 1)[:, None], 4, 1).astype(np.uint8)       # injector normally, regulator in transitions
vpre = np.full(N, 1.0); vpre[S1 - 2000:S2] = 1.1                                    # series LDO: rail raised 6 us ahead
CASES = [
    ('A  off-chip sets V (20 mV/us)', dict(), None),
    ('A  off-chip sets V (50 mV/us)', dict(slew=50e-3 / 1e-6), None),
    ('B  role swap: off-chip -> current helper in transitions, on-die regulates', dict(), reg(2), dict(role_off=sw, role_on=role_on_B)),
    ('B2 permanent: off-chip current helper, on-die regulator', dict(), reg(2), dict(role_off=ones)),
    ('C  on-die FIVR only (series), 100 nH 4ph 100 MHz fc 20 MHz', dict(), None, 'C'),
    ('D  series LDO, rail +80 mV follows off-chip', dict(series_ldo=True, rail_off=0.08, A_rail=[0.1 * s for s in pdn4.SHARE]), ldo_s()),
    ('D2 series LDO, rail +80 mV pre-raised 6 us ahead', dict(series_ldo=True, rail_off=0.08, A_rail=[0.1 * s for s in pdn4.SHARE]), ldo_s(), dict(vref_off=vpre)),
    ('E  A + parallel LDO clamp (1.2 V rail)', dict(ldo_rail=True), ldo_p()),
]
def run(c):
    name, p, spec = c[:3]; kw = c[3] if len(c) > 3 else {}
    if kw == 'C':          # FIVR only: bumps feed rail nodes at 1.8 V via the series-LDO plumbing, regulator delivers all current
        return dict(name=name, note='see FIVR-only results in dvfs3 (pdn3); series-FIVR plumbing not in pdn4')
    return case(name, p, spec, **kw)
if __name__ == '__main__':
    import multiprocessing as mp
    with mp.get_context('spawn').Pool(8) as P: res = P.map(run, CASES)
    for r in res:
        if 'note' in r: print('%-78s %s' % (r['name'], r['note'])); continue
        print('%-78s up %7.0f ns  down %7.0f ns  p99.9 %6.1f mV  loss %5.1f%%  overshoot %5.1f  undershoot %5.1f mV  area %.3f' % (
            r['name'], r['up_ns'], r['dn_ns'], r['p999'], r['loss'], r['over'], r['under'], r['area']))
    json.dump(res, open('dvfs4.json', 'w'), indent=1, default=float)
