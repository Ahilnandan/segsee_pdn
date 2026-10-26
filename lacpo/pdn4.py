"""pdn4: real-unit PDN for the MediumBoom tile (45 nm, 1.1 V, 3 ns) with role-switchable regulators.

Off-chip path: buck power stage (Vin1 3.3 V, L1 2.2 uH, R1 40 mOhm, averaged) -> board node vo (C1 47 uF, 3 mOhm) ->
board/package R_b, L_b -> package node pk (package decap 1 uF, 100 pH, 20 mOhm) -> bumps (R_p, L_p split by region area)
-> die region nodes (or, in the series-LDO architecture, -> on-die rail nodes -> LDO -> die nodes).
Die region r: intrinsic C (ESR), explicit decap A_x[r] at 12 nF/mm^2 with series damping R_d[r], active decap Ca[r],
on-die regulator (injector or regulator role), LDO, load current. Grid between regions R_g, L_g = R_g x 10 ps.
Inductor resistance: on-die spiral 0.3 Ohm/nH; package air-core 10 mOhm + 2 mOhm/nH."""
from __future__ import print_function
import ctypes, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
LIB = ctypes.CDLL(HERE + 'libpdn4.so')
NR = 4; SHARE = np.array([0.39, 0.22, 0.33, 0.06]); LEAK_MA = np.array([22.1, 9.6, 20.8, 4.7])
DIE_MM2 = 4.0; C_INTR = 3e-9; DECAP_DENS = 12e-9; HV_DENS = 5e-9; IND_MM2_PER_NH = 0.04; CP_DENS = 3e-9   # 2.2 V reservoir: stacked thin oxide
VNOM = 1.1; T_CYC = 3e-9; SUB = 12; DT = T_CYC / SUB
EDGES = [(0, 1), (0, 2), (1, 3), (2, 3)]
_sh = np.array([0.6, 1.6, 1.9, 1.6, 1.3, 1.1, 0.9, 0.8, 0.7, 0.7, 0.7, 0.7]); SHAPE = _sh / _sh.mean()

class Glob(ctypes.Structure):
    _fields_ = [(n, ctypes.c_double) for n in ('Vin1', 'L1', 'R1', 'fc_off', 'Cb', 'slew', 'f_share', 'ff_C')] + \
               [('sense_rail', ctypes.c_int), ('rail_off', ctypes.c_double)]
class Dev(ctypes.Structure):
    _fields_ = [('has_reg', ctypes.c_int)] + [(n, ctypes.c_double) for n in ('Vin', 'L', 'Rs', 'Ipk', 'band', 'th_b')] + \
               [('dly_b', ctypes.c_int), ('dwell', ctypes.c_int), ('t_taper', ctypes.c_double), ('phases', ctypes.c_int)] + \
               [(n, ctypes.c_double) for n in ('fsw', 'fc_on', 'Cr', 'Esw', 'boost_on')] + \
               [('ldo', ctypes.c_int)] + [(n, ctypes.c_double) for n in ('V_ldo', 'R_pass', 'gm', 'tau_ldo', 'off_ldo', 'boost_ldo', 'Iq_ldo')] + \
               [('ldo_series', ctypes.c_int), ('ad_on', ctypes.c_int)] + [(n, ctypes.c_double) for n in ('G', 'tau_a', 'f_lp', 'Pq')] + \
               [('ad_br', ctypes.c_int), ('cp', ctypes.c_int)] + [(n, ctypes.c_double) for n in ('I_ch', 'R_cp', 'eta_cp', 'Vmul')]

DEFAULT = dict(Vin1=3.3, L1=2.2e-6, R1=40e-3, C1=47e-6, R_c1=3e-3, fc_off=100e3, slew=20e-3 / 1e-6, f_share=20e3,
               L_b=2e-9, R_b=3e-3, C_pk=1e-6, L_pk=100e-12, R_pk=20e-3,
               L_p=0.5e-9, R_p=50e-3, R_g=20e-3, A_x=[0.0] * NR, R_d=[0.0] * NR, A_rail=[0.0] * NR,
               c_intr=C_INTR, tau_intr=0.3e-9, tau_x=0.1e-9, series_ldo=False, rail_off=0.0, c_rail_intr=0.3e-9,
               cp_rail=False, A_cp=[0.0] * NR, V_cp0=2.15, ldo_rail=False, L_b2=2e-9, R_b2=5e-3, C_pk2=0.47e-6, L_p2=1.0e-9, R_p2=100e-3)   # separate LDO rail: fewer bumps

def dcr(L, on_die):
    return 0.3 * L / 1e-9 if on_die else 0.010 + 0.002 * L / 1e-9

def build(p, Ca=None, R_amp=1.0):
    Ca = Ca if Ca is not None else [0.0] * NR
    nodes = ['vo', 'voc', 'pk', 'pkc', 'pk2', 'pk2c'] + ['d%d' % r for r in range(NR)] + ['ci%d' % r for r in range(NR)] + \
            ['cx%d' % r for r in range(NR)] + ['bp%d' % r for r in range(NR)] + ['rn%d' % r for r in range(NR)] + ['rc%d' % r for r in range(NR)]
    ni = {n: i for i, n in enumerate(nodes)}; NN = len(nodes)
    Cn = np.full(NN, 1e-12); G = np.zeros((NN, NN)); inds = []; caps2 = []
    def res(a, b, R):
        gg = 1.0 / R; ia, ib = ni[a], ni[b]; G[ia, ia] += gg; G[ib, ib] += gg; G[ia, ib] -= gg; G[ib, ia] -= gg
    res('vo', 'voc', p['R_c1']); Cn[ni['voc']] += p['C1']
    inds.append(('vo', 'pk', p['L_b'], p['R_b']))
    inds.append(('pk', 'pkc', p['L_pk'], p['R_pk'])); Cn[ni['pkc']] += p['C_pk']
    ser = p['series_ldo']
    for r in range(NR):
        s = SHARE[r]; tgt = ('rn%d' % r) if ser else ('d%d' % r)
        inds.append(('pk', tgt, p['L_p'] / s, p['R_p'] / s))
        ci = p['c_intr'] * s; Cn[ni['ci%d' % r]] += ci; res('d%d' % r, 'ci%d' % r, p['tau_intr'] / ci)
        cx = p['A_x'][r] * DECAP_DENS
        if cx > 0: Cn[ni['cx%d' % r]] += cx; res('d%d' % r, 'cx%d' % r, max(p['tau_x'] / cx, 1e-6) + p['R_d'][r])
        else: res('d%d' % r, 'cx%d' % r, 1e3)
        if p['cp_rail']:                                   # isolated on-die reservoir node, filled by the charge pump
            cr = p['A_cp'][r] * CP_DENS + 1e-12; Cn[ni['rc%d' % r]] += cr; res('rn%d' % r, 'rc%d' % r, max(p['tau_x'] / cr, 1e-6))
        elif ser or p['ldo_rail']:
            cr = p['c_rail_intr'] * s + p['A_rail'][r] * DECAP_DENS; Cn[ni['rc%d' % r]] += cr; res('rn%d' % r, 'rc%d' % r, max(p['tau_x'] / cr, 1e-6))
            if p['ldo_rail']: inds.append(('pk2', 'rn%d' % r, p['L_p2'] / s, p['R_p2'] / s))
        else:
            res('rn%d' % r, 'pk', 1e6); res('rn%d' % r, 'rc%d' % r, 1e6)
        if Ca[r] > 0:
            caps2.append(('d%d' % r, 'bp%d' % r, Ca[r])); inds.append(('SRC%d' % r, 'bp%d' % r, 1e-12, R_amp))
        else: res('bp%d' % r, 'pk', 1e6)
    for a, b in EDGES: inds.append(('d%d' % a, 'd%d' % b, p['R_g'] * 10e-12, p['R_g']))
    if p['ldo_rail']:
        inds.append(('VRAIL', 'pk2', p['L_b2'], p['R_b2'])); inds.append(('pk2', 'pk2c', p['L_pk'], p['R_pk'])); Cn[ni['pk2c']] += p['C_pk2']
    else:
        res('pk2', 'pk', 1e6); res('pk2c', 'pk', 1e6)
    NL = len(inds); nx = NN + NL; nu = 2 + 4 * NR
    E = np.zeros((nx, nx)); A = np.zeros((nx, nx)); Bm = np.zeros((nx, nu))
    E[:NN, :NN] = np.diag(Cn); A[:NN, :NN] = -G
    for a, b, C in caps2:
        ia, ib = ni[a], ni[b]; E[ia, ia] += C; E[ib, ib] += C; E[ia, ib] -= C; E[ib, ia] -= C
    adbr = [0] * NR
    for k, (a, b, L, R) in enumerate(inds):
        j = NN + k; E[j, j] = L; A[j, j] = -R
        if a == 'VRAIL': Bm[j, 1 + 4 * NR] = 1.0
        elif isinstance(a, str) and a.startswith('SRC'): Bm[j, 1 + 2 * NR + int(a[3:])] = 1.0; adbr[int(a[3:])] = j
        else: A[j, ni[a]] += 1.0; A[ni[a], j] -= 1.0
        A[j, ni[b]] -= 1.0; A[ni[b], j] += 1.0
    Bm[ni['vo'], 0] = 1.0
    for r in range(NR):
        Bm[ni['d%d' % r], 1 + r] = -1.0; Bm[ni['d%d' % r], 1 + NR + r] = 1.0; Bm[ni['rn%d' % r], 1 + 3 * NR + r] = -1.0
    Ml = E - A * DT / 2; Pm = np.linalg.solve(Ml, E + A * DT / 2); Qm = np.linalg.solve(Ml, Bm * DT / 2)
    vidx = np.array([ni['d%d' % r] for r in range(NR)], np.int32)
    ridx = np.array([ni['rn%d' % r] if (ser or p['ldo_rail'] or p['cp_rail']) else -1 for r in range(NR)], np.int32)
    def dc(v_die, Iload_A, rail_off=0.0, v_ldo=1.2):
        """operating point with die (or rails) pinned; returns x0 and the off-chip DC current"""
        u = np.zeros(nu); u[0] = Iload_A.sum(); u[1:1 + NR] = Iload_A; u[1 + 2 * NR:1 + 3 * NR] = 0.9; u[1 + 4 * NR] = v_ldo if p['ldo_rail'] else 0.0
        if ser: u[1 + NR:1 + 2 * NR] = Iload_A; u[1 + 3 * NR:1 + 4 * NR] = Iload_A
        Ad = A.copy(); rhs = -Bm.dot(u)
        Ad[np.arange(NN), np.arange(NN)] -= 1e-12
        pins = [(ni['d%d' % r], v_die) for r in range(NR)] if ser else [(ni['vo'], None)]
        if ser:
            for idx, v in pins: Ad[idx, :] = 0; Ad[idx, idx] = 1; rhs[idx] = v
            # rails float: their level is set by vo; pin vo so that mean rail = v_die + rail_off (two passes)
        vo_guess = v_die + rail_off + 0.06
        for _ in range(3):
            A2 = Ad.copy(); r2 = rhs.copy(); iv = ni['vo']; A2[iv, :] = 0; A2[iv, iv] = 1; r2[iv] = vo_guess
            x = np.linalg.solve(A2, r2)
            sensed = np.mean([x[ni['rn%d' % r]] for r in range(NR)]) if ser else np.mean([x[vidx[r]] for r in range(NR)])
            vo_guess += (v_die + rail_off) - sensed
        return x, Iload_A.sum()
    return Pm, Qm, vidx, ridx, ni['vo'], dc, nx, adbr

def ripple_mV(spec_r, Cr):
    """analytic output ripple (mV pk-pk) of an N-phase interleaved buck in regulator role"""
    g = spec_r['reg']; N = g['phases']; Vin = g['Vin']; D = VNOM / Vin
    dI = (Vin - VNOM) * D / (g['L'] * g['fsw'])
    if N > 1:
        m = np.floor(N * D); K = N * (D - m / N) * ((m + 1) / N - D) / (D * (1 - D))
    else: K = 1.0
    return 1e3 * K * dI / (8 * N * g['fsw'] * Cr), dI

def make(spec, p, adbr):
    """spec: per-region dict or None: reg=dict(role 1|2, L, on_die, Vin, Ron, phases, fsw, fc, Ipk, th, boost, dwell, taper)
       ldo=dict(series, V, R_pass, gm, tau, off, boost, Iq)   ad=dict(Ca, G, amp)"""
    D = (Dev * NR)()
    for r in range(NR):
        d = D[r]; d.Vin = 1.8; d.L = 2e-9; d.Rs = 1.0; d.t_taper = 30e-9; d.tau_a = 1e-9; d.f_lp = 1e6; d.phases = 1
        d.tau_ldo = 1e-9; d.R_pass = 1.0
        s = spec[r] if spec else None
        if not s: continue
        if 'reg' in s:
            g = s['reg']; d.has_reg = int(g['role']); d.Vin = g['Vin']; d.L = g['L']; N = int(g.get('phases', 1)); d.phases = N
            d.Rs = g['Ron'] + dcr(g['L'], g.get('on_die', False)); d.Ipk = g.get('Ipk', 0.1); d.band = g.get('band', 0.25 * d.Ipk)
            d.th_b = g.get('th', 0.04); d.dly_b = int(g.get('dly', 4)); d.dwell = int(g.get('dwell', 2)); d.t_taper = g.get('taper', 30e-9)
            d.fsw = g.get('fsw', 100e6); d.fc_on = g.get('fc', 10e6)
            d.Cr = p['c_intr'] * SHARE[r] + p['A_x'][r] * DECAP_DENS
            d.Esw = 5e-12 * (d.Vin / 1.8) ** 2 / g['Ron']; d.boost_on = g.get('boost', 0.0)
        if 'ldo' in s:
            l = s['ldo']; d.ldo = 1; d.ldo_series = int(l.get('series', 0)); d.V_ldo = l.get('V', 1.2); d.R_pass = l['R_pass']
            d.gm = l['gm']; d.tau_ldo = l.get('tau', 1e-9); d.off_ldo = l.get('off', 0.02); d.boost_ldo = l.get('boost', 0.0); d.Iq_ldo = l.get('Iq', 0.5e-3)
            if 'cp' in l: c = l['cp']; d.cp = 1; d.I_ch = c['I_ch']; d.R_cp = c.get('R_cp', 5.0); d.eta_cp = c.get('eta', 0.9); d.Vmul = c.get('Vmul', 2.0)
        if 'ad' in s:
            a = s['ad']; d.ad_on = 1; d.G = a['G']; d.tau_a = 1e-9 / a.get('amp', 1.0); d.Pq = 1.8e-3 * a['G'] / 4 * a.get('amp', 1.0); d.ad_br = int(adbr[r])
    return D

def area(p, spec):
    a = dict(decap=float(np.sum(p['A_x'])) + float(np.sum(p.get('A_rail', [0] * NR))), active=0.0, reg=0.0, ldo=0.0,
             cp_res=float(np.sum(p.get('A_cp', [0] * NR))))
    for s in (spec or []):
        if not s: continue
        if 'ad' in s: a['active'] += s['ad']['Ca'] / HV_DENS
        if 'reg' in s:
            g = s['reg']; N = g.get('phases', 1)
            a['reg'] += 0.005 + N * 0.002 / g['Ron'] + (IND_MM2_PER_NH * g['L'] / 1e-9 * N if g.get('on_die', False) else 0.0)
        if 'ldo' in s:
            a['ldo'] += 0.005 + 0.002 / s['ldo']['R_pass']
            if 'cp' in s['ldo']: a['ldo'] += 0.003
    a['total'] = sum(a.values()); return a

def predict_triggers(Ihat, th_mV, lead, g=None):
    if g is None:
        import replay3; g = replay3.nominal_g()
    tot = Ihat.sum(1).astype(float); x = tot - np.convolve(tot, np.ones(256) / 256, 'same')
    v = np.convolve(x, g)[:len(x)]; fut = np.zeros(len(v))
    for k in range(1, lead + 1): fut[:-k] = np.minimum(fut[:-k], v[k:])
    return np.repeat((fut < -th_mV)[:, None], NR, 1).astype(np.uint8)

def simulate(I_mA, p=None, spec=None, pred=None, vref=None, role_off=None, role_on=None, warm=2000, want_v=False, vref_off=None):
    p = dict(DEFAULT, **(p or {}))
    Ca = [(spec[r]['ad']['Ca'] if (spec and spec[r] and 'ad' in spec[r]) else 0.0) for r in range(NR)]
    Pm, Qm, vidx, ridx, vo_idx, dc, nx, adbr = build(p, Ca)
    pad = lambda arr, fill: np.vstack([np.repeat(fill, warm, 0), arr]) if arr.ndim == 2 else np.concatenate([np.repeat(fill, warm), arr])
    I = np.ascontiguousarray(pad(I_mA.astype(np.float32), I_mA[:1].astype(np.float32)), np.float32); ncyc = len(I)
    Im = I[warm:warm + 3000].mean(0) * 1e-3
    v0 = VNOM if vref is None else float(vref[0])
    vl = max([spec[r]['ldo'].get('V', 1.2) for r in range(NR) if spec and spec[r] and 'ldo' in spec[r] and not spec[r]['ldo'].get('series', 0)] or [1.2])
    x0, i0 = dc(v0, Im, p['rail_off'], vl)
    if p['cp_rail']:                                         # reservoir starts charged
        for r in range(NR): x0[ridx[r]] = p['V_cp0']; x0[ridx[r] + NR] = p['V_cp0']
    D = make(spec, p, adbr)
    Gl = Glob(Vin1=p['Vin1'], L1=p['L1'], R1=p['R1'], fc_off=p['fc_off'], Cb=p['C1'] + p['C_pk'], slew=p['slew'], f_share=p['f_share'],
              ff_C=p['C1'] + p['C_pk'] + (p['c_intr'] + DECAP_DENS * sum(p['A_x'])), sense_rail=int(p['series_ldo']), rail_off=p['rail_off'])
    pr = np.zeros((ncyc, NR), np.uint8) if pred is None else np.ascontiguousarray(pad(pred.astype(np.uint8), np.zeros((1, NR), np.uint8)), np.uint8)
    vrc = None if vref is None else np.ascontiguousarray(pad(np.asarray(vref, float), np.asarray(vref[:1], float)), np.float64)
    vro = None if vref_off is None else np.ascontiguousarray(pad(np.asarray(vref_off, float), np.asarray(vref_off[:1], float)), np.float64)
    ro = None if role_off is None else np.ascontiguousarray(pad(np.asarray(role_off, np.uint8), np.asarray(role_off[:1], np.uint8)), np.uint8)
    rn = None if role_on is None else np.ascontiguousarray(pad(np.asarray(role_on, np.uint8), np.asarray(role_on[:1], np.uint8)), np.uint8)
    vmin = np.zeros((ncyc, NR), np.float32); vavg = np.zeros((ncyc, NR), np.float32); ion = np.zeros((ncyc, NR), np.float32); st = np.zeros(14)
    x0 = np.ascontiguousarray(x0, np.float64)
    dp = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)); NULLD = ctypes.POINTER(ctypes.c_double)(); NULLU = ctypes.POINTER(ctypes.c_ubyte)()
    up = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)) if a is not None else NULLU
    LIB.run.restype = ctypes.c_int
    LIB.run(ctypes.c_int(ncyc), ctypes.c_int(SUB), ctypes.c_int(nx), ctypes.c_int(NR), dp(np.ascontiguousarray(Pm)), dp(np.ascontiguousarray(Qm)),
            vidx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), ridx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), ctypes.c_int(vo_idx),
            I.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), dp(np.ascontiguousarray(LEAK_MA, np.float64)), dp(np.ascontiguousarray(SHAPE, np.float64)),
            up(pr), up(ro), up(rn), ctypes.byref(Gl), D, (dp(vrc) if vrc is not None else NULLD), (dp(vro) if vro is not None else NULLD), ctypes.c_double(VNOM), ctypes.c_double(DT),
            ctypes.c_double(i0), dp(x0), vmin.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), vavg.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ion.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), dp(st))
    vmin, vavg, ion = vmin[warm:], vavg[warm:], ion[warm:]
    ref = VNOM if vref is None else np.asarray(vref)[:, None]
    # analytic switching ripple of regulator-role converters (steady): half pk-pk added to droop
    rip = np.zeros(NR)
    for r in range(NR):
        s = spec[r] if spec else None
        if s and 'reg' in s and int(s['reg']['role']) == 2:
            Cr = p['c_intr'] * SHARE[r] + p['A_x'][r] * DECAP_DENS; rip[r] = ripple_mV(s, Cr)[0]
    dr = (ref - vmin) * 1e3 + rip[None, :] / 2; dw = dr.max(1)
    E = st[0]; f = lambda x: float(x) / E if E > 0 else 0.0
    loss = dict(off=f(st[1]), reg=f(st[2] + st[3]), inj=f(st[6] - st[7] + st[8]), ldo=f(st[9] + st[12]), ad=f(st[10] + st[11]), cp=f(st[13]))
    m = dict(p99=float(np.percentile(dw, 99)), p999=float(np.percentile(dw, 99.9)), max=float(dw.max()),
             reg_p999=[float(np.percentile(dr[:, r], 99.9)) for r in range(NR)], over_max=float(((vavg - ref) * 1e3).max()),
             loss=loss, loss_total=sum(v for k, v in loss.items() if k != 'off'), ripple_mV=rip.tolist(), E_load=E,
             mean_v=float(vavg.mean()), ion_mean=float(ion.sum(1).mean()))
    if want_v: m['vmin'] = vmin; m['vavg'] = vavg; m['ion'] = ion
    return m
