"""Real-unit PDN for the MediumBoom tile (45 nm, 1.1 V, 3 ns): network builder (descriptor MNA -> trapezoidal),
C++ kernel wrapper (libpdn3.so), droop / energy / area metrics. All values SI unless the name says otherwise.

Topology (one big off-chip path + optional on-die parallel path per region):
  off-chip VR (integral loop on mean die V, crossover fc_vr) -- board/package R_b,L_b -- package node pk
  pk: package decap (ESL, ESR) + bulk; pk -> each die region r through bumps (R_p/s_r, L_p/s_r)
  die region r (2x2 floorplan R0 R1 / R2 R3, grid R_g, L_g = R_g x 10 ps): intrinsic C (ESR), explicit decap
  area A_x[r] (12 nF/mm^2, series damping R_d[r]), load current, on-die devices (small buck from 1.8 V, active decap)."""
from __future__ import print_function
import ctypes, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
LIB = ctypes.CDLL(HERE + os.environ.get('PDN3_LIB', 'libpdn3.so'))
NEWSIG = os.environ.get('PDN3_LIB', 'libpdn3.so') != 'libpdn3.so'
NR = 4; SHARE = np.array([0.39, 0.22, 0.33, 0.06]); LEAK_MA = np.array([22.1, 9.6, 20.8, 4.7])
DIE_MM2 = 4.0; C_INTR = 3e-9; DECAP_DENS = 12e-9; HV_DENS = 5e-9; IND_MM2_PER_NH = 0.04
VNOM = 1.1; T_CYC = 3e-9; SUB = 12; DT = T_CYC / SUB
EDGES = [(0, 1), (0, 2), (1, 3), (2, 3)]
_sh = np.array([0.6, 1.6, 1.9, 1.6, 1.3, 1.1, 0.9, 0.8, 0.7, 0.7, 0.7, 0.7]); SHAPE = _sh / _sh.mean()   # intra-cycle profile (front-loaded)

class Dev(ctypes.Structure):
    _fields_ = [('buck_mode', ctypes.c_int), ('Vin', ctypes.c_double), ('L', ctypes.c_double), ('Rs', ctypes.c_double),
                ('Ipk', ctypes.c_double), ('band', ctypes.c_double), ('th_b', ctypes.c_double), ('dly_b', ctypes.c_int),
                ('use_pred_b', ctypes.c_int), ('dwell', ctypes.c_int), ('t_taper', ctypes.c_double),
                ('fsw', ctypes.c_double), ('fc_b', ctypes.c_double), ('Esw', ctypes.c_double),
                ('ad_on', ctypes.c_int), ('G', ctypes.c_double), ('tau_a', ctypes.c_double), ('f_lp', ctypes.c_double),
                ('Pq', ctypes.c_double), ('ad_br', ctypes.c_int), ('Cr', ctypes.c_double), ('no_int', ctypes.c_int)]

DEFAULT = dict(L_b=2e-9, R_b=3e-3, C_pk=1e-6, L_pk=100e-12, R_pk=20e-3, C_blk=47e-6, R_blk=5e-3,
               L_p=0.5e-9, R_p=50e-3, R_g=20e-3, A_x=[0.0] * NR, R_d=[0.0] * NR, fc_vr=100e3,
               die_mm2=DIE_MM2, c_intr=C_INTR, tau_intr=0.3e-9, tau_x=0.1e-9, slew_vr=20e-3 / 1e-6, fivr_only=False)

def build(p, Ca=None, R_amp=1.0):
    """returns P, Q, vidx (die nodes + VR branch current index), dc solver, nx, ad branch indices.
    Ca[r] > 0 adds an active-decap capacitor between die node d_r and bottom-plate node bp_r, which is driven by the
    amplifier voltage input u[1+2R+r] through R_amp (an inductor branch with 1 pH so the drive is a state)."""
    Ca = Ca if Ca is not None else [0.0] * NR
    nodes = ['pk', 'pkc', 'blk'] + ['d%d' % r for r in range(NR)] + ['ci%d' % r for r in range(NR)] + ['cx%d' % r for r in range(NR)] + ['bp%d' % r for r in range(NR)]
    ni = {n: i for i, n in enumerate(nodes)}; NN = len(nodes)
    Cn = np.full(NN, 1e-12)                                   # parasitic floor on every node
    G = np.zeros((NN, NN)); inds = []                          # (from, to, L, R) ; from None = VR source
    def res(a, b, R):
        g = 1.0 / R; ia, ib = ni[a], ni[b]; G[ia, ia] += g; G[ib, ib] += g; G[ia, ib] -= g; G[ib, ia] -= g
    inds.append((None, 'pk', p['L_b'], p['R_b']))
    inds.append(('pk', 'pkc', p['L_pk'], p['R_pk'])); Cn[ni['pkc']] += p['C_pk']
    res('pk', 'blk', p['R_blk']); Cn[ni['blk']] += p['C_blk']
    sh = SHARE
    for r in range(NR):
        inds.append(('pk', 'd%d' % r, p['L_p'] / sh[r], (1e6 if p['fivr_only'] else p['R_p'] / sh[r])))
        ci = p['c_intr'] * sh[r]; Cn[ni['ci%d' % r]] += ci; res('d%d' % r, 'ci%d' % r, p['tau_intr'] / ci)
        cx = p['A_x'][r] * DECAP_DENS
        if cx > 0:
            Cn[ni['cx%d' % r]] += cx; res('d%d' % r, 'cx%d' % r, max(p['tau_x'] / cx, 1e-6) + p['R_d'][r])
        else:
            res('d%d' % r, 'cx%d' % r, 1e3)
    for a, b in EDGES: inds.append(('d%d' % a, 'd%d' % b, p['R_g'] * 10e-12, p['R_g']))
    caps2 = []                                                  # capacitors between two non-ground nodes
    adbr = [-1] * NR
    for r in range(NR):
        if Ca[r] > 0:
            caps2.append(('d%d' % r, 'bp%d' % r, Ca[r])); adbr[r] = len(inds); inds.append(('SRC%d' % r, 'bp%d' % r, 1e-12, R_amp))
        else:
            res('bp%d' % r, 'pk', 1e6)
    NL = len(inds); nx = NN + NL; nu = 1 + 3 * NR
    E = np.zeros((nx, nx)); A = np.zeros((nx, nx)); Bm = np.zeros((nx, nu))
    E[:NN, :NN] = np.diag(Cn); A[:NN, :NN] = -G
    for a, b, C in caps2:
        ia, ib = ni[a], ni[b]; E[ia, ia] += C; E[ib, ib] += C; E[ia, ib] -= C; E[ib, ia] -= C
    for k, (a, b, L, R) in enumerate(inds):
        j = NN + k; E[j, j] = L; A[j, j] = -R
        if a is None: Bm[j, 0] = 1.0
        elif isinstance(a, str) and a.startswith('SRC'): Bm[j, 1 + 2 * NR + int(a[3:])] = 1.0
        else: A[j, ni[a]] += 1.0; A[ni[a], j] -= 1.0
        A[j, ni[b]] -= 1.0; A[ni[b], j] += 1.0
    for r in range(NR):
        Bm[ni['d%d' % r], 1 + r] = -1.0; Bm[ni['d%d' % r], 1 + NR + r] = 1.0
    Ml = E - A * DT / 2; Pm = np.linalg.solve(Ml, E + A * DT / 2); Qm = np.linalg.solve(Ml, Bm * DT / 2)
    vidx = np.array([ni['d%d' % r] for r in range(NR)] + [NN], np.int32)
    def dc(vr, Iload_A):
        u = np.zeros(nu); u[0] = vr; u[1:1 + NR] = Iload_A; u[1 + 2 * NR:] = 0.9
        Ad = A.copy(); Ad[np.arange(NN), np.arange(NN)] -= 1e-12      # tiny leak so floating cap pairs have a DC solution
        return np.linalg.solve(Ad, -Bm.dot(u))
    return Pm, Qm, vidx, dc, nx, [NN + b if b >= 0 else 0 for b in adbr]

def make_devs(spec, adbr):
    """spec: list per region of dicts (or None). keys:
       buck=dict(L, Ron, Ipk, th, dly, pred, band, dwell, taper)  transient injector on the parallel 1.8 V path
       reg=dict(L, Ron, fsw, fc)                                  PWM-regulating buck (DVFS / FIVR mode)
       ad=dict(Ca, G, tau_a, f_lp, Pq)                            Meng/Saleh active decap"""
    D = (Dev * NR)()
    for r in range(NR):
        s = spec[r] if spec else None; d = D[r]
        d.Vin = 1.8; d.L = 2e-9; d.Rs = 1.0; d.t_taper = 30e-9; d.tau_a = 1e-9; d.f_lp = 1e6; d.G = 0.0
        if not s: continue
        if 'buck' in s:
            b = s['buck']; d.buck_mode = 1; d.L = b['L']; d.Vin = b.get('Vin', 1.8); d.Rs = b['Ron'] + b.get('Rdcr', 0.3 * b['L'] / 1e-9)
            d.Ipk = b['Ipk']; d.band = b.get('band', 0.25 * b['Ipk']); d.th_b = b['th']; d.dly_b = int(b.get('dly', 4))
            d.use_pred_b = int(b.get('pred', 0)); d.dwell = int(b.get('dwell', 4)); d.t_taper = b.get('taper', 30e-9)
            d.Esw = 5e-12 / b['Ron']
        if 'reg' in s:
            g = s['reg']; d.buck_mode = 2; d.L = g['L']; d.Rs = g['Ron'] + 0.3 * g['L'] / 1e-9; d.fsw = g['fsw']; d.fc_b = g['fc']
            d.Esw = 5e-12 / g['Ron']; d.Cr = g.get('Cr', 2e-9); d.no_int = int(g.get('no_int', 0)); d.Vin = g.get('Vin', 1.8)
        if 'ad' in s:
            a = s['ad']; d.ad_on = 1; d.G = a['G']; d.tau_a = a.get('tau_a', 1e-9); d.f_lp = a.get('f_lp', 1e6)
            d.Pq = a.get('Pq', 1.8e-3 * a['G'] / 4); d.ad_br = int(adbr[r])
    return D

def area(p, spec):
    a = dict(decap=float(np.sum(p['A_x'])), active=0.0, reg=0.0)
    for s in (spec or []):
        if not s: continue
        if 'ad' in s: a['active'] += s['ad']['Ca'] / HV_DENS
        for key in ('buck', 'reg'):
            if key in s: a['reg'] += IND_MM2_PER_NH * s[key]['L'] / 1e-9 * s[key].get('phases', 1)
    a['total'] = a['decap'] + a['active'] + a['reg']; return a

def predict_triggers(Ihat, th_mV, lead, g=None):
    """per-region boolean trigger: forecast die-voltage deviation from the issue-stream estimate (mA) below -th within lead cycles"""
    if g is None:
        import replay3; g = replay3.nominal_g()
    tot = Ihat.sum(1).astype(float); x = tot - np.convolve(tot, np.ones(256) / 256, 'same')
    v = np.convolve(x, g)[:len(x)]
    fut = np.full(len(v), 0.0)
    for k in range(1, lead + 1): fut[:-k] = np.minimum(fut[:-k], v[k:])
    trig = (fut < -th_mV)
    return np.repeat(trig[:, None], NR, 1).astype(np.uint8)

def simulate(I_mA, p=None, spec=None, pred=None, vref=None, warm=2000, want_v=False):
    """I_mA: (ncyc, 4) region currents. Returns metrics dict (droop in mV relative to VNOM)."""
    p = dict(DEFAULT, **(p or {}))
    Ca = [(spec[r]['ad']['Ca'] if (spec and spec[r] and 'ad' in spec[r]) else 0.0) for r in range(NR)]
    Pm, Qm, vidx, dc, nx, adbr = build(p, Ca)
    I_mA = np.ascontiguousarray(np.vstack([np.repeat(I_mA[:1], warm, 0), I_mA]), np.float32); ncyc = len(I_mA)
    Im = I_mA[warm:warm + 3000].mean(0) * 1e-3
    v0 = VNOM if vref is None else float(vref[0])
    if p.get('fivr_only'):      # die fed only by on-die regulators: start from the connected network's operating point
        _, _, _, dcn, _, _ = build(dict(p, fivr_only=False), Ca); x0 = dcn(v0, Im); vr0 = v0
    else:
        x0 = dc(v0, Im); vd = x0[vidx[:NR]].mean(); vr0 = v0 + (v0 - vd); x0 = dc(vr0, Im)
    D = make_devs(spec, adbr)
    pr = np.zeros((ncyc, NR), np.uint8) if pred is None else np.ascontiguousarray(np.vstack([np.zeros((warm, NR), np.uint8), pred]), np.uint8)
    vrc = None if vref is None else np.ascontiguousarray(np.concatenate([np.full(warm, vref[0]), vref]), np.float64)
    vmin = np.zeros((ncyc, NR), np.float32); vavg = np.zeros((ncyc, NR), np.float32); st = np.zeros(9)
    x0 = np.ascontiguousarray(x0, np.float64)
    f = LIB.run; f.restype = ctypes.c_int
    dp = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    rc = f(ctypes.c_int(ncyc), ctypes.c_int(SUB), ctypes.c_int(nx), ctypes.c_int(NR), dp(np.ascontiguousarray(Pm)), dp(np.ascontiguousarray(Qm)),
           vidx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), I_mA.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
           dp(np.ascontiguousarray(LEAK_MA, np.float64)), dp(np.ascontiguousarray(SHAPE, np.float64)),
           pr.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)), D,
           (dp(vrc) if vrc is not None else ctypes.POINTER(ctypes.c_double)()), ctypes.c_double(VNOM), ctypes.c_double(p['fc_vr']),
           *([ctypes.c_double(p['slew_vr'])] if NEWSIG else []), ctypes.c_double(DT), ctypes.c_double(vr0), dp(x0), vmin.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
           vavg.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), dp(st))
    vmin = vmin[warm:]; vavg = vavg[warm:]
    ref = VNOM if vref is None else np.asarray(vref)[:, None]
    dr = (ref - vmin) * 1e3; dw = dr.max(1)
    Eload = st[1]
    m = dict(p99=float(np.percentile(dw, 99)), p999=float(np.percentile(dw, 99.9)), max=float(dw.max()),
             reg_p999=[float(np.percentile(dr[:, r], 99.9)) for r in range(NR)],
             over_max=float(((vavg - ref) * 1e3).max()),
             loss_buck=float(st[2] - st[3] + st[4]) / Eload if Eload > 0 else 0.0, loss_ad=float(st[5] + st[6]) / Eload if Eload > 0 else 0.0,
             buck_duty=st[7] / (ncyc * SUB), f_sw_MHz=st[8] / (ncyc * T_CYC) / 1e6, E_load=Eload, mean_v=float(vavg.mean()))
    if want_v: m['vmin'] = vmin; m['vavg'] = vavg
    return m
