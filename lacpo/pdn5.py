"""pdn5: tiled (8x8, 250 um) PDN of the MediumBoom tile — builder (sparse descriptor MNA -> trapezoidal), steppable
C++ simulator wrapper (libpdn5.so), frequency-domain impedance, area model.

Network (real units, see docs/protocol-checklist.md Phase 0):
 off-chip buck (averaged, current into board node vo) -> board C1 47 uF (ESR) -> board R_b/L_b -> package node pk
 (package decap 1 uF, ESL/ESR) -> package plane L_plane/R_plane -> bump field bf -> per-tile bump branch (VDD bumps per
 tile) -> tile die node d_t. Die nodes joined by grid resistors R_edge (4-neighbour). Tile: intrinsic C (cell share),
 optional damped decap (area, series R), optional charge-pump bank + LDO (reservoir node rn_t with C_res), optional
 on-die buck. Loads per tile from tilecurrents / closed-loop replay."""
import ctypes, json, os, numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spl
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
LIB = ctypes.CDLL(HERE + os.environ.get('PDN5_LIB', 'libpdn5.so'))
FP = json.load(open(HERE + 'floorplan.json')); N = FP['n']; T = N * N
VNOM = 1.1; T_CYC = 3e-9; SUB = 12; DT = T_CYC / SUB
DECAP_DENS = 12e-9; CP_DENS = 3e-9; DIE_MM2 = FP['die_mm2']
_sh = np.array([0.6, 1.6, 1.9, 1.6, 1.3, 1.1, 0.9, 0.8, 0.7, 0.7, 0.7, 0.7]); SHAPE = _sh / _sh.mean()
LEAK_T = np.load(HERE + 'tile_leak_mA.npy')
CELL_SHARE = LEAK_T / LEAK_T.sum()

class Tile(ctypes.Structure):
    _fields_ = [('bank', ctypes.c_int)] + [(n, ctypes.c_double) for n in ('R_pass', 'gm', 'tau', 'off', 'I_ch', 'eta', 'Vmul')] + \
               [('ridx', ctypes.c_int), ('buck', ctypes.c_int)] + [(n, ctypes.c_double) for n in ('Vin', 'L', 'Rs', 'Ipk', 'band', 'th', 't_taper', 'fc', 'Cr')] + \
               [('phases', ctypes.c_int), ('dwell', ctypes.c_int)]

DEFAULT = dict(Vin1=3.3, L1=2.2e-6, R1=40e-3, C1=47e-6, R_c1=3e-3, fc_off=100e3, slew=20e-3 / 1e-6, f_share=20e3,
               L_b=2e-9, R_b=3e-3, C_pk=1e-6, L_pk=100e-12, R_pk=20e-3, L_plane=0.45e-9, R_plane=40e-3,
               vdd_bumps=40, L_bump=0.25e-9, R_bump=20e-3, R_edge=0.08, c_intr=3e-9,
               A_x=None, R_d=None, bank=None, buck=None, tau_x=0.1e-9, tau_intr=0.3e-9, V_cp0=2.15)

def build(p):
    p = dict(DEFAULT, **p); A_x = p['A_x'] if p['A_x'] is not None else np.zeros(T); R_d = p['R_d'] if p['R_d'] is not None else np.zeros(T)
    bank = p['bank'] or {}; nodes = ['vo', 'voc', 'pk', 'pkc', 'bf'] + ['d%d' % t for t in range(T)] + ['ci%d' % t for t in range(T)]
    for t in range(T):
        if A_x[t] > 0: nodes.append('cx%d' % t)
    for t in bank: nodes += ['rn%d' % t, 'rc%d' % t]
    ni = {n: i for i, n in enumerate(nodes)}; NN = len(nodes)
    Cn = np.full(NN, 1e-12); G = sp.lil_matrix((NN, NN)); inds = []
    def res(a, b, R):
        g = 1.0 / R; ia, ib = ni[a], ni[b]; G[ia, ia] += g; G[ib, ib] += g; G[ia, ib] -= g; G[ib, ia] -= g
    res('vo', 'voc', p['R_c1']); Cn[ni['voc']] += p['C1']
    inds.append(('vo', 'pk', p['L_b'], p['R_b'])); inds.append(('pk', 'pkc', p['L_pk'], p['R_pk'])); Cn[ni['pkc']] += p['C_pk']
    inds.append(('pk', 'bf', p['L_plane'], p['R_plane']))
    nb = p['vdd_bumps'] / T
    for t in range(T):
        inds.append(('bf', 'd%d' % t, p['L_bump'] / nb, p['R_bump'] / nb))
        ci = p['c_intr'] * CELL_SHARE[t]; Cn[ni['ci%d' % t]] += ci; res('d%d' % t, 'ci%d' % t, p['tau_intr'] / ci)   # intrinsic C with its ESR
        if A_x[t] > 0:
            cx = A_x[t] * DECAP_DENS; Cn[ni['cx%d' % t]] += cx; res('d%d' % t, 'cx%d' % t, p['tau_x'] / cx + R_d[t])
        r, c = divmod(t, N)
        if c + 1 < N: res('d%d' % t, 'd%d' % (t + 1), p['R_edge'])
        if r + 1 < N: res('d%d' % t, 'd%d' % (t + N), p['R_edge'])
    for t, b in bank.items():
        cr = b['C_res']; Cn[ni['rc%d' % t]] += cr; res('rn%d' % t, 'rc%d' % t, max(p['tau_x'] / cr, 1e-3))
    NL = len(inds); nx = NN + NL; nu = 1 + 3 * T
    E = sp.lil_matrix((nx, nx)); A = sp.lil_matrix((nx, nx)); Bm = sp.lil_matrix((nx, nu))
    E[:NN, :NN] = sp.diags(Cn); A[:NN, :NN] = -G
    for k, (a, b, L, R) in enumerate(inds):
        j = NN + k; E[j, j] = L; A[j, j] = -R; A[j, ni[a]] += 1; A[ni[a], j] -= 1; A[j, ni[b]] -= 1; A[ni[b], j] += 1
    Bm[ni['vo'], 0] = 1
    for t in range(T):
        Bm[ni['d%d' % t], 1 + t] = -1; Bm[ni['d%d' % t], 1 + T + t] = 1
        if t in bank: Bm[ni['rn%d' % t], 1 + 2 * T + t] = -1
    E = E.tocsc(); A = A.tocsc(); Bm = Bm.tocsc()
    return dict(p=p, E=E, A=A, B=Bm, ni=ni, nx=nx, nu=nu, NN=NN, didx=np.array([ni['d%d' % t] for t in range(T)], np.int32))

def dc(net, v_die, I_A):
    """operating point: die mean at v_die with load I_A (per tile), off-chip current = sum; reservoirs at V_cp0."""
    A, B, ni, NN = net['A'].tolil(), net['B'], net['ni'], net['NN']
    for i in range(NN): A[i, i] -= 1e-12
    u = np.zeros(net['nu']); u[0] = I_A.sum(); u[1:1 + T] = I_A
    vo_set = v_die + 0.02
    for _ in range(3):
        A2 = A.copy(); rhs = -B.dot(u); iv = ni['vo']; A2[iv, :] = 0; A2[iv, iv] = 1; rhs[iv] = vo_set
        x = spl.spsolve(A2.tocsc(), rhs); vo_set += v_die - x[net['didx']].mean()
    for k in ni:
        if k.startswith('rn') or k.startswith('rc'): x[ni[k]] = net['p']['V_cp0']
    return x

def impedance(net, f, weights=None):
    """|Z(f)| from a load distributed by `weights` (default: cell share) to the weighted mean die voltage."""
    w = CELL_SHARE if weights is None else weights; E, A, B = net['E'], net['A'], net['B']
    u = np.zeros(net['nu']); u[1:1 + T] = -w; b = B.dot(u); out = []
    for fk in f:
        M = (2j * np.pi * fk) * E - A; x = spl.spsolve(M.tocsc(), b); out.append(abs(x[net['didx']].dot(w)))
    return np.array(out)

class Sim(object):
    def __init__(self, p, Im_A=None):
        net = build(p); self.net = net; P = net['p']
        Ml = (net['E'] - net['A'] * DT / 2).tocoo(); Mr = (net['E'] + net['A'] * DT / 2).tocoo(); Bq = (net['B'] * DT / 2).tocoo()
        Im_A = Im_A if Im_A is not None else LEAK_T * 1e-3 * 3
        x0 = dc(net, VNOM, Im_A)
        tiles = (Tile * T)()
        for t in range(T):
            tl = tiles[t]; tl.ridx = -1; tl.tau = 1e-9; tl.R_pass = 1.0; tl.L = 2e-9; tl.phases = 1; tl.t_taper = 30e-9
            if P['bank'] and t in P['bank']:
                b = P['bank'][t]; tl.bank = 1; tl.R_pass = b['R_pass']; tl.gm = b['gm']; tl.tau = b['tau']; tl.off = b['off']
                tl.I_ch = b['I_ch']; tl.eta = b.get('eta', 0.9); tl.Vmul = b.get('Vmul', 2.0); tl.ridx = net['ni']['rn%d' % t]
            if P['buck'] and t in P['buck']:
                g = P['buck'][t]; tl.buck = int(g['role']); tl.Vin = g['Vin']; tl.L = g['L']; tl.Rs = g['Ron'] + 0.010 + 0.002 * g['L'] / 1e-9
                tl.Ipk = g.get('Ipk', 0.05); tl.band = 0.25 * tl.Ipk; tl.th = g.get('th', 0.04); tl.dwell = int(g.get('dwell', 2))
                tl.fc = g.get('fc', 20e6); tl.Cr = P['c_intr'] * CELL_SHARE[t] + (P['A_x'][t] if P['A_x'] is not None else 0) * DECAP_DENS + 1e-12
                tl.phases = int(g.get('phases', 1))
        ip = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_int)); dp = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        self._keep = []
        def arr(m):
            i, j, v = np.ascontiguousarray(m.row, np.int32), np.ascontiguousarray(m.col, np.int32), np.ascontiguousarray(m.data, np.float64)
            self._keep += [i, j, v]; return ctypes.c_int(len(v)), ip(i), ip(j), dp(v)
        x0 = np.ascontiguousarray(x0, np.float64); self._keep += [x0, net['didx']]
        LIB.sim_new.restype = ctypes.c_void_p
        self.h = LIB.sim_new(ctypes.c_int(net['nx']), ctypes.c_int(T), ctypes.c_int(SUB), ctypes.c_double(DT), ctypes.c_double(VNOM),
                             *arr(Ml), *arr(Mr), *arr(Bq), ctypes.c_int(net['nu']), ip(net['didx']), ctypes.c_int(net['ni']['vo']), dp(x0), tiles,
                             ctypes.c_double(P['Vin1']), ctypes.c_double(P['L1']), ctypes.c_double(P['R1']), ctypes.c_double(P['fc_off']),
                             ctypes.c_double(P['C1'] + P['C_pk']), ctypes.c_double(P['slew']), ctypes.c_double(P['f_share']), ctypes.c_double(float(Im_A.sum())))
        self.tiles = tiles; self.stats = np.zeros(6)
        assert LIB.sim_info(ctypes.c_void_p(self.h)) == 0, 'sparse LU failed'
    def run(self, I_mA, vref=None, role_off=None, role_on=None, arm=None, want_res=False):
        I = np.ascontiguousarray(I_mA, np.float32); n = len(I)
        vmin = np.zeros((n, T), np.float32); vavg = np.zeros((n, T), np.float32); vres = np.zeros((n, T), np.float32) if want_res else None
        f = lambda a, t: None if a is None else np.ascontiguousarray(a, t)
        vr, ro, rn, ar = f(vref, np.float64), f(role_off, np.uint8), f(role_on, np.uint8), f(arm, np.uint8)
        P = lambda a, c: a.ctypes.data_as(ctypes.POINTER(c)) if a is not None else ctypes.POINTER(c)()
        LIB.sim_run(ctypes.c_void_p(self.h), ctypes.c_int(n), P(I, ctypes.c_float), P(np.ascontiguousarray(LEAK_T, np.float64), ctypes.c_double),
                    P(np.ascontiguousarray(SHAPE, np.float64), ctypes.c_double), P(vr, ctypes.c_double), P(ro, ctypes.c_ubyte), P(rn, ctypes.c_ubyte),
                    P(ar, ctypes.c_ubyte), P(vmin, ctypes.c_float), P(vavg, ctypes.c_float), P(vres, ctypes.c_float), P(self.stats, ctypes.c_double))
        return vmin, vavg, vres
    def __del__(self):
        try: LIB.sim_free(ctypes.c_void_p(self.h))
        except Exception: pass

def area(p):
    a = float(np.sum(p.get('A_x') if p.get('A_x') is not None else 0))
    for t, b in (p.get('bank') or {}).items(): a += b['C_res'] / CP_DENS + 0.002 / b['R_pass'] + 0.004
    for t, g in (p.get('buck') or {}).items(): a += 0.005 + g.get('phases', 1) * 0.002 / g['Ron']
    return a
