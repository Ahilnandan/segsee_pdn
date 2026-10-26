"""2-D on-die power grid with distributed micro-regulators (Pulsim, switching level).

Die = 4x4 tiles (replay2 tiles=True floorplan). Each tile: Cdie/16 to ground, its own load current source,
its own package bump; neighbouring tiles joined by grid segments (Rg, Lg). Package: common plane segment from
the off-chip bulk + per-tile bumps, split so that a stiff grid reproduces the lumped pdnchain package (507 pH,
79.6 mOhm) exactly. Off-chip 4-phase buck carries DC (its loop nulls the FIVR's average current, as case 'D').
FIVR: N phases from 1.8 V, natural-sampled double-edge PWM, type-3 control. Placement:
  'central'     all phases feed one tile (conventional single-point regulator)
  'distributed' phase k feeds tile k mod 16 (micro-regulator slices spread over the die, Intel-style bridges)
Control:
  'global'  one compensator on the mean tile voltage drives every phase
  'local'   each phase has its own compensator on ITS tile voltage (fully distributed regulation)
Per-core-share scaling as pdnchain (same totals: inductance, switch area, decap)."""
import numpy as np
import pulsim as p
from pdnchain import _buck, PI, VNOM, L_PKG, R_PKG, OFF_L, OFF_DCR, OFF_C, OFF_ESR, OFF_FI
import t3design

class T3:
    """type-3: wI/s * prod (1+s/wz)/(1+s/wp), explicit delay line (same realisation as pdnchain ctl='t3')"""
    def __init__(s, wI, fz, fp, Td, dt):
        s.wI = wI; s.wz = [2 * np.pi * f for f in fz]; s.wp = [2 * np.pi * f for f in fp]
        s.xi = 0.0; s.x = [0.0, 0.0]; s.dl = [0.0] * max(1, int(round(Td / dt))); s.dp = 0
    def step(s, err, dts, dnom, sat):
        e = s.dl[s.dp]; s.dl[s.dp] = err; s.dp = (s.dp + 1) % len(s.dl)
        xi = s.xi + s.wI * e * dts; u = xi
        for j in range(2):
            r = s.wp[j] / s.wz[j]
            s.x[j] += (1.0 - np.exp(-s.wp[j] * dts)) * (u - s.x[j]); u = r * u + (1.0 - r) * s.x[j]
        c = dnom + u
        if not ((c > sat and e > 0) or (c < 0.0 and e < 0)): s.xi = xi
        return min(sat, max(0.0, c))

def simulate_grid(It_mA, place='distributed', ctl='global', Rg=0.2, Lg=0.3e-9, fsw=35e6, nph=16, L_ph=448e-9,
                  ron=0.2, dcr=None, cdie=5e-9, plane_frac=0.5, dt=50e-12, cyc_s=1e-9, warmup=3000, wI=None, fz=1e6):
    It = np.asarray(It_mA, float) * 1e-3
    It = np.concatenate([np.repeat(It[:1], warmup, 0), It]); n = len(It); NT = 16
    dcr = 0.2 * L_ph / 112e-9 if dcr is None else dcr
    b = p.CircuitBuilder()
    b.add_voltage_source('V12', 'v12', 'gnd', 12.0); b.add_voltage_source('V18', 'v18', 'gnd', 1.8)
    i0 = float(It[0].sum()); nph_off = 4; sw = []
    vb0 = VNOM + R_PKG * i0
    sw += _buck(b, 'o', 'v12', 'vb', nph_off, OFF_L, 2e-3, OFF_DCR, i0)
    b.add_resistor('Resr', 'vb', 'cb', OFF_ESR); b.add_capacitor('Cb', 'cb', 'gnd', OFF_C); b.set_initial('Cb', vb0)
    # package: common plane segment + 16 bumps (parallel combination == lumped package when the grid is stiff)
    Lc, Rc = plane_frac * L_PKG, plane_frac * R_PKG
    Lb, Rb = (1 - plane_frac) * L_PKG * NT, (1 - plane_frac) * R_PKG * NT
    b.add_inductor('Lpl', 'vb', 'vpm', Lc, i0=i0); b.add_resistor('Rpl', 'vpm', 'vpk', Rc)
    tn = [f't{k}' for k in range(NT)]
    for k in range(NT):
        b.add_inductor(f'Lb{k}', 'vpk', f'bm{k}', Lb, i0=i0 / NT); b.add_resistor(f'Rb{k}', f'bm{k}', tn[k], Rb)
        b.add_capacitor(f'C{k}', tn[k], 'gnd', cdie / NT); b.set_initial(f'C{k}', VNOM)
        b.add_current_source(f'I{k}', tn[k], 'gnd', 0.0)
    for r in range(4):
        for c in range(4):
            k = 4 * r + c
            for k2 in ((k + 1) if c < 3 else None, (k + 4) if r < 3 else None):
                if k2 is None: continue
                b.add_inductor(f'Lg{k}_{k2}', tn[k], f'g{k}_{k2}', Lg, i0=0.0); b.add_resistor(f'Rg{k}_{k2}', f'g{k}_{k2}', tn[k2], Rg)
    # FIVR phases
    tgt = [5] * nph if place == 'central' else [k % NT for k in range(nph)]
    for k in range(nph):
        sw += _buck(b, f'f{k}_', 'v18', f'fo{k}', 1, L_ph, ron, dcr, 0.0)
        b.add_resistor(f'Rgr{k}', f'fo{k}', tn[tgt[k]], 20e-3 * nph)
    svn = list(b.state_var_names()); swi = {nm: b.switch_index_of(nm) for nm in sw}
    tid = [b.node_id_of(t) for t in tn]; vbid = b.node_id_of('vb')
    Lf = [svn.index(f'I(f{k}_L0)') for k in range(nph)]; Lo = [svn.index(f'I(oL{k})') for k in range(nph_off)]
    nunk, nsw = b.graph.num_branches, b.graph.num_switches
    # controller design on the lumped-equivalent parallel plant (same totals); conservative crossover <= fsw/1.8
    if wI is None:
        d = t3design.best_par(fsw / 1.8, L_ph=L_ph, nph=nph, ron=ron, dcr=dcr, C=cdie)
        wI, fz = d['wI'], d['fz']
    ncl = 1 if ctl == 'global' else nph
    wIk = wI   # each local loop drives 1/N of the current but sees ~1/N of the plant gain: same total loop gain
    C = [T3(wIk, (fz, fz), (600e6, 1.2e9), 0.6e-9, dt) for _ in range(ncl)]
    Tp, Tpo = 1.0 / fsw, 1.0 / 500e3
    Dn = VNOM / 1.8
    st = dict(cf=[Dn] * nph, co=vb0 / 12.0, iref=i0, ILf=0.0, v=[VNOM] * NT, prev=0.0)
    mo = max(1, int(round((Tpo / nph_off) / dt))); ob = dict(buf=None, sum=0.0, p=0)
    mf = max(1, int(round(Tp / dt))); fb = dict(buf=[0.0] * mf, sum=0.0, p=0)
    Rk_off = 2 * np.pi * OFF_FI * (OFF_L / nph_off) / 12.0
    VMIN = np.full((n, NT), 9.0); VMAX = np.full((n, NT), -9.0)
    def obs(t, x):
        v = [float(x[i]) for i in tid]; st['v'] = v
        kk = int(t / cyc_s + 1e-9)
        if kk < n:
            a = np.array(v); np.minimum(VMIN[kk], a, out=VMIN[kk]); np.maximum(VMAX[kk], a, out=VMAX[kk])
        dts = t - st['prev']; st['prev'] = t
        if dts <= 0: return
        # FIVR total inductor current, averaged over one switching period -> split loop
        ilf = sum(float(x[j]) for j in Lf)
        fb['sum'] += ilf - fb['buf'][fb['p']]; fb['buf'][fb['p']] = ilf; fb['p'] = (fb['p'] + 1) % mf; st['ILf'] = fb['sum'] / mf
        ilo = sum(float(x[j]) for j in Lo)
        if ob['buf'] is None: ob['buf'] = [ilo] * mo; ob['sum'] = ilo * mo
        ob['sum'] += ilo - ob['buf'][ob['p']]; ob['buf'][ob['p']] = ilo; ob['p'] = (ob['p'] + 1) % mo
        st['iref'] += 2 * np.pi * 20e3 * st['ILf'] * dts
        st['co'] = min(0.95, max(0.0, float(x[vbid]) / 12.0 + Rk_off * (st['iref'] - ob['sum'] / mo)))
        if ctl == 'global':
            c = C[0].step(VNOM - float(np.mean(v)), dts, Dn, 0.98); st['cf'] = [c] * nph
        else:
            st['cf'] = [C[k].step(VNOM - v[tgt[k]], dts, Dn, 0.98) for k in range(nph)]
    def swfn(t):
        m = p.SwitchStateMask(nsw)
        for k in range(nph_off):
            ph = (t / Tpo + k / nph_off) % 1.0; on = ph < st['co']
            m.set(swi[f'oHS{k}'], on); m.set(swi[f'oLS{k}'], not on)
        for k in range(nph):
            ph = (t / Tp + k / nph) % 1.0; on = st['cf'][k] > 1.0 - abs(2.0 * ph - 1.0)
            m.set(swi[f'f{k}_HS0'], on); m.set(swi[f'f{k}_LS0'], not on)
        return m
    def bex(t):
        k = min(int(t / cyc_s + 1e-9), n - 1); out = [0.0] * nunk
        row = It[k]
        for j in range(NT): out[tid[j]] = row[j]
        return out
    p.simulate(b, t_end=n * cyc_s, dt=dt, b_extra_fn=bex, switch_fn=swfn, step_observer=obs, start_from_dc_op=False, store_every=10**9)
    return dict(vmin=VMIN[warmup:n - 1], vmax=VMAX[warmup:n - 1])
