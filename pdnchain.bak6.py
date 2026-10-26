"""Real switching power-delivery chains in Pulsim, driven by per-cycle core current
traces (mA per 1 ns cycle) from the validated replay (replay2.run()['In']).

Absolute voltages (V_core nominal 0.9 V). Every converter is a real switched
synchronous buck (ideal switches with Ron, real inductors/caps), PWM generated
per phase from a carrier, closed loop with a discrete PI.

Component values are a PER-CORE SHARE of a real multi-core rail (a VR built for
tens of amps, of which this one BOOM core draws ~0.03-0.43 A): inductances
scaled up / capacitances scaled down by the same factor so ripple ratios and
time constants match a real design. The die-side model (507 pH / 79.6 mOhm /
5 nF) is identical to every earlier result in this project so case A is
directly comparable to replay2.droop().

  case 'A'  off-chip high-current buck direct to the core rail:
            12V --4ph buck 500kHz--> bulk C --pkg L,R--> vdie (5nF) --> load
  case 'B'  off-chip multiphase feeding an on-die FIVR:
            12V --4ph buck 500kHz--> 1.8V bulk --pkg L,R--> v18 (C_in)
                --4ph interleaved FIVR (fsw, L_ph)--> R_grid --> vdie (5nF) --> load

Optional (novelty test, case B): scheduler->FIVR current feed-forward. The
issue stage knows the current shape of every uop it issues 1-3 cycles before
execute (power-token ROM, BACT_CAND amp/dur). `ff_lead` cycles of that
prediction are fed into the FIVR duty as an inductor-current target, so the
FIVR starts ramping before vdie moves instead of after.
"""
import numpy as np
import pulsim as p

VNOM = 0.9
L_PKG, R_PKG, C_DIE = 507e-12, 79.6e-3, 5e-9
# off-chip VR: per-core share (K=250) of a real 100 A 4-phase VR -- 150 nH / 0.5 mOhm per phase,
# 2 mF bulk with 0.5 mOhm effective ESR. Scaling L, R up and C down by the SAME K keeps f0, Q and
# ripple ratio of the real design (an earlier inconsistent choice gave Q~29 and a slow instability).
OFF_K = 250.0
OFF_L, OFF_DCR, OFF_C, OFF_ESR = 150e-9 * OFF_K, 0.5e-3 * OFF_K, 2e-3 / OFF_K, 0.5e-3 * OFF_K
OFF_FI = 100e3    # off-chip INNER (average-current) loop crossover ~ fsw/5
OFF_FV = 30e3     # off-chip OUTER voltage loop crossover (real multiphase VRs: tens-100s kHz)

def _buck(b, tag, vin_node, vout_node, nph, L, ron, rdcr, i0_total):
    sws = []
    for k in range(nph):
        sw = f'{tag}sw{k}'; mid = f'{tag}m{k}'
        b.add_switch(f'{tag}HS{k}', vin_node, sw, g_on=1.0/ron, g_off=1e-9); sws.append(f'{tag}HS{k}')
        b.add_switch(f'{tag}LS{k}', sw, 'gnd', g_on=1.0/ron, g_off=1e-9);    sws.append(f'{tag}LS{k}')
        b.add_inductor(f'{tag}L{k}', sw, mid, L, i0=i0_total/nph)
        b.add_resistor(f'{tag}R{k}', mid, vout_node, rdcr)
    return sws

class PI:
    def __init__(s, kp, ki, lo, hi, i0):
        s.kp, s.ki, s.lo, s.hi, s.i = kp, ki, lo, hi, i0
    def step(s, e, dt):
        ni = s.i + s.ki * e * dt
        out = s.kp * e + ni
        if (out > s.hi and e > 0) or (out < s.lo and e < 0):
            out = s.kp * e + s.i
        else:
            s.i = ni
        return min(s.hi, max(s.lo, out))

def simulate_chain(I_mA, case='A', pred_mA=None, ffsig_mA=None, fsw_fivr=140e6, nph_fivr=4, L_fivr=4e-9, c_in18=10e-9,
                   ff_lead=0, ff_gain=1.0, cyc_s=1e-9, dt=50e-12, fivr_bw=None, ideal_in=False, ctl='vm', fci_div=8.0, fcv_div=20.0, kp_vm=0.02, fint=2e6, warmup=1500, cdie=C_DIE, bank=None, sched=None, nmax=None, sig_gain=1.0, sig_mode='both',
                   pwm='saw', ron_f=20e-3, dcr_f=30e-3, t3=None):
    if sched is not None:
        first = float(next(sched)); I_mA = np.full(nmax, first)     # filled cycle-by-cycle by the scheduler
    I_A = np.asarray(I_mA, float) * 1e-3
    P_A = I_A if pred_mA is None else np.asarray(pred_mA, float) * 1e-3
    # warm-up prefix at the first sample's current, discarded from the outputs (start-up transient)
    I_A = np.concatenate([np.full(warmup, I_A[0]), I_A]); P_A = np.concatenate([np.full(warmup, P_A[0]), P_A])
    FF = None if ffsig_mA is None else np.concatenate([np.zeros(warmup), np.asarray(ffsig_mA, float) * 1e-3])
    n = len(I_A); i0 = I_A[0]
    b = p.CircuitBuilder()
    b.add_voltage_source('V12', 'v12', 'gnd', 12.0)
    fsw_off, nph_off = 500e3, 4
    sw_names = []
    hasF = case in ('B', 'D')
    # case 'D' (PARALLEL / hybrid supply): the off-chip buck feeds the die DIRECTLY through the package
    # (carries the DC load current, no double conversion) AND an on-die FIVR from the 1.8 V rail feeds
    # the same die node (fast loop owns the voltage). The off-chip loop nulls the FIVR's average current.
    if case in ('A', 'D'):
        L_off, Cb, esr, dcr_off = OFF_L, OFF_C, OFF_ESR, OFF_DCR
        v_bulk0 = VNOM + R_PKG * i0
        sw_names += _buck(b, 'o', 'v12', 'vb', nph_off, L_off, 2e-3, dcr_off, i0)
        b.add_resistor('Resr', 'vb', 'cb', esr); b.add_capacitor('Cb', 'cb', 'gnd', Cb); b.set_initial('Cb', v_bulk0)
        b.add_inductor('Lpkg', 'vb', 'vp', L_PKG, i0=i0); b.add_resistor('Rpkg', 'vp', 'vdie', R_PKG)
        b.add_capacitor('Cdie', 'vdie', 'gnd', cdie); b.set_initial('Cdie', VNOM)
        D_off0 = (v_bulk0 + i0 * (2e-3 + dcr_off) / nph_off) / 12.0
        ctl_off = PI(2 * np.pi * OFF_FV * Cb, 2 * np.pi * OFF_FV * Cb * 2 * np.pi * OFF_FV / 2, -50, 50, i0)
        sense_off = 'vdie'; vref_off = VNOM; vout_off0 = v_bulk0
    if case == 'D':
        b.add_voltage_source('V18', 'v18', 'gnd', 1.8)
        rgrid = 20e-3
        sw_names += _buck(b, 'f', 'v18', 'vf', nph_fivr, L_fivr, ron_f, dcr_f, 0.0)
        b.add_resistor('Rgrid', 'vf', 'vdie', rgrid)
        D_f0 = VNOM / 1.8; Dnom_f = D_f0; Ki_vm = 2 * np.pi * fint / 1.8; L_eff = L_fivr / nph_fivr
        st_split = {'iref': i0}
        if ctl == 't3':
            import t3design
            T3 = dict(fc=70e6, fz=(3e6, 3e6), fp=(600e6, 1.2e9), Td=0.6e-9); T3.update(t3 or {})
            if 'wI' in T3:
                dsg = t3design.design_par_gain(T3['wI'], fz=T3['fz'], fp=T3['fp'], Td=T3['Td'], L_ph=L_fivr, nph=nph_fivr,
                                               C=cdie, ron=ron_f, dcr=dcr_f, rgrid=rgrid)
            else:
                dsg = t3design.design_par(fc=T3['fc'], fz=T3['fz'], fp=T3['fp'], Td=T3['Td'], L_ph=L_fivr, nph=nph_fivr,
                                          C=cdie, ron=ron_f, dcr=dcr_f, rgrid=rgrid)
            dsg = {k: v for k, v in dsg.items() if not isinstance(v, list)}
            t3s = dict(wI=dsg['wI'], wz=[2 * np.pi * f for f in T3['fz']], wp=[2 * np.pi * f for f in T3['fp']],
                       xi=0.0, x=[0.0, 0.0], dl=[0.0] * max(1, int(round(T3['Td'] / dt))), dp=0, design=dsg)
    if case == 'B':
        eff = 0.9
        i18 = i0 * VNOM / 1.8 / eff
        L_off, Cb, esr, dcr_off = OFF_L, OFF_C, OFF_ESR, OFF_DCR
        v18b0 = 1.8 + R_PKG * i18
        if ideal_in:
            b.add_voltage_source('V18', 'v18', 'gnd', 1.8)
        else:
            sw_names += _buck(b, 'o', 'v12', 'vb', nph_off, L_off, 2e-3, dcr_off, i18)
            b.add_resistor('Resr', 'vb', 'cb', esr); b.add_capacitor('Cb', 'cb', 'gnd', Cb); b.set_initial('Cb', v18b0)
            b.add_inductor('Lpkg', 'vb', 'vp', L_PKG, i0=i18); b.add_resistor('Rpkg', 'vp', 'v18', R_PKG)
            b.add_capacitor('Cin', 'v18', 'gnd', c_in18); b.set_initial('Cin', 1.8)
        rgrid = 20e-3
        sw_names += _buck(b, 'f', 'v18', 'vf', nph_fivr, L_fivr, ron_f, dcr_f, i0)
        b.add_resistor('Rgrid', 'vf', 'vdie', rgrid)
        b.add_capacitor('Cdie', 'vdie', 'gnd', cdie); b.set_initial('Cdie', VNOM)
        D_off0 = (v18b0 + i18 * (2e-3 + dcr_off) / nph_off) / 12.0
        ctl_off = PI(2 * np.pi * OFF_FV * Cb, 2 * np.pi * OFF_FV * Cb * 2 * np.pi * OFF_FV / 2, -50, 50, i18)
        sense_off = 'v18'; vref_off = 1.8; vout_off0 = v18b0
        D_f0 = (VNOM + i0 * (rgrid + (ron_f + dcr_f) / nph_fivr)) / 1.8
        # CURRENT-MODE control (inner inductor-current loop damps the high-Q output LC;
        # a voltage-mode PI oscillated at ~55 MHz against f0~71 MHz, Q~14 -- see docs).
        L_eff = L_fivr / nph_fivr
        fci = fsw_fivr / fci_div                   # inner current-loop crossover
        fcv = fivr_bw or fsw_fivr / fcv_div        # outer voltage-loop crossover
        Rk = 2 * np.pi * fci * L_eff / 1.8         # duty per amp of current error
        Kpv = 2 * np.pi * fcv * cdie              # amps per volt
        ctl_f = PI(Kpv, Kpv * 2 * np.pi * fcv / 4.0, -5.0, 5.0, i0)   # outputs inductor-current reference
        Dnom_f = D_f0; Ki_vm = 2 * np.pi * fint / 1.8
        if ctl == 't3':
            # faithful FIVR loop (Burton APEC'14 class): type-3 compensator designed on the averaged
            # plant (t3design.py), wI solved for the target UGF; explicit loop delay Td on the sensed vdie
            import t3design
            T3 = dict(fc=78e6, fz=(10e6, 10e6), fp=(400e6, 800e6), Td=0.6e-9); T3.update(t3 or {})
            dsg = t3design.design(fc=T3['fc'], fz=T3['fz'], fp=T3['fp'], Td=T3['Td'], L_ph=L_fivr, nph=nph_fivr,
                                  C=cdie, ron=ron_f, dcr=dcr_f, rgrid=rgrid)
            t3s = dict(wI=dsg['wI'], wz=[2 * np.pi * f for f in T3['fz']], wp=[2 * np.pi * f for f in T3['fp']],
                       xi=0.0, x=[0.0, 0.0], dl=[0.0] * max(1, int(round(T3['Td'] / dt))), dp=0, design=dsg)
    b.add_current_source('Iload', 'vdie', 'gnd', 0.0)
    if bank:
        vb0 = VNOM + bank.get('head', 0.25)
        b.add_capacitor('Cbank', 'bank', 'gnd', bank['cb'] * 1e-9); b.set_initial('Cbank', vb0)
        b.add_switch('SWrel', 'bank', 'vdie', g_on=1.0 / bank.get('ron', 0.05), g_off=1e-9); sw_names.append('SWrel')
        rsrc = 'v18' if hasF else 'v12'      # recharge drawn from the FIVR input rail (B) / 12 V (A)
        CP = bank.get('circuit')
        if CP:
            # CIRCUIT-LEVEL storing path: 3:2 series-parallel switched-capacitor pump lifting overshoot charge
            # from vdie into the bank.  ph1: Cf1,Cf2 in series across vdie-gnd (absorb);  ph2: Cf1||Cf2 stacked
            # on vdie -> top at vdie+vdie/2 ~1.35 V > vbank (deliver).  Real switches with Ron.
            cf = CP.get('cf', 0.1) * 1e-9; rp = CP.get('ron', 0.5)
            for j in (1, 2):
                b.add_capacitor(f'Cf{j}', f'pa{j}', f'pb{j}', cf); b.set_initial(f'Cf{j}', VNOM / 2)
            PS1 = [('P1a', 'vdie', 'pa1'), ('P1b', 'pb1', 'pa2'), ('P1c', 'pb2', 'gnd')]
            PS2 = [('P2a', 'pa1', 'bank'), ('P2b', 'pa2', 'bank'), ('P2c', 'pb1', 'vdie'), ('P2d', 'pb2', 'vdie')]
            for nm, a_, b_ in PS1 + PS2:
                b.add_switch(nm, a_, b_, g_on=1.0 / rp, g_off=1e-9); sw_names.append(nm)
            # real recharge path: switch + resistor from the 1.8 V rail (replaces the ideal 25 mA source)
            b.add_resistor('Rrch', rsrc, 'rchm', CP.get('rch', (1.8 - vb0) / (bank.get('ir', 25.0) * 1e-3)))
            b.add_switch('SWrch', 'rchm', 'bank', g_on=1.0 / 0.5, g_off=1e-9); sw_names.append('SWrch')

    idx = {nm: b.node_id_of(nm) for nm in ('vdie', sense_off)}
    nunk = b.graph.num_branches; nsw = b.graph.num_switches
    # switch indices in insertion order
    sw_index = {nm: b.switch_index_of(nm) for nm in sw_names}
    svn = list(b.state_var_names())
    Lidx = [svn.index(f'I(fL{k})') for k in range(nph_fivr)] if hasF else []
    Oidx = [svn.index(f'I(oL{k})') for k in range(nph_off)] if not (case == 'B' and ideal_in) else []
    vbid = b.node_id_of('vb') if Oidx else None
    Rk_off = 2 * np.pi * OFF_FI * (OFF_L / nph_off) / 12.0
    mo = max(1, int(round((1.0 / 500e3 / nph_off) / dt)))   # average over the 4-phase ripple period
    st_o = {'buf': [None] * mo, 'sum': 0.0, 'p': 0, 'init': False}
    v18id = b.node_id_of('v18') if hasF else None
    bkid = b.node_id_of('bank') if bank else None
    rsid = b.node_id_of(rsrc) if bank else None
    CP = bank.get('circuit') if bank else None
    if CP:
        cpn = {nm: b.node_id_of(nm) for nm in ('pa1', 'pb1', 'pa2', 'pb2', 'rchm')}
        cpn['vdie'] = None; cpn['gnd'] = None
        cp_on = dict(en=False, ph=0, shunt=False)
        cpg = 1.0 / CP.get('ron', 0.5)
    Lpkg_i = svn.index('I(Lpkg)') if 'I(Lpkg)' in svn else None
    Z0 = np.sqrt((L_fivr / nph_fivr if case == 'B' else (L_fivr / nph_fivr * L_PKG / (L_fivr / nph_fivr + L_PKG)) if case == 'D' else L_PKG) / cdie)
    sch = {'k': warmup, 'done': None, 'res': None, 'SIG': 0.0, 'nsig': 0}
    st = {'vdie': VNOM, 'voff': vref_off, 't_last_off': 0.0, 't_last_f': 0.0,
          'd_off': [D_off0] * nph_off, 'd_f': [None] * (nph_fivr if hasF else 0),
          'c_off': D_off0, 'c_f': (D_f0 if hasF else 0)}
    Tp_off = 1.0 / fsw_off
    st['po'] = [1.0] * nph_off; st['pf'] = [1.0] * nph_fivr
    if bank: st['vbank'] = vb0
    st['ILs'] = i0 if case != 'D' else 0.0
    if hasF:
        Tp_f = 1.0 / fsw_fivr; st['d_f'] = [D_f0] * nph_fivr
        st['vint'] = 0.0
        mavg = max(1, int(round(Tp_f / dt))); st['ilbuf'] = [st['ILs']] * mavg; st['ilsum'] = st['ILs'] * mavg; st['ilp'] = 0
        # running estimate of total FIVR inductor current = load current filtered by the output LC is unknown;
        # feed-forward uses the change in PREDICTED load current over the next ff_lead cycles.
    prev_t = [0.0]

    VMIN = np.full(n, 9.0); VMAX = np.full(n, -9.0); VSUM = np.zeros(n); VCNT = np.zeros(n)
    V18MIN = np.full(n, 9.0); IFV = np.zeros(n); IFV2 = np.zeros(n)
    dts_ok = lambda t: True
    def obs(t, x):
        st['vdie'] = float(x[idx['vdie']]); st['voff'] = float(x[idx[sense_off]])
        if bank: st['vbank'] = float(x[bkid])
        if CP and dts_ok(t):
            ddt = t - st.get('tcp', t); st['tcp'] = t
            V = lambda nm: (0.0 if nm == 'gnd' else st['vdie'] if nm == 'vdie' else st['vbank'] if nm == 'bank' else float(x[cpn[nm]]))
            if cp_on['en'] and ddt > 0:
                # current drawn OUT of vdie by the pump switches, and delivered INTO the bank
                if cp_on['ph'] == 1: i_vd = (st['vdie'] - V('pa1')) * cpg; i_bk = 0.0
                else:
                    i_vd = (st['vdie'] - V('pb1')) * cpg * -1 + (st['vdie'] - V('pb2')) * cpg * -1
                    i_vd = -i_vd   # current leaving vdie into pb1/pb2 (positive = drawn from vdie)
                    i_bk = (V('pa1') - st['vbank']) * cpg + (V('pa2') - st['vbank']) * cpg
                st['Ein'] = st.get('Ein', 0.0) + i_vd * st['vdie'] * ddt
                st['Estored'] = st.get('Estored', 0.0) + i_bk * st['vbank'] * ddt
            if st.get('rch_on') and ddt > 0:
                ir_ = (float(x[rsid]) - st['vbank']) / (CP.get('rch', (1.8 - vb0) / (bank.get('ir', 25.0) * 1e-3)) + 0.5)
                st['Erch'] = st.get('Erch', 0.0) + ir_ * float(x[rsid]) * ddt
        if case == 'A' and Lpkg_i is not None: st['ILs'] = float(x[Lpkg_i])
        kk = int(t / cyc_s + 1e-9)
        if kk < n:
            v = st['vdie']
            if v < VMIN[kk]: VMIN[kk] = v
            if v > VMAX[kk]: VMAX[kk] = v
            VSUM[kk] += v; VCNT[kk] += 1
            if hasF: IFV[kk] += st['ILs']; IFV2[kk] += st['ILs'] ** 2
            if v18id is not None:
                w = float(x[v18id])
                if w < V18MIN[kk]: V18MIN[kk] = w
        dts = t - prev_t[0]; prev_t[0] = t
        if dts <= 0: return
        if Oidx:
            ilo = sum(float(x[j]) for j in Oidx)
            if not st_o['init']: st_o['buf'] = [ilo] * mo; st_o['sum'] = ilo * mo; st_o['init'] = True
            bb = st_o['buf']; st_o['sum'] += ilo - bb[st_o['p']]; bb[st_o['p']] = ilo; st_o['p'] = (st_o['p'] + 1) % mo
            if case == 'D':
                # current SPLIT loop: off-chip slowly takes over whatever average current the FIVR supplies
                st_split['iref'] += 2 * np.pi * 20e3 * st['ILs'] * dts; iref_o = st_split['iref']
            else:
                iref_o = ctl_off.step(vref_off - st['voff'], dts)
            st['c_off'] = min(0.95, max(0.0, float(x[vbid]) / 12.0 + Rk_off * (iref_o - st_o['sum'] / mo)))
        if hasF:
            if ctl == 'pid':
                # high-bandwidth compensator (designed: UGF ~110 MHz, PM 54 deg, GM 17 dB on the averaged
                # plant) -- a STRONG baseline, above Intel FIVR's measured 78 MHz / 40 deg (Burton 2014)
                ilr = sum(float(x[j]) for j in Lidx)
                buf = st['ilbuf']; st['ilsum'] += ilr - buf[st['ilp']]; buf[st['ilp']] = ilr
                st['ilp'] = (st['ilp'] + 1) % len(buf); st['ILs'] = st['ilsum'] / len(buf)
                e = VNOM - st['vdie']; st['vint'] += e * dts
                de = (e - st.get('eprev', e)) / dts; st['eprev'] = e
                wf = 2 * np.pi * 400e6; a = min(1.0, wf * dts)
                st['dfilt'] = st.get('dfilt', 0.0) + a * (de - st.get('dfilt', 0.0))
                c = Dnom_f * 1.8 / max(float(x[v18id]), 0.5) + 0.040 * e + 3.14e5 * st['vint'] + 1.44e-10 * st['dfilt']
                st['c_f'] = min(0.98, max(0.0, c)); return
            if ctl == 't3':
                ilr = sum(float(x[j]) for j in Lidx)
                buf = st['ilbuf']; st['ilsum'] += ilr - buf[st['ilp']]; buf[st['ilp']] = ilr
                st['ilp'] = (st['ilp'] + 1) % len(buf); st['ILs'] = st['ilsum'] / len(buf)
                dl = t3s['dl']; e = dl[t3s['dp']]; dl[t3s['dp']] = VNOM - st['vdie']; t3s['dp'] = (t3s['dp'] + 1) % len(dl)
                xi_new = t3s['xi'] + t3s['wI'] * e * dts
                u = xi_new
                for j in range(2):          # (1+s/wz)/(1+s/wp) = r + (1-r)/(1+s/wp), r = wp/wz
                    wz, wp = t3s['wz'][j], t3s['wp'][j]; r = wp / wz
                    t3s['x'][j] += (1.0 - np.exp(-wp * dts)) * (u - t3s['x'][j])
                    u = r * u + (1.0 - r) * t3s['x'][j]
                c = Dnom_f * 1.8 / max(float(x[v18id]), 0.5) + u
                if (c > 0.98 and e > 0) or (c < 0.0 and e < 0): pass     # anti-windup: hold integrator
                else: t3s['xi'] = xi_new
                st['c_f'] = min(0.98, max(0.0, c)); return
            if ctl == 'vm':
                ilr = sum(float(x[j]) for j in Lidx)
                buf = st['ilbuf']; st['ilsum'] += ilr - buf[st['ilp']]; buf[st['ilp']] = ilr
                st['ilp'] = (st['ilp'] + 1) % len(buf); st['ILs'] = st['ilsum'] / len(buf)
                # VOLTAGE-MODE: nominal duty feedforward keeps the buck's natural LC stiffness;
                # small Kp (no peaking at the output-LC resonance) + slow integrator for DC.
                e = VNOM - st['vdie']; st['vint'] += e * dts
                c = Dnom_f * 1.8 / max(float(x[v18id]), 0.5) + kp_vm * e + Ki_vm * st['vint']
                if ff_lead > 0:
                    k = int(t / cyc_s + 1e-9)
                    kf = min(n - 1, k + ff_lead); kn = min(n - 1, k)
                    # extra volt-seconds to slew L_eff by the predicted load change over the lead window
                    dI = (FF[kn] if FF is not None and kn < len(FF) else P_A[kf] - P_A[kn])
                    c += ff_gain * L_eff * dI / (1.8 * max(ff_lead * cyc_s, Tp_f))   # each phase holds a latched duty for Tp_f
                st['c_f'] = min(0.98, max(0.0, c)); return
            iref = ctl_f.step(VNOM - st['vdie'], dts) if ctl == 'full' else i0
            if ff_lead > 0:
                # scheduler->FIVR feed-forward: the issue stage's predicted load current ff_lead
                # cycles ahead (power tokens). Only the predicted CHANGE is added, so the integrator
                # still owns DC and an imperfect prediction can't cause a static error.
                k = int(t / cyc_s + 1e-9)
                kf = min(n - 1, k + ff_lead); kn = min(n - 1, k)
                iref += ff_gain * (P_A[kf] - P_A[kn])
            # AVERAGE-current sensing: one-switching-period moving average of total inductor
            # current (instantaneous sampling at D~0.5 = subharmonic instability, found by bisection)
            ilr = sum(float(x[j]) for j in Lidx)
            buf = st['ilbuf']; st['ilsum'] += ilr - buf[st['ilp']]; buf[st['ilp']] = ilr
            st['ilp'] = (st['ilp'] + 1) % len(buf); IL = st['ilsum'] / len(buf); st['ILs'] = IL
            v18 = float(x[v18id])
            c = st['vdie'] / max(v18, 0.5) + Rk * (iref - IL)
            if ctl == 'open': c = D_f0
            st['c_f'] = min(0.98, max(0.0, c))

    def swfn(t):
        m = p.SwitchStateMask(nsw)
        for k in (range(nph_off) if not (case == 'B' and ideal_in) else ()):
            ph = (t / Tp_off + k / nph_off) % 1.0
            if ph < st['po'][k]: st['d_off'][k] = st['c_off']      # latch duty at carrier wrap
            st['po'][k] = ph
            on = ph < st['d_off'][k]
            m.set(sw_index[f'oHS{k}'], on); m.set(sw_index[f'oLS{k}'], not on)
        if hasF:
            for k in range(nph_fivr):
                ph = (t / Tp_f + k / nph_fivr) % 1.0
                if pwm == 'tri':
                    # natural-sampled double-edge PWM: continuous compare against a 0..1..0 triangle
                    on = st['c_f'] > 1.0 - abs(2.0 * ph - 1.0)
                else:
                    if ph < st['pf'][k]: st['d_f'][k] = st['c_f']
                    st['pf'][k] = ph
                    on = ph < st['d_f'][k]
                m.set(sw_index[f'fHS{k}'], on); m.set(sw_index[f'fLS{k}'], not on)
        if CP:
            full = st['vbank'] >= vb0 + CP.get('full', 0.1)
            over = st['vdie'] > VNOM + bank.get('th_h', 0.01)
            if not cp_on['en'] and over and not full:
                cp_on['en'] = True; cp_on['t0'] = t
            elif cp_on['en'] and (full or not over):
                # finish the current clock period so the flying caps end a full ph1+ph2 cycle
                if ((t - cp_on['t0']) * CP.get('fclk', 1e9)) % 1.0 < 0.02 or full: cp_on['en'] = False
            ph = 1 if ((t - cp_on.get('t0', 0.0)) * CP.get('fclk', 1e9)) % 1.0 < 0.5 else 2
            cp_on['ph'] = ph if cp_on['en'] else 0
            for nm in ('P1a', 'P1b', 'P1c'): m.set(sw_index[nm], cp_on['ph'] == 1)
            for nm in ('P2a', 'P2b', 'P2c', 'P2d'): m.set(sw_index[nm], cp_on['ph'] == 2)
            st['rch_on'] = st['vbank'] < vb0
            m.set(sw_index['SWrch'], st['rch_on'])
            cp_on['shunt'] = CP.get('mode') == 'hybrid' and (full or st['vdie'] > VNOM + CP.get('th2', 0.04))
        if bank:
            m.set(sw_index['SWrel'], (st['vdie'] < VNOM - bank.get('th', 0.02)) and (st['vbank'] > st['vdie']))
        return m

    vid = idx['vdie']
    def bex(t):
        k = int(t / cyc_s + 1e-9)
        if k >= n: k = n - 1
        if sched is not None:
            while sch['k'] < k:
                sch['k'] += 1
                if sch['done'] is None and sch['k'] >= warmup:
                    tok = 0.0
                    if bank and bank.get('tok_model'):
                        tm = bank['tok_model']; kk = sch['k'] - 1
                        st.setdefault('tokm', tm['Qb'])
                        lo = I_A[max(0, kk - 5):kk].min() * 1e3 if kk > 0 else I_A[0] * 1e3
                        draw = max(0.0, I_A[kk] * 1e3 - (lo + 100.0))
                        st['tokm'] = min(tm['Qb'], st['tokm'] + tm['Ir']) - min(draw, st['tokm'])
                        tok = st['tokm']
                    elif bank:
                        floor = min(st['vdie'], VNOM - bank.get('th', 0.02))
                        tok = max(0.0, bank['cb'] * 1e-9 * (st['vbank'] - floor)) * 1e12
                    # supply-state signal (mA): surplus supply current + voltage deviation / Z0 (resistive damper target)
                    iterm = (st['ILs'] - I_A[sch['k'] - 1]) if sig_mode in ('both', 'i') else 0.0
                    vterm = (st['vdie'] - VNOM) / Z0 if sig_mode in ('both', 'v') else 0.0
                    SIG = sig_gain * (iterm + vterm) * 1e3
                    try: I_A[sch['k']] = float(sched.send((tok, SIG))) * 1e-3
                    except StopIteration as e: sch['done'] = sch['k']; sch['res'] = e.value; I_A[sch['k']:] = 27.6e-3
        out = [0.0] * nunk; out[vid] = I_A[k]
        if CP:
            if cp_on['shunt']:
                ex = st['vdie'] - VNOM - bank.get('th_h', 0.01)
                if ex > 0:
                    sink = min(0.3, ex / Z0); out[vid] += sink
                    st['Ediss'] = st.get('Ediss', 0.0) + sink * st['vdie'] * dt
            return out
        hyb_full = bank and bank.get('hybrid') and st['vbank'] >= vb0 + 0.1
        if bank and (bank.get('sink_dissipate') or hyb_full):
            ex = st['vdie'] - VNOM - bank.get('th_h', 0.01)
            if ex > 0:
                sink = min(0.3, ex / Z0); out[vid] += sink
                st['Ediss'] = st.get('Ediss', 0.0) + sink * st['vdie'] * dt
        if bank and (bank.get('harvest_v') or bank.get('hybrid')) and st['vbank'] < vb0 + 0.1:
            # BIDIRECTIONAL bank: absorb SURPLUS (overshoot above Vnom+th_h) into the reservoir instead of
            # burning it in dummy ops; boost transfer at eta, energy-correct
            ex = st['vdie'] - VNOM - bank.get('th_h', 0.01)
            if ex > 0:
                sink = min(0.3, ex / Z0)
                g = bank.get('eta', 0.85) * st['vdie'] / st['vbank']
                out[vid] += sink; out[bkid] -= g * sink; st.setdefault('hq', 0.0); st['hq'] += g * sink * dt
                st['Ein'] = st.get('Ein', 0.0) + sink * st['vdie'] * dt; st['Estored'] = st.get('Estored', 0.0) + g * sink * st['vbank'] * dt
        if bank and st['vbank'] < vb0:
            ir = bank.get('ir', 25.0) * 1e-3
            out[bkid] -= ir; out[rsid] += ir
            if bank.get('harvest') and k > 5:
                sink = min(0.15, max(0.0, I_A[k-5:k].max() - bank.get('thf', 150.0) * 1e-3 - I_A[k]))
                if sink > 0:
                    g = bank.get('eta', 0.85) * st['vdie'] / st['vbank']
                    out[vid] += sink; out[bkid] -= g * sink
        return out

    kw_sim = {}
    if sched is not None: kw_sim['should_continue'] = lambda: sch['done'] is None or sch['k'] < sch['done'] + 300
    res = p.simulate(b, t_end=n * cyc_s, dt=dt, b_extra_fn=bex, switch_fn=swfn, step_observer=obs,
                     start_from_dc_op=False, store_every=100000, **kw_sim)
    w = warmup
    if sched is not None: n = min(n, sch['k'])
    out = dict(vmin=VMIN[w:n-1], vmax=VMAX[w:n-1], vavg=VSUM[w:n-1] / np.maximum(VCNT[w:n-1], 1))
    if hasF:
        out['v18min'] = V18MIN[w:n-1]; cnt = np.maximum(VCNT[w:n-1], 1)
        out['ifivr'] = IFV[w:n-1] / cnt; out['ifivr_rms'] = float(np.sqrt((IFV2[w:n-1].sum()) / cnt.sum()))
    if sched is not None: out['sched'] = sch['res']
    out['E'] = {k: st.get(k, 0.0) for k in ('Ediss', 'Ein', 'Estored', 'Erch')}
    if ctl == 't3' and hasF: out['t3'] = {k: float(v) for k, v in t3s['design'].items()}
    return out
    S = np.asarray(res.states)
    sub = int(round(cyc_s / dt))
    vdie = S[:, vid]
    # per-cycle minimum (droop) and maximum (overshoot) so sub-cycle ripple troughs are not missed
    m = (len(vdie) - 1) // sub
    blk = vdie[1:m * sub + 1].reshape(m, sub)
    out = dict(vmin=blk.min(1), vmax=blk.max(1), vavg=blk.mean(1))
    if case == 'B':
        v18 = S[:, b.node_id_of('v18')][1:m * sub + 1].reshape(m, sub)
        out['v18min'] = v18.min(1)
    return out
