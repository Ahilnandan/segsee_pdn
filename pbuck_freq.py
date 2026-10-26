"""Analytical frequency-domain model: one BOOM core fed by TWO supplies in parallel, real (unscaled) units.

  Supply 1 (direct path):  board VR1 3.3 V -> 0.9 V  ->  board trace  ->  package decap  ->  package-to-die  ->  DIE
  Supply 2 (mini-buck):    board VR2 1.8 V (stiff)   ->  on-die switches -> small air-core inductor(s) -> DIE

Everything is linear and small-signal except where stated (ripple is a large-signal periodic waveform; the
on-demand injector is analysed in the time domain). Each section implements one numbered equation of
docs/parallel-buck-derivation.md and prints the numbers used there; figures go to docs/figs/.
"""
import json, numpy as np, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

FIG = os.path.join(os.path.dirname(__file__), 'docs', 'figs'); os.makedirs(FIG, exist_ok=True)
f = np.logspace(2, 9.7, 6000); s = 2j * np.pi * f; w = 2 * np.pi * f

# ------------------------------------------------------------------ parameters (real units)
P = dict(
    VDD=0.9,
    # load (one BOOM core, current model +-30-40%): from replay traces
    I_avg=0.20, I_peak=0.43,
    # die
    C_die=5e-9, R_die=10e-3,
    # package-to-die (project's die-side model, kept for comparability -- provenance unverified)
    L_p=507e-12, R_p=79.6e-3,
    # package decap
    C_pk=1e-6, ESL_pk=100e-12, ESR_pk=20e-3,
    # board trace VR1 -> package
    L_b=2e-9, R_b=3e-3,
    # VR1 (board buck, 3.3 V -> 0.9 V, sized for ~1 A)
    Vin1=3.3, f1=2e6, L1=2.2e-6, DCR1=40e-3, C1=47e-6, ESR1=3e-3, ESL1=0.5e-9, fc1=100e3,
    # supply 2 / mini-buck
    Vin2=1.8, L_mb=2e-9, DCR_mb=20e-3, Ron_mb=0.10, Td=0.6e-9,
    # switch energy per phase per switching cycle at Ron=0.1 ohm (gate drive of a ~10 mm-wide 22nm switch pair, see doc)
    Esw_ref=65e-12, Ron_ref=0.10,
)
def par(*Z):
    return 1.0 / sum(1.0 / z for z in Z)

# ------------------------------------------------------------------ (1) VR1 closed-loop output impedance
def Z_vr1(p=P):
    Zl = s * p['L1'] + p['DCR1']                         # inductor branch (VR1's own switch node is a voltage source)
    Zc = 1.0 / (s * p['C1']) + p['ESR1'] + s * p['ESL1']  # output capacitor
    Zol = par(Zl, Zc)
    T = (2 * np.pi * p['fc1']) / s                       # loop gain approximated as integrator crossing at fc1 (type-III shaped)
    return Zol / (1.0 + T)
# ------------------------------------------------------------------ (2) direct path seen from the die (excluding die C)
def Z_direct(p=P):
    Z1 = Z_vr1(p) + s * p['L_b'] + p['R_b']
    Zpk = 1.0 / (s * p['C_pk']) + p['ESR_pk'] + s * p['ESL_pk']
    Z2 = par(Z1, Zpk)
    return Z2 + s * p['L_p'] + p['R_p']
Zdie = lambda p=P: 1.0 / (s * p['C_die']) + p['R_die']
# ------------------------------------------------------------------ (3) mini-buck branch admittance (continuous mode, small signal)
def Gc_type3(fz, fp1, fp2, wI):
    return wI / s * (1 + s / (2 * np.pi * fz)) ** 2 / ((1 + s / (2 * np.pi * fp1)) * (1 + s / (2 * np.pi * fp2)))
def Y_buck(N, L, Rph, Gc=None, p=P):
    """N interleaved phases sharing one duty command. i_b = N (Vin*d - v)/(sL+R), d = -Gc*e^{-sTd}*v (+ fixed duty)
       => Y_b = N (1 + Vin*Gc*e^{-sTd}) / (sL + R).  Gc=None -> open loop (fixed duty): a plain inductive path."""
    loop = 0 if Gc is None else p['Vin2'] * Gc * np.exp(-s * p['Td'])
    return N * (1 + loop) / (s * L + Rph)
def Z_total(Yb=None, p=P):
    Y = 1 / Zdie(p) + 1 / Z_direct(p)
    if Yb is not None: Y = Y + Yb
    return 1 / Y
def loop_gain(N, L, Rph, Gc, p=P):
    """loop seen by the mini-buck controller: T = Vin*Gc*e^{-sTd} * N/(sL+R) * Z_shunt, Z_shunt = die C || direct path"""
    Zsh = 1 / (1 / Zdie(p) + 1 / Z_direct(p))
    return p['Vin2'] * Gc * np.exp(-s * p['Td']) * N / (s * L + Rph) * Zsh
def margins(T):
    mag = np.abs(T); ph = np.unwrap(np.angle(T)) * 180 / np.pi
    ic = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    pm = [float(((180 + ph[i]) + 180) % 360 - 180) for i in ic]
    ig = np.where((ph[:-1] > -180) & (ph[1:] <= -180))[0]
    gm = float(-20 * np.log10(mag[ig[0]])) if len(ig) else float('inf')
    return [float(f[i]) for i in ic], pm, gm
def design_loop(N, L, Rph, fcap, p=P):
    """largest gain (lowest peak |Z|) type-3 whose crossings all lie below fcap with PM>=45 deg, GM>=6 dB"""
    best = None
    for fz in (2e6, 5e6, 10e6, 20e6):
        for fpm in (2.0, 4.0):
            for wI in np.logspace(4, 10, 160):
                Gc = Gc_type3(fz, fpm * fcap, 2 * fpm * fcap, wI)
                xs, pm, gm = margins(loop_gain(N, L, Rph, Gc, p))
                if not xs or max(xs) > fcap or min(pm) < 45 or gm < 6: continue
                Zp = np.abs(Z_total(Y_buck(N, L, Rph, Gc, p), p)).max()
                if best is None or Zp < best['Zpk']:
                    best = dict(fz=fz, fp=fpm * fcap, wI=float(wI), xs=xs, pm=pm, gm=gm, Zpk=float(Zp), Gc=Gc)
    return best
# ------------------------------------------------------------------ (4) interleaved ripple
def ripple(N, L, fsw, D, p=P, kmax=60):
    """per-phase triangle ripple dI = Vin D(1-D)/(L fsw). Fourier: phase current harmonic k has amplitude
       a_k = dI * 2*|sin(pi k D)| / (pi^2 k^2 D (1-D)) (peak). Interleaving N phases: sum_j e^{-i2pi k j/N} = N if N|k else 0.
       Output ripple voltage harmonic = N a_k |Z_shunt(k fsw)| for N|k."""
    dI = p['Vin2'] * D * (1 - D) / (L * fsw)
    Zsh = lambda ff: 1 / (1 / (1 / (2j * np.pi * ff * p['C_die']) + p['R_die']) + 1 / Z_direct_at(ff, p))
    harm = []
    for k in range(1, kmax + 1):
        if k % N: continue
        ak = dI * 2 * abs(np.sin(np.pi * k * D)) / (np.pi ** 2 * k ** 2 * D * (1 - D))
        harm.append((k * fsw, N * ak, N * ak * abs(Zsh(k * fsw))))
    vpp = 2 * sum(h[2] for h in harm)                     # conservative (all harmonics in phase)
    return dI, harm, vpp
def Z_direct_at(ff, p=P):
    global s
    s0 = s; s = np.array([2j * np.pi * ff]); z = Z_direct(p)[0]; s = s0; return z
# ------------------------------------------------------------------ (5) losses, continuous mode
def losses_cont(N, L, fsw, Idc, p=P, Ron=None):
    Ron = p['Ron_mb'] if Ron is None else Ron
    D = p['VDD'] / p['Vin2']; dI = p['Vin2'] * D * (1 - D) / (L * fsw)
    Irms2 = (Idc / N) ** 2 + dI ** 2 / 12
    Pc = N * Irms2 * (Ron + p['DCR_mb'])
    Esw = p['Esw_ref'] * p['Ron_ref'] / Ron                 # gate charge scales with switch width ~ 1/Ron
    Psw = N * fsw * Esw
    return dict(dI=dI, Pcond=Pc, Psw=Psw, Ptot=Pc + Psw)
# ------------------------------------------------------------------ (6) on-demand injector (time domain, closed form)
def injector(L, dI_load, p=P, Ron=None, cyc=1e-9):
    Ron = p['Ron_mb'] if Ron is None else Ron
    slew = (p['Vin2'] - p['VDD']) / L                      # A/s while high side on
    t_up = dI_load / slew; t_dn = dI_load * L / p['VDD']   # freewheel down through low side
    # charge delivered in one event (triangle-ish up, hold, down) and its conduction + 2 transitions of switching
    return dict(slew_A_per_ns=slew * 1e-9, t_up_ns=t_up * 1e9, t_down_ns=t_dn * 1e9,
                E_sw_event_pJ=2 * p['Esw_ref'] * p['Ron_ref'] / Ron * 1e12 / 2)

if __name__ == '__main__':
    R = {}
    Zd = Z_direct(); Zn = Z_total()
    # ---- direct-path landmarks
    ipk = np.argmax(np.abs(Zn)); R['direct_only'] = dict(Zpk_mOhm=float(np.abs(Zn).max() * 1e3), f_pk_MHz=float(f[ipk] / 1e6),
        Z_1MHz_mOhm=float(np.abs(Zn[np.argmin(abs(f - 1e6))]) * 1e3), Z_100kHz_mOhm=float(np.abs(Zn[np.argmin(abs(f - 1e5))]) * 1e3))
    # ---- mini-buck variants
    cases = {}
    for N in (1, 2, 4):
        for fsw in (140e6, 500e6):
            Rph = P['Ron_mb'] + P['DCR_mb']
            Yo = Y_buck(N, P['L_mb'], Rph)                  # open loop (fixed duty): passive inductive parallel path
            Zo = Z_total(Yo)
            d = design_loop(N, P['L_mb'], Rph, fsw / 5)      # conventional crossover limit fsw/5
            Zc = Z_total(Y_buck(N, P['L_mb'], Rph, d['Gc'])) if d else None
            dI, harm, vpp = ripple(N, P['L_mb'], fsw, P['VDD'] / P['Vin2'])
            lz = losses_cont(N, P['L_mb'], fsw, Idc=0.0); la = losses_cont(N, P['L_mb'], fsw, Idc=P['I_avg'])
            key = f'N{N}_{int(fsw/1e6)}MHz'
            cases[key] = dict(N=N, fsw=fsw, Zpk_openloop_mOhm=float(np.abs(Zo).max() * 1e3),
                              Zpk_closed_mOhm=(d['Zpk'] * 1e3 if d else None), loop=(None if not d else dict(fz=d['fz'], fp=d['fp'], wI=d['wI'], xs=d['xs'], pm=d['pm'], gm=d['gm'])),
                              ripple_dI_A=dI, ripple_vpp_mV=vpp * 1e3, first_ripple_harm_MHz=harm[0][0] / 1e6 if harm else None,
                              loss_idle_mW=lz['Ptot'] * 1e3, loss_carryDC_mW=la['Ptot'] * 1e3, Psw_mW=la['Psw'] * 1e3)
            if N == 4 and fsw == 500e6: Zsave = (Zo, Zc)
    R['minibuck'] = cases
    R['injector'] = {f'L{int(L*1e9)}nH': injector(L, 0.25) for L in (1e-9, 2e-9, 5e-9)}
    R['P_load_mW'] = P['VDD'] * P['I_avg'] * 1e3
    # ---- load spectrum and linear droop prediction for dhrystone (policy A)
    try:
        It = np.load('It_dhrystone_A.npy'); I = It.sum(1)
        x = (I - I.mean()) * 1e-3; nF = len(x); F = np.fft.rfft(x); ff = np.fft.rfftfreq(nF, 1e-9)
        def droop(Zfun):
            global s
            s0 = s; s = 2j * np.pi * np.maximum(ff, 1.0); Z = Zfun(); s = s0
            v = np.fft.irfft(F * Z, nF); return -v * 1e3
        vd = droop(lambda: Z_total())
        R['dhry_linear_droop_direct_p999_mV'] = float(np.percentile(vd, 99.9))
        R['load_spectrum'] = dict(frac_power_above_10MHz=float((np.abs(F[ff > 1e7]) ** 2).sum() / (np.abs(F[1:]) ** 2).sum()),
                                  frac_power_50_200MHz=float((np.abs(F[(ff > 5e7) & (ff < 2e8)]) ** 2).sum() / (np.abs(F[1:]) ** 2).sum()))
        # figure: load spectrum
        plt.figure(figsize=(7, 3.2)); band = ff > 1e4
        plt.loglog(ff[band], np.abs(F[band]) / nF * 2 * 1e3, lw=0.4); plt.xlabel('frequency (Hz)'); plt.ylabel('|I(f)| (mA)')
        plt.title('dhrystone core current spectrum (1 ns cycles)'); plt.tight_layout(); plt.savefig(f'{FIG}/load_spectrum.png', dpi=140); plt.close()
    except FileNotFoundError:
        pass
    # ---- figures: impedance
    plt.figure(figsize=(7, 4))
    plt.loglog(f, np.abs(Zd) * 1e3, label='direct path only (no die C)', lw=1)
    plt.loglog(f, np.abs(Zn) * 1e3, label='die: direct path || die C', lw=1.6)
    plt.loglog(f, np.abs(Zsave[0]) * 1e3, '--', label='+ 4-phase mini-buck, fixed duty (2 nH)')
    if Zsave[1] is not None: plt.loglog(f, np.abs(Zsave[1]) * 1e3, label='+ 4-phase mini-buck @500 MHz, closed loop (fc<=100 MHz)')
    plt.xlabel('frequency (Hz)'); plt.ylabel('|Z| seen by the core (mOhm)'); plt.ylim(1, 1e4); plt.grid(True, which='both', alpha=0.25)
    plt.legend(fontsize=7); plt.tight_layout(); plt.savefig(f'{FIG}/impedance.png', dpi=140); plt.close()
    json.dump(R, open('pbuck_freq.json', 'w'), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float))

# ------------------------------------------------------------------ second pass: fair optimisation + droop per config
def best_cont_loss(N, L, Idc=0.0, p=P):
    """minimise continuous-mode loss over fsw and switch size (Ron): P(N,f,Ron) = N[dI^2/12 (Ron+DCR) + (Idc/N)^2(Ron+DCR) + f Eref Rref/Ron]"""
    best = None
    for fsw in np.logspace(7.5, 9.3, 160):
        for Ron in np.logspace(-2.3, 0.3, 160):
            l = losses_cont(N, L, fsw, Idc, p, Ron)
            if best is None or l['Ptot'] < best['P']: best = dict(P=l['Ptot'], fsw=fsw, Ron=Ron, dI=l['dI'], Pc=l['Pcond'], Psw=l['Psw'])
    return best
def droop_lin(Y_fn, It_path='It_dhrystone_A.npy'):
    global s
    It = np.load(It_path); I = It.sum(1); x = (I - I.mean()) * 1e-3; nF = len(x)
    F = np.fft.rfft(x); ff = np.fft.rfftfreq(nF, 1e-9)
    s0 = s; s = 2j * np.pi * np.maximum(ff, 1.0); Z = Z_total(Y_fn()); s = s0
    v = -np.fft.irfft(F * Z, nF) * 1e3
    return float(np.percentile(v, 99.9)), float(v.max())
if __name__ == '__main__':
    print('\n==== SUMMARY ====')
    print(f"load power {P['VDD']*P['I_avg']*1e3:.0f} mW (0.9 V x 0.20 A avg)")
    print(f"direct path only: peak |Z| {R['direct_only']['Zpk_mOhm']:.0f} mOhm at {R['direct_only']['f_pk_MHz']:.0f} MHz; "
          f"|Z| @1 MHz {R['direct_only']['Z_1MHz_mOhm']:.1f}, @100 kHz {R['direct_only']['Z_100kHz_mOhm']:.1f} mOhm")
    print('dhrystone linear droop p99.9 / max (mV):')
    print('  direct only                      ', droop_lin(lambda: None))
    for N in (1, 2, 4):
        Rph = P['Ron_mb'] + P['DCR_mb']
        print(f'  + {N}-ph 2 nH mini-buck fixed duty  ', droop_lin(lambda: Y_buck(N, P['L_mb'], Rph)))
        for fsw in (140e6, 500e6):
            d = design_loop(N, P['L_mb'], Rph, fsw / 5)
            if d: print(f'  + {N}-ph closed loop, fc<=fsw/5={fsw/5e6:.0f} MHz', droop_lin(lambda: Y_buck(N, P['L_mb'], Rph, Gc_type3(d['fz'], d['fp'], 2 * d['fp'], d['wI']))), f"loop xings {[round(x/1e6,1) for x in d['xs']]} MHz PM {[round(x) for x in d['pm']]} GM {d['gm']:.0f} dB")
    print('ripple vs duty (N=1,2,4; 2 nH; 160 MHz):')
    for N in (1, 2, 4):
        for D in (0.48, 0.50, 0.52):
            dI, harm, vpp = ripple(N, 2e-9, 160e6, D)
            print(f'   N{N} D{D:.2f}: phase dI {dI:.2f} A pp, die ripple <= {vpp*1e3:6.2f} mV pp (first harmonic {harm[0][0]/1e6 if harm else 0:.0f} MHz)')
    print('continuous-mode loss, optimised over fsw and switch size (idle = carrying no DC; then carrying 0.2 A):')
    for L in (1e-9, 2e-9, 5e-9):
        for N in (1, 2, 4):
            b0 = best_cont_loss(N, L, 0.0); b1 = best_cont_loss(N, L, 0.2)
            print(f'   L {L*1e9:.0f} nH N{N}: idle {b0["P"]*1e3:5.1f} mW ({100*b0["P"]/0.18:4.1f}% of core) at fsw {b0["fsw"]/1e6:5.0f} MHz Ron {b0["Ron"]:.3f} ohm dI {b0["dI"]:.2f} A | '
                  f'carrying 0.2 A: {b1["P"]*1e3:5.1f} mW')
