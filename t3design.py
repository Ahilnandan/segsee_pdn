"""Type-3 voltage-mode compensator for the faithful FIVR plant (case 'C').

Plant (averaged, per-core share): 16-phase buck, L_ph per phase -> L_eff = L_ph/16, die C,
series R_eff = (ron+dcr)/16 + R_grid, current-source load.  Gvd(s) = Vin / (L C s^2 + R C s + 1).
Natural-sampled double-edge (triangle 0..1) PWM: modulator gain 1 per unit duty; a lumped
loop delay Td (sensing + comparator + driver) is included.
Gc(s) = wI/s * (1+s/wz1)(1+s/wz2) / ((1+s/wp1)(1+s/wp2)); wI is solved for the target UGF.
Target: Intel FIVR (Burton APEC'14) ~78 MHz UGF, ~40 deg PM.
"""
import numpy as np

VIN = 1.8
def plant(L_ph=112e-9, nph=16, C=5e-9, ron=0.2, dcr=0.2, rgrid=20e-3):
    Le = L_ph / nph; Re = (ron + dcr) / nph + rgrid
    return Le, C, Re

def comp(fz1, fz2, fp1, fp2):
    wz1, wz2, wp1, wp2 = (2 * np.pi * f for f in (fz1, fz2, fp1, fp2))
    return lambda s: (1 + s / wz1) * (1 + s / wz2) / (s * (1 + s / wp1) * (1 + s / wp2))

def design(fc=78e6, fz=(14e6, 14e6), fp=(400e6, 600e6), Td=0.3e-9, **pk):
    Le, C, Re = plant(**pk)
    G = lambda s: VIN / (Le * C * s * s + Re * C * s + 1) * np.exp(-s * Td)
    H = comp(*fz, *fp)
    sc = 2j * np.pi * fc
    wI = 1.0 / abs(G(sc) * H(sc))
    f = np.logspace(5, 10, 200000); s = 2j * np.pi * f
    T = wI * G(s) * H(s)
    mag = np.abs(T); ph = np.unwrap(np.angle(T)) * 180 / np.pi
    i_c = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    ugf = f[i_c[-1]]; pm = 180 + ph[i_c[-1]]
    while pm > 180: pm -= 360
    i_g = np.where((ph[:-1] > -180) & (ph[1:] <= -180))[0]
    gm = -20 * np.log10(mag[i_g[0]]) if len(i_g) else np.inf
    # closed-loop output impedance (A/V): Zout = Zol / (1 + T)
    Zol = (Re + s * Le) / (Le * C * s * s + Re * C * s + 1)
    Zcl = np.abs(Zol / (1 + T))
    return dict(wI=wI, ugf=ugf, pm=pm, gm=gm, n_cross=len(i_c), Zpk=Zcl.max(), fZpk=f[Zcl.argmax()],
                Zol_pk=np.abs(Zol).max(), f0=1 / (2 * np.pi * np.sqrt(Le * C)), Le=Le, Re=Re)

if __name__ == '__main__':
    for fz in (8e6, 10e6, 14e6, 20e6):
        for fp in ((300e6, 600e6), (400e6, 800e6), (250e6, 1e9)):
            for Td in (0.3e-9, 0.6e-9):
                d = design(fz=(fz, fz), fp=fp, Td=Td)
                print(f"fz {fz/1e6:4.0f} fp {fp[0]/1e6:4.0f}/{fp[1]/1e6:4.0f} Td {Td*1e9:.1f} | UGF {d['ugf']/1e6:5.1f} PM {d['pm']:5.1f} GM {d['gm']:5.1f} "
                      f"xings {d['n_cross']} | Zcl pk {d['Zpk']*1e3:6.1f} mOhm @ {d['fZpk']/1e6:5.1f} MHz (open {d['Zol_pk']*1e3:6.0f}) f0 {d['f0']/1e6:.1f}")

# ---- PARALLEL (direct + FIVR) plant: vdie fed by the FIVR (L_eff, R_eff) AND the package path (Lpkg, Rpkg)
# from the off-chip bulk node, which is AC-stiff (mF-class bulk, per-core share) -> treat as AC ground.
def design_par(fc=78e6, fz=(10e6, 10e6), fp=(400e6, 800e6), Td=0.6e-9, Lpkg=507e-12, Rpkg=79.6e-3, **pk):
    Le, C, Re = plant(**pk)
    def G(s):
        Zp = 1.0 / (s * C + 1.0 / (s * Lpkg + Rpkg))          # die C || package branch
        return VIN * Zp / (s * Le + Re + Zp) * np.exp(-s * Td)
    H = comp(*fz, *fp); sc = 2j * np.pi * fc
    wI = 1.0 / abs(G(sc) * H(sc))
    f = np.logspace(5, 10, 200000); s = 2j * np.pi * f
    T = wI * G(s) * H(s); mag = np.abs(T); ph = np.unwrap(np.angle(T)) * 180 / np.pi
    ic = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    pms = [((180 + ph[i]) + 180) % 360 - 180 for i in ic]
    ig = np.where((ph[:-1] > -180) & (ph[1:] <= -180))[0]
    gm = -20 * np.log10(mag[ig[0]]) if len(ig) else np.inf
    Zol = 1.0 / (s * C + 1.0 / (s * Lpkg + Rpkg) + 1.0 / (s * Le + Re))
    Zcl = np.abs(Zol / (1 + T))
    Zc_only = np.abs(1.0 / (s * C + 1.0 / (s * Lpkg + Rpkg)))
    return dict(wI=wI, xings=[float(f[i]) for i in ic], pms=pms, gm=gm, Zpk=Zcl.max(), fZpk=f[Zcl.argmax()],
                Zol_pk=np.abs(Zol).max(), Zdirect_pk=Zc_only.max(), fdirect=f[Zc_only.argmax()])

if __name__ == '__main__' and 'par' in __import__('sys').argv:
    for fc in (78e6, 150e6, 250e6, 400e6):
        for fz in (10e6, 30e6, 60e6):
            for Td in (0.3e-9, 0.6e-9):
                try: d = design_par(fc=fc, fz=(fz, fz), Td=Td)
                except Exception as e: print(fc, fz, 'fail', e); continue
                print(f"fc {fc/1e6:4.0f} fz {fz/1e6:3.0f} Td {Td*1e9:.1f} | crossings {[round(x/1e6,1) for x in d['xings']]} PM {[round(p,1) for p in d['pms']]} GM {d['gm']:5.1f} | "
                      f"Zcl pk {d['Zpk']*1e3:6.1f} mOhm @ {d['fZpk']/1e6:5.1f} | direct-only pk {d['Zdirect_pk']*1e3:6.1f} @ {d['fdirect']/1e6:.0f} MHz")

def design_par_gain(wI, fz=(3e6, 3e6), fp=(600e6, 1.2e9), Td=0.6e-9, Lpkg=507e-12, Rpkg=79.6e-3, **pk):
    """parallel plant, compensator gain given directly (the target-UGF parametrisation is ill-posed here:
    the package path makes |G| tiny below ~100 MHz and resonant near 150 MHz)."""
    Le, C, Re = plant(**pk)
    f = np.logspace(5, 10, 200000); s = 2j * np.pi * f
    Zp = 1.0 / (s * C + 1.0 / (s * Lpkg + Rpkg))
    G = VIN * Zp / (s * Le + Re + Zp) * np.exp(-s * Td)
    T = wI * G * comp(*fz, *fp)(s); mag = np.abs(T); ph = np.unwrap(np.angle(T)) * 180 / np.pi
    ic = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    pms = [float(((180 + ph[i]) + 180) % 360 - 180) for i in ic]
    ig = np.where((ph[:-1] > -180) & (ph[1:] <= -180))[0]
    gm = float(-20 * np.log10(mag[ig[0]])) if len(ig) else np.inf
    Zol = 1.0 / (s * C + 1.0 / (s * Lpkg + Rpkg) + 1.0 / (s * Le + Re))
    Zcl = np.abs(Zol / (1 + T))
    return dict(wI=wI, fcmax=float(f[ic[-1]]) if len(ic) else 0.0, pmin=min(pms) if pms else 180.0, gm=gm,
                Zpk=float(Zcl.max()), fZpk=float(f[Zcl.argmax()]))

def best_par(fcap, **kw):
    """largest-benefit stable design whose highest crossover <= fcap, PM>=40, GM>=6 dB"""
    best = None
    for wI in np.logspace(4, 8, 120):
        for fz in (1e6, 3e6, 10e6):
            d = design_par_gain(wI, fz=(fz, fz), **kw)
            if d['fcmax'] <= fcap and d['pmin'] >= 40 and d['gm'] >= 6:
                if best is None or d['Zpk'] < best['Zpk']: best = dict(d, fz=fz)
    return best
