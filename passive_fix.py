"""Step 1 of the fix: passive detuning / damping of the die-C x package-L resonance, real units, exact linear droop on traces.
Options: (a) more die capacitance, (b) lower package inductance (more power bumps), (c) damped decap: extra C_d with series R_d
in parallel with the die (no DC path -> no static loss). Droop = -IFFT(Z_total(f) I(f)) on the per-cycle current (LTI, exact)."""
import json, numpy as np
import pbuck_freq as M
MIM = 12e-9    # assumed on-die MIM density: 12 nF per mm^2 (high-density MIM class; deep-trench would be ~10x denser)
WL = ('qsort', 'dhrystone', 'mm')
I = {wl: np.load(f'It_{wl}_A.npy').astype(float).sum(1) for wl in WL}
def droops(p, Cd=0.0, Rd=0.0):
    out = {}
    for wl in WL:
        x = (I[wl] - I[wl].mean()) * 1e-3; nF = len(x); F = np.fft.rfft(x); ff = np.fft.rfftfreq(nF, 1e-9)
        s = 2j * np.pi * np.maximum(ff, 1.0); s0 = M.s; M.s = s
        Y = 1 / M.Zdie(p) + 1 / M.Z_direct(p)
        if Cd > 0: Y = Y + 1 / (1 / (s * Cd) + Rd)
        M.s = s0
        v = -np.fft.irfft(F / Y, nF) * 1e3
        out[wl] = (round(float(np.percentile(v, 99.9)), 1), round(float(v.max()), 1))
    return out
def zpk(p, Cd=0.0, Rd=0.0):
    s = M.s; Y = 1 / M.Zdie(p) + 1 / M.Z_direct(p)
    if Cd > 0: Y = Y + 1 / (1 / (s * Cd) + Rd)
    Z = np.abs(1 / Y); i = Z.argmax(); return round(float(Z[i] * 1e3)), round(float(M.f[i] / 1e6))
R = {}
base = dict(M.P)
print(f"{'option':52s} | {'Zpk mOhm @MHz':>14s} | {'extra area':>10s} | p99.9 / max droop (mV): qsort, dhrystone, mm")
def row(name, p, Cd=0.0, Rd=0.0, area=0.0):
    z = zpk(p, Cd, Rd); d = droops(p, Cd, Rd); R[name] = dict(Zpk=z, droop=d, area_mm2=area)
    print(f"{name:52s} | {z[0]:6d} @{z[1]:4d} | {area:7.2f} mm2 | " + ', '.join(f'{d[w][0]:5.1f}/{d[w][1]:5.1f}' for w in WL))
row('baseline: 5 nF die, 507 pH package', base)
for C in (10e-9, 20e-9, 40e-9):
    row(f'(a) die C {C*1e9:.0f} nF', dict(base, C_die=C), area=(C - 5e-9) / MIM)
for L in (250e-12, 125e-12):
    row(f'(b) package L {L*1e12:.0f} pH (x{507/L/1e12:.0f} power bumps)', dict(base, L_p=L, R_p=base['R_p'] * L / 507e-12))
# (c) damped decap: classic optimum R_d ~ Z0 = sqrt(L/C_total); sweep C_d and R_d
Z0 = np.sqrt(base['L_p'] / base['C_die'])
print(f'   characteristic impedance Z0 = sqrt(L/C) = {Z0*1e3:.0f} mOhm')
best = None
for Cd in (2.5e-9, 5e-9, 10e-9):
    for Rd in (0.1, 0.2, 0.3, 0.45, 0.7):
        d = droops(base, Cd, Rd); m = np.mean([d[w][0] for w in WL])
        if best is None or m < best[0]: best = (m, Cd, Rd)
        if Rd in (0.2, 0.45): row(f'(c) damped decap C_d {Cd*1e9:.1f} nF + R_d {Rd*1e3:.0f} mOhm', base, Cd, Rd, area=Cd / MIM)
print(f'   best damped decap: C_d {best[1]*1e9:.1f} nF, R_d {best[2]*1e3:.0f} mOhm (mean p99.9 {best[0]:.1f} mV)')
row(f'(c*) best damped decap', base, best[1], best[2], area=best[1] / MIM)
# iso-area comparison: same 5 nF extra either as plain cap or as damped decap; and damped decap with the EXISTING 5 nF split
row('(a) iso-area: +5 nF plain die cap', dict(base, C_die=10e-9), area=5e-9 / MIM)
row('(c) iso-area: +5 nF damped (R_d best)', base, 5e-9, best[2], area=5e-9 / MIM)
for frac in (0.3, 0.5):
    p = dict(base, C_die=5e-9 * (1 - frac))
    bb = min(((np.mean([v[0] for v in droops(p, 5e-9 * frac, Rd).values()]), Rd) for Rd in (0.1, 0.2, 0.3, 0.45, 0.7, 1.0)))
    row(f'(c) ZERO extra area: {100*frac:.0f}% of the 5 nF as damped branch (R_d {bb[1]*1e3:.0f} mOhm)', p, 5e-9 * frac, bb[1], area=0.0)
# combination: damped + lower package L
row('(b)+(c) 250 pH package + best damped decap', dict(base, L_p=250e-12, R_p=base['R_p'] * 250 / 507), best[1], best[2], area=best[1] / MIM)
json.dump(R, open('passive_fix.json', 'w'), indent=1)
