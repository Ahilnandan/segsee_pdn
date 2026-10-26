"""Decap discharge/recharge on a load step (0.5 nH / 50 mOhm package, 8.6 nF die decap with 35 mOhm ESR):
load +150 mA at 20 ns, back at 220 ns. Shows decap voltage/charge, decap current, package current."""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
L, R, C, ESR, V0, I0 = 0.5e-9, 50e-3, 8.6e-9, 35e-3, 1.1, 0.10
dt = 0.01e-9; t = np.arange(0, 400e-9, dt)
il = np.where((t > 20e-9) & (t < 220e-9), I0 + 0.15, I0)
iL = I0; vC = V0 - 0.0; vs = V0 + I0 * R          # supply sized so the die sits at 1.1 V at the starting load
VD, ICAP, IPKG, Q = [], [], [], []
for k in range(len(t)):
    vd = vC + ESR * (iL - il[k]); icap = iL - il[k]       # current INTO the decap (negative = decap discharging)
    iL += (vs - vd - R * iL) / L * dt; vC += icap / C * dt
    VD.append(vd); ICAP.append(icap); IPKG.append(iL); Q.append(C * vC)
VD, ICAP, IPKG, Q = map(np.array, (VD, ICAP, IPKG, Q)); tn = t * 1e9
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK, 'axes.facecolor': SURF, 'figure.facecolor': SURF})
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ax[0].plot(tn, il * 1e3, color=INK2, lw=1.4, label='core load'); ax[0].plot(tn, IPKG * 1e3, color='#2a78d6', lw=2, label='package current (from supply)')
ax[0].set_ylabel('mA'); ax[0].legend(frameon=False, fontsize=8, loc='upper right'); ax[0].set_title('(a) the supply current ramps up slowly, overshoots, rings, settles', loc='left', fontsize=9.5)
ax[1].plot(tn, ICAP * 1e3, color='#eb6834', lw=2); ax[1].axhline(0, color=INK2, lw=0.8)
ax[1].text(25, -120, 'decap DISCHARGING\n(gives charge)', fontsize=8, color='#eb6834'); ax[1].text(30, 60, 'decap RECHARGING', fontsize=8, color='#eb6834')
ax[1].set_ylabel('mA into decap'); ax[1].set_title('(b) decap current: out during the dip, back in afterwards (net ≈ 0)', loc='left', fontsize=9.5)
ax[2].plot(tn, (VD - V0) * 1e3, color='#1baf7a', lw=2); ax[2].set_ylabel('die voltage − 1.1 V (mV)')
ax2 = ax[2].twinx(); ax2.set_ylim(*(np.array(ax[2].get_ylim()) * 1e-3 * C * 1e9 + C * V0 * 1e9)); ax2.set_ylabel('decap charge (nC)', color=INK2)
ax[2].set_xlabel('time (ns)'); ax[2].set_title('(c) decap voltage/charge: dips, refills, settles 7.5 mV lower (extra I×R of the package)', loc='left', fontsize=9.5)
for a in ax: a.grid(True, color=GRID, lw=0.5); [a.spines[s].set_visible(False) for s in ('top', 'right')]
fig.savefig('/home/ahilnandan/Documents/PMIC/pdn_work/docs/figs/decap_charge.png', dpi=140, bbox_inches='tight')
k1 = np.argmin(VD[:20000]); print('stored charge %.2f nC; min dip %.1f mV at %.1f ns; charge given %.3f nC (%.1f%%); net charge into decap 20-200ns %.4f nC'
      % (C * V0 * 1e9, (VD[k1] - V0) * 1e3, tn[k1], (C * V0 - Q[k1]) * 1e9, 100 * (C * V0 - Q[k1]) / (C * V0), ICAP[2000:20000].sum() * dt * 1e9))
