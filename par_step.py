"""Validate case D (direct + FIVR in parallel) vs A (direct only) and C (FIVR only) on a bursty step load."""
import numpy as np, time
from pdnchain import simulate_chain
C = dict(nph_fivr=16, L_fivr=112e-9, ron_f=0.2, dcr_f=0.2, pwm='tri', ctl='t3', fsw_fivr=140e6)
n = 30000; I = np.full(n, 100.0)
for k0 in range(2000, 28000, 3000): I[k0:k0 + 600] = 400.0     # 300 mA bursts
I[20000:26000] = 250.0                                           # long DC shift: does the split loop hand it to the off-chip?
for name, kw in [('A direct only', dict(case='A')), ('C FIVR only', dict(case='B', ideal_in=True, **C)), ('D direct+FIVR', dict(case='D', **C))]:
    t0 = time.time(); o = simulate_chain(I, warmup=3000, **kw)
    v = 1e3 * (0.9 - o['vmin'])
    print(f"{name:15s} droop max {v.max():6.1f} mV p99.9 {np.percentile(v,99.9):6.1f} | over {1e3*(o['vmax'].max()-0.9):6.1f} | "
          f"DC err {1e3*(o['vavg'][1000:1900].mean()-0.9):+.2f} / {1e3*(o['vavg'][25000:25900].mean()-0.9):+.2f} mV | {time.time()-t0:.0f}s")
    if 't3' in o: print('   loop:', {k: (round(v, 3) if abs(v) < 1e3 else f'{v:.3g}') for k, v in o['t3'].items()})
