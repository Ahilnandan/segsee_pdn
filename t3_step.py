"""Validate the faithful FIVR (case C preset): stability, DC accuracy, step droop vs averaged-model prediction."""
import sys, numpy as np, time
from pdnchain import simulate_chain
C = dict(case='B', nph_fivr=16, L_fivr=112e-9, ron_f=0.2, dcr_f=0.2, pwm='tri', ctl='t3', ideal_in=True)
n = 3000
I = np.full(n, 100.0); I[1000:2000] = 400.0
for fsw in (140, 70):
    for kw in (dict(), dict(ctl='vm', kp_vm=0.0, fint=0.5e6)):
        cfg = dict(C, fsw_fivr=fsw * 1e6, **kw)
        t0 = time.time(); o = simulate_chain(I, **cfg); el = time.time() - t0
        v = o['vavg']
        print(f"fsw {fsw} ctl {cfg['ctl']}: DC@100mA {1e3*v[900:1000].mean():.2f} mV-avg  DC@400 {1e3*v[1900:2000].mean():.2f} | "
              f"step droop {1e3*(0.9-o['vmin'][1000:1300].min()):.1f} mV  release over {1e3*(o['vmax'][2000:2300].max()-0.9):.1f} mV | "
              f"ripple(pp, steady) {1e3*(o['vmax'][800:1000].max()-o['vmin'][800:1000].min()):.2f} mV | {el:.0f}s")
        if 't3' in o: print('   design', {k: round(v, 4) if abs(v) < 1e3 else f'{v:.3g}' for k, v in o['t3'].items()})
