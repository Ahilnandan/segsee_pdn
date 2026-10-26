import json, numpy as np, par_exp
from ident import kernel
CF4 = dict(par_exp.CF, L_fivr=448e-9)
kw = dict(case='D', fsw_fivr=35e6, t3=par_exp.dpar(35e6, False, 448e-9), **CF4)
ker = json.load(open('kernels.json'))
for name, extra in [('D35L4', dict(cdie=5e-9)), ('D35L4bank', dict(cdie=4.5e-9, bank=dict(cb=0.5, ir=25.0)))]:
    h, sd = kernel(dict(kw, **extra)); ker[name] = h.tolist()
    print(name, f"peak {h.min():.4f} mV/mA @ {h.argmin()} step-min {np.cumsum(h).min():.3f} noise {sd[:20].mean()/abs(h.min()):.2f}", flush=True)
json.dump(ker, open('kernels.json', 'w'))
