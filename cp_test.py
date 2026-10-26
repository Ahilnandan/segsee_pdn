import numpy as np, importlib, sys
import pdnchain, pdnchain_bak4
n = 4000; I = np.full(n, 100.0)
for k0 in range(1000, 3600, 500): I[k0:k0+150] = 400.0
base = dict(case='B', fsw_fivr=140e6, L_fivr=2e-9, cdie=4.5e-9, ctl='vm')
# regression: behavioral bidirectional bank identical to the pre-edit file
a = pdnchain.simulate_chain(I, bank=dict(cb=0.5, ir=25.0, harvest_v=True), **base)
b = pdnchain_bak4.simulate_chain(I, bank=dict(cb=0.5, ir=25.0, harvest_v=True), **base)
print('regression max |dv| =', np.abs(a['vmin'] - b['vmin']).max())
def rep(tag, o):
    v = 1e3 * (0.9 - o['vmin']); E = o['E']
    print(f"{tag:14s} droop {v.max():6.1f} p99.9 {np.percentile(v,99.9):6.1f} over {1e3*(o['vmax'].max()-0.9):6.1f} | "
          f"Ein {E['Ein']*1e12:7.1f} pJ Estored {E['Estored']*1e12:7.1f} pJ (eta {E['Estored']/max(E['Ein'],1e-30):.2f}) Ediss {E['Ediss']*1e12:7.1f} Erch {E['Erch']*1e12:7.1f}")
rep('none', pdnchain.simulate_chain(I, bank=dict(cb=0.5, ir=25.0), **base))
rep('beh-store', a)
rep('beh-dissip', pdnchain.simulate_chain(I, bank=dict(cb=0.5, ir=25.0, sink_dissipate=True), **base))
rep('beh-hybrid', pdnchain.simulate_chain(I, bank=dict(cb=0.5, ir=25.0, hybrid=True), **base))
for cf in (0.05, 0.1, 0.2):
    c2 = dict(base, cdie=(4.5 - 2 * cf) * 1e-9)
    rep(f'cir-store {cf}', pdnchain.simulate_chain(I, bank=dict(cb=0.5, ir=25.0, circuit=dict(cf=cf, mode='store')), **c2))
    rep(f'cir-hyb {cf}', pdnchain.simulate_chain(I, bank=dict(cb=0.5, ir=25.0, circuit=dict(cf=cf, mode='hybrid')), **c2))
