"""Energy / efficiency of direct-only (A), FIVR-only (C) and parallel (D) supplies on full traces.
Simulated: conduction losses (Ron, DCR, package R, grid R) via source energies.
Modelled: FIVR switching loss = nph * fsw * E_cyc, E_cyc calibrated so FIVR-only @140 MHz is 90% efficient
(Intel FIVR class, Burton APEC'14); off-chip VR switching/driver loss = 5% of its input (-> ~90% total with
its simulated conduction); the 1.8 V rail is made off-chip at 90%."""
import json, numpy as np
from multiprocessing import Pool
import par_exp
CF = par_exp.CF
CFG = [('A direct-only', dict(case='A'), 0),
       ('C FIVR-only 140', dict(case='B', ideal_in=True, fsw_fivr=140e6, t3=dict(fc=78e6), **CF), 140e6),
       ('C FIVR-only 35', dict(case='B', ideal_in=True, fsw_fivr=35e6, t3=dict(fc=35e6 / 1.8), **CF), 35e6),
       ('D parallel 140', dict(case='D', fsw_fivr=140e6, t3=par_exp.dpar(140e6, False), **CF), 140e6),
       ('D parallel 35', dict(case='D', fsw_fivr=35e6, t3=par_exp.dpar(35e6, False), **CF), 35e6),
       ('D parallel 70 L2x', dict(case='D', fsw_fivr=70e6, t3=par_exp.dpar(70e6, False, 224e-9), **dict(CF, L_fivr=224e-9)), 70e6),
       ('D parallel 35 L4x', dict(case='D', fsw_fivr=35e6, t3=par_exp.dpar(35e6, False, 448e-9), **dict(CF, L_fivr=448e-9)), 35e6)]
def job(a):
    wl, (name, kw, fsw) = a
    from sweep import setup
    from cosim import run_gen2
    from pdnchain import simulate_chain
    inst, prods, cons, age, pin = setup(wl)
    n = int(max(x[7] for x in inst) - min(x[7] for x in inst)) * 2 + 20000
    g = run_gen2(inst, prods, cons, age, pin, 'A', 100, 0, 25.0)
    o = simulate_chain(None, sched=g, nmax=n, cdie=5e-9, warmup=3000, **kw)
    E = o['EN']; T = o['T_meas']
    v = 1e3 * (0.9 - o['vmin'])
    return dict(wl=wl, name=name, fsw=fsw, T=T, E12=-E['E12'], E18=-E['E18'], Eload=E['Eload'], nsw_f=E['nsw_f'],
                p999=float(np.percentile(v, 99.9)))
if __name__ == '__main__':
    jobs = [(wl, c) for wl in ('qsort', 'dhrystone', 'mm') for c in CFG]
    with Pool(15) as p: R = p.map(job, jobs)
    json.dump(R, open('par_energy.json', 'w'), indent=1)
    # calibrate FIVR switching energy per phase-cycle on FIVR-only @140 MHz (eta_FIVR = 0.90 averaged over workloads)
    c140 = [r for r in R if r['name'] == 'C FIVR-only 140']
    Psw = np.mean([r['Eload'] / 0.90 - r['E18'] for r in c140]) / np.mean([r['T'] for r in c140])
    Ecyc = Psw / (16 * 140e6)
    print(f'calibrated FIVR switching energy {Ecyc*1e12:.2f} pJ per phase-cycle -> {Psw*1e3:.2f} mW at 140 MHz, 16 ph', flush=True)
    for r in R:
        pfsw = 16 * r['fsw'] * Ecyc * r['T'] if r['fsw'] else 0.0
        Ein = r['E12'] / 0.95 + (r['E18'] + pfsw) / 0.90 if r['name'][0] != 'A' or True else 0
        if r['name'][0] == 'A': Ein = r['E12'] / 0.95
        if r['name'][0] == 'C': Ein = (r['E18'] + pfsw) / 0.90
        r['eta'] = r['Eload'] / Ein; r['Psw_f_mW'] = pfsw / r['T'] * 1e3; r['Pload_mW'] = r['Eload'] / r['T'] * 1e3
        r['chatter'] = r['nsw_f'] / max(1.0, 2 * 16 * r['fsw'] * r['T']) if r['fsw'] else 0
    for r in sorted(R, key=lambda r: (r['wl'], r['name'])):
        print(f"{r['wl']:9s} {r['name']:16s} | load {r['Pload_mW']:6.1f} mW | system eta {100*r['eta']:5.1f}% | FIVR sw {r['Psw_f_mW']:5.2f} mW | "
              f"chatter x{r['chatter']:.2f} | p99.9 {r['p999']:6.1f} mV")
    json.dump(R, open('par_energy.json', 'w'), indent=1)
