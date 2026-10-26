"""Closed-loop scheduler <-> circuit co-simulation.

The scheduler is replay2.run() itself, turned into a generator at import time: its idealized
token-bucket bank block (open-loop model) is replaced by
        tok = yield float(I[t])
so each cycle the circuit (pulsim_chargepump's validated bank circuit + optional harvest path)
consumes the load current the scheduler just produced, simulates that cycle, and sends back
tokens = charge the REAL reservoir can deliver right now, from SENSED voltages:
        tok = Cb * (v_bank - v_release_floor),  v_release_floor = min(v_die, -TH_RELEASE)
This removes the model/circuit mismatch found in harvest_circuit.jsonl (open-loop tokens credited
harvested charge the comparator-triggered bank didn't deliver when needed).
"""
import inspect, textwrap, sys, numpy as np, pulsim as p
import replay2
from pulsim_chargepump import build_circuit, V_HEAD, TH_RELEASE_V, CB_F

def _make_gen():
    src = textwrap.dedent(inspect.getsource(replay2.run))
    a = src.index("        tok=min(float(Qb),tok+Ir)")
    b = src.index("        t+=1", a)
    src = src[:a] + "        In[t]=I[t]\n        tok = yield float(I[t])\n" + src[b:]
    src = src.replace("def run(", "def run_gen(", 1)
    ns = dict(vars(replay2)); exec(compile(src, 'replay2.run_gen', 'exec'), ns)
    return ns['run_gen']
run_gen = _make_gen()

def _make_gen2():
    # same scheduler, feedback = (sensed tokens, supply-state signal)
    src = textwrap.dedent(inspect.getsource(replay2.run))
    a = src.index("        tok=min(float(Qb),tok+Ir)"); b = src.index("        t+=1", a)
    src = src[:a] + "        In[t]=I[t]\n        tok, SIG = yield float(I[t])\n" + src[b:]
    src = src.replace("def run(", "def run_gen2(", 1)
    ns = dict(vars(replay2)); exec(compile(src, 'replay2.run_gen2', 'exec'), ns)
    return ns['run_gen2']
run_gen2 = _make_gen2()

def cosim(inst, prods, cons, age, pin, pol='TD', TH=100, ir_mA=25.0, harvest=False, THF=150.0, eta=0.85,
          SMAX=150.0, cyc_s=1e-9, sub=4, dummy=None):
    dt = cyc_s / sub
    tr = np.array([x[7] for x in inst]); nmax = int(tr.max() - tr.min()) * 2 + 20000
    gen = run_gen(inst, prods, cons, age, pin, pol, TH, 0, ir_mA, dummy=dummy, THF=THF)
    st = {'k': -1, 'I': float(next(gen)) * 1e-3, 'vd': 0.0, 'vb': V_HEAD, 'done': None, 'hist': [0.0] * 5, 'hq': 0.0}
    b = build_circuit(st['I'])
    vid, bid = b.node_id_of('vdie'), b.node_id_of('bank')
    nunk, nsw = b.graph.num_branches, b.graph.num_switches
    DR = []
    def tokens():
        floor = min(st['vd'], -TH_RELEASE_V)
        return max(0.0, CB_F * (st['vb'] - floor)) * 1e12
    def advance(k):
        # cycle boundary: scheduler gets SENSED tokens, returns next cycle's load current
        while st['k'] < k:
            st['k'] += 1
            if st['k'] == 0: continue
            st['hist'] = st['hist'][1:] + [st['I'] * 1e3]
            if st['done'] is None:
                try: st['I'] = float(gen.send(tokens())) * 1e-3
                except StopIteration as e: st['done'] = e.value; st['I'] = replay2.FLOOR * 1e-3
            DR.append(-9.0)
    def obs(t, x):
        st['vd'] = float(x[vid]); st['vb'] = float(x[bid])
        if DR and -st['vd'] > DR[-1]: DR[-1] = -st['vd']
    def sw(t):
        m = p.SwitchStateMask(nsw); m.set(0, (st['vd'] < -TH_RELEASE_V) and (st['vb'] > st['vd'])); return m
    def bex(t):
        advance(int(t / cyc_s + 1e-9))
        out = [0.0] * nunk; out[vid] = st['I']
        out[bid] = -(ir_mA * 1e-3) if st['vb'] < V_HEAD else 0.0
        if harvest and st['vb'] < V_HEAD:
            sink = min(SMAX, max(0.0, max(st['hist']) - THF - st['I'] * 1e3)) * 1e-3
            if sink > 0:
                g = eta * (0.9 + st['vd']) / (0.9 + st['vb'])
                out[vid] += sink; out[bid] -= g * sink; st['hq'] += g * sink * dt
        return out
    def keep_going():
        return st['done'] is None or st['k'] < st['done']['makespan'] + 200
    p.simulate(b, t_end=nmax * cyc_s, dt=dt, b_extra_fn=bex, switch_fn=sw, step_observer=obs,
               start_from_dc_op=False, store_every=10**7, should_continue=keep_going)
    R = st['done']
    d = np.array(DR) * 1e3
    return dict(droop=float(d.max()), p999=float(np.percentile(d, 99.9)), slow=100 * (R['makespan'] / R['treal'].max() - 1),
                nhold=R['nhold'], hq_nC=st['hq'] * 1e9)

if __name__ == '__main__':
    import json
    from multiprocessing import Pool
    from sweep import setup
    def job(a):
        wl, ir, h = a
        inst, prods, cons, age, pin = setup(wl)
        r = cosim(inst, prods, cons, age, pin, ir_mA=ir, harvest=h)
        r.update(wl=wl, ir=ir, harv=h); return r
    jobs = [(wl, ir, h) for wl in ('qsort', 'dhrystone', 'mm') for ir in (25, 10) for h in (False, True)]
    with Pool(8) as pl, open('cosim.jsonl', 'w') as f:
        for r in pl.imap_unordered(job, jobs): f.write(json.dumps(r) + '\n'); f.flush(); print(r, flush=True)
