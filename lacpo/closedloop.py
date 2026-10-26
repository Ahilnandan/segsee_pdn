"""Closed-loop co-simulation: dependency-aware issue replay (real ready sets, MediumBoom ports) <-> tiled PDN (pdn5),
stepped together every clock cycle. The issue policy sees: the issue-queue current forecast (data-agnostic kernels),
per-tile die-voltage sensors (delay d cycles, quantised q mV) and bank reservoir levels; it can reorder, hold, insert
dummies, pre-arm the LDO banks (and block their recharge), and move the off-chip reference along a load line.

Policies:  'A' baseline (must be identical to the open-loop replay).  'TOK' token protocol (see token_policy()).
Per-uop current: workload5 attribution (exact in the baseline); background per tile time-warped by the running slip."""
from __future__ import print_function
import os, sys, collections, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay3, pdn5, workload5 as W5
B = replay3.B; HERE = replay3.HERE; K = replay3.K; T = pdn5.T

def pdn_impulse(design=None, n=48):
    """mV of mean die voltage per mA of total load current, 1-cycle pulse (from the tiled network, frequency solve)."""
    net = pdn5.build(design or dict(A_x=0.3 * pdn5.CELL_SHARE)); f = np.fft.rfftfreq(1024, pdn5.T_CYC)
    import scipy.sparse.linalg as spl
    u = np.zeros(net['nu']); u[1:1 + T] = pdn5.CELL_SHARE; b = net['B'].dot(u); Z = np.zeros(len(f), complex)
    for k, fk in enumerate(f):
        if k == 0: continue
        x = spl.spsolve(((2j * np.pi * fk) * net['E'] - net['A']).tocsc(), b); Z[k] = x[net['didx']].dot(pdn5.CELL_SHARE)
    # 1 mA, one-cycle rectangular pulse -> voltage (cycle-averaged)
    P = np.fft.rfft(np.r_[1e-3, np.zeros(1023)]); h = np.fft.irfft(P * Z, 1024) * 1e3
    return h[:n]

class CL(object):
    def __init__(self, w):
        rp = replay3.Replay(w); rp.calibrate(); self.rp = rp; self.w = w
        z = np.load(B + w + '_attr5.npz'); assert len(z['t']) == rp.m and (z['t'] == rp.t).all()
        self.A5 = z['A'].astype(np.float32); self.BG = z['BG'].astype(np.float32); self.n = len(self.BG)
        self.hp = np.load(HERE + 'hpos_blocks.npy').sum(0)                       # (classes, K) predicted mA per uop
        self.charge = self.hp.sum(1); self.WB = W5.WB
        alu = rp.cls == 0; self.dummyA = self.A5[alu].mean(0) if alu.any() else np.zeros_like(self.A5[0])
        self.g = pdn_impulse()

    def run(self, design, pol='A', prm=None, a=0, b=None, warm=600, sensor_delay=2, q_mV=5.0, record_t=False):
        rp = self.rp; prm = dict(dict(th=8.0, kv=0.3, kb=4.0, maxd=4, H=8, th_arm=10.0, arm_len=3, arm_duty=0.2, rel_cap=1, reo=1,
                                      dummy_th=None, R_LL=0.0, vg=0.0, vlow=40.0, hold=1), **(prm or {}))
        b = self.n if b is None else b
        sel = np.where((rp.t >= a) & (rp.t < b))[0]; loc = {g_: k for k, g_ in enumerate(sel)}; m = len(sel)
        iq, fu = rp.iq[sel], rp.fu[sel]; cls = rp.cls[sel]; age = rp.age[sel]; tt = rp.t[sel] - a
        r = (np.where(rp.pin[sel], rp.t[sel], rp.r0[sel]) - a).astype(np.int64); r = np.maximum(r, 0)
        prods = [[loc[p] for p in set(rp.prods[g_]) if p in loc] for g_ in sel]; cons = [[] for _ in range(m)]
        for i, ps in enumerate(prods):
            for p_ in ps: cons[p_].append(i)
        pend = np.array([len(ps) for ps in prods]); delay = np.zeros(m, np.int64); tnew = np.full(m, -1, np.int64); held = np.zeros(m, np.int64)
        wake = collections.defaultdict(list)
        for i in range(m):
            if pend[i] == 0: wake[int(r[i])].append(i)
        L = (b - a) + 20000; ring = np.zeros((L + K, self.A5.shape[1])); Ihat = np.zeros(L + K); BG = self.BG[a:b]
        sim = pdn5.Sim(design, Im_A=(BG[:min(3000, len(BG))].mean(0) + 0.5 * self.A5[:, :, :].sum() / max(1, self.n) * 0) * 1e-3)
        g = self.g; Kg = len(g); H = int(prm['H'])
        vq = collections.deque([np.full(T, pdn5.VNOM)] * (sensor_delay + 1), maxlen=sensor_delay + 1)
        bq = collections.deque([1.0] * (sensor_delay + 1), maxlen=sensor_delay + 1)
        has_bank = bool(design.get('bank')); src = 0; slip = 0.0; ready = []; done = 0; t = 0
        vmin_all = []; arm_left = 0; arm_ewma = 0.0; E_ref = 0.0; Ilp = None; recov = 0; nhold = 0; nreo = 0; ndum = 0; narm = 0; mu = 0.0
        I_ref = BG.sum(1).mean() + self.A5[sel].sum() / max(1, (b - a))
        orc = bool(prm.get('oracle', 0))
        kv = (lambda i: self.A5[sel[i]].sum(0)) if orc else (lambda i: self.hp[cls[i]])
        def forecast(extra):
            fut = Ihat[t:t + H + 1].copy()
            for (i_, k0) in extra: fut[k0:k0 + K][:H + 1 - k0] += kv(i_)[:H + 1 - k0]
            lo = max(0, t - Kg); past = Ihat[lo:t] - mu
            seq = np.concatenate([np.zeros(Kg - (t - lo)), past, fut - mu])
            return np.array([g[::-1].dot(seq[j + 1:j + 1 + Kg]) for j in range(H)])
        while done < m or t < (b - a):
            if t >= L - 2: raise RuntimeError('overflow')
            if t in wake: ready.extend(wake.pop(t))
            pick = []
            vs = vq[0]; margin = (vs.min() - (pdn5.VNOM - prm['vlow'] * 1e-3)) * 1e3; bank_frac = bq[0]
            if ready:
                ready.sort(key=lambda i: age[i]); base, _ = replay3.port_fill(ready, fu, iq); pick = base
                if pol == 'TOK' and base:
                    th_eff = float(np.clip(prm['th'] + prm['kv'] * np.clip(margin, -20, 20) + prm['kb'] * (bank_frac if has_bank else 0.0), 2.0, 40.0))
                    fb = forecast([(i, 0) for i in base])
                    if fb.min() < -prm['th_arm'] and has_bank and arm_ewma < prm['arm_duty']: arm_left = int(prm['arm_len'])
                    if fb.min() < -th_eff:
                        if prm['reo']:
                            light = sorted(ready, key=lambda i: (kv(i).sum(), age[i]))
                            forced = [i for i in ready if held[i] >= prm['maxd']]
                            sel2, _ = replay3.port_fill(forced + [i for i in light if i not in forced], fu, iq)
                            if len(sel2) == len(base):
                                f2 = forecast([(i, 0) for i in sel2])
                                if f2.min() > fb.min(): pick, fb = sel2, f2; nreo += 1
                        if fb.min() < -th_eff and prm['hold']:
                            keep = [i for i in pick if held[i] >= prm['maxd']]
                            for i in sorted([i for i in pick if held[i] < prm['maxd']], key=lambda i: kv(i).sum()):
                                if forecast([(j, 0) for j in keep + [i]]).min() >= -th_eff: keep.append(i)
                            if len(keep) < len(pick): nhold += 1; recov = 4
                            pick = keep
                    elif recov > 0 and len(pick) > prm['rel_cap']:          # release limiter right after holds
                        pick = sorted(pick, key=lambda i: (-held[i], age[i]))[:max(1, int(prm['rel_cap']))]
                    if recov > 0: recov -= 1
            if pol == 'TOK' and prm['dummy_th'] is not None:
                fd = forecast([(i, 0) for i in pick])
                if fd.max() > prm['dummy_th']:
                    used = {q: [False] * len(p) for q, p in replay3.PORTS.items()}; replay3.port_fill(pick, fu, iq, used)
                    for _ in range(sum(1 for x in used[0] if not x)): ring[t:t + K] += self.dummyA.T; Ihat[t:t + K] += self.hp[0]; ndum += 1
            for i in pick:
                ready.remove(i); tnew[i] = t; done += 1; ring[t:t + K] += self.A5[sel[i]].T; Ihat[t:t + K] += kv(i)
                d = t - tt[i]; slip += (max(d, 0) - slip) / 32.0
                for c in cons[i]:
                    delay[c] = max(delay[c], d); pend[c] -= 1
                    if pend[c] == 0: wake[int(max(r[c] + max(delay[c], 0), t + 1))].append(c)
            for i in ready:
                if i not in pick and pol != 'A': held[i] += 1
            mu += (Ihat[t] - mu) / 256.0
            # background time-warp and this cycle's tile currents
            tgt = t - int(round(slip)); src = max(0, tgt) if t == 0 else max(src, min(src + 1, tgt)); src = min(src, len(BG) - 1)
            I_t = BG[src] + ring[t].dot(self.WB)
            vref = pdn5.VNOM
            if prm['R_LL'] > 0:
                tot = I_t.sum(); Ilp = tot if Ilp is None else Ilp + (tot - Ilp) / 300.0
                vref = pdn5.VNOM + prm['vg'] * 1e-3 - prm['R_LL'] * (Ilp - I_ref) * 1e-3
            arm = np.full(T, 1 if arm_left > 0 else 0, np.uint8); narm += int(arm_left > 0); arm_ewma += (float(arm_left > 0) - arm_ewma) / 256.0
            arm_left = max(0, arm_left - 1)
            if t >= warm: E_ref += I_t.sum() * 1e-3 * pdn5.VNOM * pdn5.T_CYC
            vmin, vavg, vres = sim.run(I_t[None, :].astype(np.float32), vref=np.array([vref]), arm=arm[None, :], want_res=has_bank)
            vq.append(np.round(vmin[0] / (q_mV * 1e-3)) * q_mV * 1e-3)
            if has_bank:
                vr = vres[0][[t_ for t_ in design['bank']]]; bq.append(float(np.clip((vr.mean() - 1.3) / (2.15 - 1.3), 0, 1)))
            vmin_all.append(vmin[0]); t += 1
        V = np.array(vmin_all); d = (pdn5.VNOM - V[warm:]) * 1e3; worst = d.max(1)
        st = sim.stats; E = st[0]
        exec_new = int(tnew.max()) if m else 0; exec_old = int(tt.max()) if m else 1
        res = dict(p999=float(np.percentile(worst, 99.9)), max=float(worst.max()), p99=float(np.percentile(worst, 99)),
                   tile_spread=float(np.percentile(d, 99.9, axis=0).ptp()), slowdown=(exec_new - exec_old) / max(1, exec_old),
                   loss=float((st[2] + st[3] + st[4] + st[5]) / E), E=E,
                   energy_over=float((E + st[2] + st[3] + st[4] + st[5]) / max(E_ref * (len(vmin_all)) / max(1, len(vmin_all) - warm), 1e-30) - 1), nhold=nhold, nreo=nreo, ndum=ndum, arm_frac=narm / max(1, t))
        if record_t: res['tnew'] = tnew; res['sel'] = sel
        return res
