"""policy6: Phase 5 closed-loop engine. Circuit-state-aware replay (main contribution) + supporting mechanisms.
Same issue replay and tiled PDN as closedloop6; every new mechanism is a parameter that defaults to OFF, and with all
defaults the run is identical to closedloop6 (gate). Mechanisms:
  A  timing errors (timing6) on the datapath at guardband G_mV; detected errors squash the uop and every dependent that
     already issued; replay modes:  'razor' = re-issue after rep_pen cycles (state-blind baseline)
                                    'state' = replay buffer released by circuit state (forecast, sensed voltage vs the
                                              unit's failure voltage, bank charge, refill lookahead), spaced rep_space
                                              apart; forward progress: release at rep_Dmax, second failure -> protected
                                              by a stretched clock cycle
     A4 error feedback: error rate tightens the token threshold (fb_err), arms banks (err_arm), per-class error counts
        weight the reorder choice (cls_w)
  B2 adaptive clock stretch: detector on the worst-tile voltage (threshold cs_thr mV below nominal, noise, latency
     cs_lat); a trigger lengthens the next cs_len clock periods by (1 + cs_s) (phase-picking stretch, not pulse masking):
     each stretched cycle's dynamic charge is spread over (1 + cs_s) T via a fractional time accumulator (extra PDN steps
     at the reduced current). Uops whose worst-voltage cycle is stretched cannot fail. cs_pred: predictive stretch of the
     cycle where the forecast shows the dip; forced-safe replays stretch their execute cycle.
  B4 scheduler v2 terms (tighten the token threshold or add violations): edge-rate limit on predicted issue current
     (edge_dI over edge_w cycles), second-droop window budget (win_dI), charge debt from the delayed off-chip buck
     current (debt_k), refill lookahead from cache events (look_k, look_arm), 77 MHz tone penalty (tone_thr, tone_k),
     bank readiness (ready_min, ready_k).
  ROB (rob = entries, 0 = off): in-order commit (a uop retires once it and every older uop completed and passed its
     error check, completion = issue + 1 + LAT); a uop cannot dispatch (become ready) until it is within rob entries of
     the oldest unretired uop. Without it a delayed / replayed uop never stalls younger work (slowdown underestimated).
  flush=1: recovery squashes the failed uop and EVERY younger issued uop and freezes issue for rep_pen cycles (Razor /
     resilient-core style); flush=0: selective replay (failed uop + issued dependents only).
  avoid_vm (mV): predictive error avoidance - a picked uop is held (up to maxd) when the forecast voltage over its execute
     cycles (sensed worst V + forecast change) falls below its unit's failure voltage + avoid_vm.
  fc_mode: 'pred' = forecast from predicted issue-coupled current only (original); 'meas' = past from the MEASURED total
     current (activity-counter estimate, one cycle late) and future = committed predicted issue current + last measured
     background, net of the measured LDO-bank current. fc_log=1 records the final pick's forecast for scoring.
  fc_kernel: 'default' = kernel of the default 0.3 mm^2 undamped network (original); 'design' = kernel of this design's
     linear network (decap, damping; the active banks are NOT in it). fc_obs (gain): observer - residual between the sensed worst-tile
     voltage (delayed) and the now-cast made for that cycle corrects the forecast, decaying 0.7^(j + delay).
  fc_mode='fc7': learned forecaster (forecast7.Forecaster7, weights fc7_design / fc7_lag0 / fc7_tag): measured current,
     bank current, sensor, cache-pipeline flags, learned per-class issue history, next two greedy selects. Candidate sets
     are evaluated exactly as f_pre + sum of per-class effect vectors.
  lock_k (cycles): regional heavy-op lockout (power-aware steering rule) - when a heavy uop (predicted charge >= lock_q
     mA*cycles: LD, ST, MUL, FPU, CVT) issues in one region (issue queue: int / mem / fp), heavy uops in the OTHER regions
     wait lock_k cycles (per-uop cap maxd), so the regions do not draw heavy current at the same time.
  Baseline: Razor + error-rate voltage scaling (evs_gain, evs_target): rail reference raised by a slow loop on the error
     rate (errors per 1k cycles), clamped 0..80 mV; the off-chip buck follows at its slew limit."""
from __future__ import print_function
import os, sys, collections, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay3, pdn6 as pdn5, timing6 as TM, closedloop6 as C6
B = replay3.B; K = replay3.K; T = pdn5.T
DEF = dict(th=8.0, kv=0.3, kb=4.0, maxd=4, H=8, th_arm=10.0, arm_len=3, arm_duty=0.2, rel_cap=1, reo=1, dummy_th=None,
           R_LL=0.0, vg=0.0, vlow=40.0, hold=1, oracle=0,
           G_mV=None, W=TM.W, rep_mode='none', rep_pen=5, rep_Dmax=24, rep_th=4.0, rep_space=2, rep_bmin=0.3, rep_vm=5.0,
           seed=1, fb_err=0.0, cls_w=0.0, err_arm=0,
           cs_on=0, cs_thr=None, cs_lat=2, cs_len=4, cs_s=0.25, cs_noise=2.0, cs_pred=None,
           edge_dI=None, edge_w=3, win_dI=None, win_k=0.2, debt_k=0.0, look_k=0.0, look_arm=0, tone_thr=None, tone_k=0.0,
           ready_min=0.0, ready_k=0.0, evs_gain=0.0, evs_target=0.0, rob=0, flush=0, avoid_vm=None, fc_mode='pred', fc_log=0, fc_kernel='default', fc_obs=0.0, fc7_design='banks', fc7_lag0=0, fc7_tag='train7', lock_k=None, lock_q=45.0)
W0 = 2 * np.pi * 77e6 * pdn5.T_CYC; TONE = np.exp(-1j * W0 * np.arange(28))

class Engine(C6.CL):
    def __init__(self, w):
        C6.CL.__init__(self, w)
        z = np.load(B + w + '_cache.npz'); n = self.n
        ed = (z['dc_mshr'][:n] > 0).astype(float); ei = (z['ic_miss'][:n] > 0).astype(float)
        # causal lookahead: a refill is due within the next 8 cycles if a D$ miss happened in the last 14 cycles or an
        # I$ miss in the last 22 (refill lead ~6 / ~14 cycles, 8-beat burst)
        cd = np.convolve(ed, np.ones(15))[:n]; ci = np.convolve(ei, np.ones(23))[:n]
        self.soon = (cd + ci) > 0

    def kernel(self, design):
        arr = lambda k: np.asarray(design.get(k) if design.get(k) is not None else [0.0], float).round(12).tobytes()
        key = (arr('A_x'), arr('R_d'), tuple(sorted((design.get('bank') or {}).keys())), repr(sorted((k, v) for k, v in design.items() if k not in ('A_x', 'R_d', 'bank', 'buck'))))
        cache = self.__dict__.setdefault('_gk', {})
        if key not in cache: cache[key] = C6.pdn_impulse(design)
        return cache[key]

    def run(self, design, pol='A', prm=None, a=0, b=None, warm=600, sensor_delay=2, q_mV=5.0, record=False):
        rp = self.rp; prm = dict(DEF, **(prm or {}))
        b = self.n if b is None else b
        sel = np.where((rp.t >= a) & (rp.t < b))[0]; loc = {g_: k for k, g_ in enumerate(sel)}; m = len(sel)
        iq, fu = rp.iq[sel], rp.fu[sel]; cls = rp.cls[sel]; age = rp.age[sel]; tt = rp.t[sel] - a
        r = (np.where(rp.pin[sel], rp.t[sel], rp.r0[sel]) - a).astype(np.int64); r = np.maximum(r, 0)
        prods = [[loc[p] for p in set(rp.prods[g_]) if p in loc] for g_ in sel]; cons = [[] for _ in range(m)]
        for i, ps in enumerate(prods):
            for p_ in ps: cons[p_].append(i)
        pend = np.array([len(ps) for ps in prods]); delay = np.zeros(m, np.int64); tnew = np.full(m, -1, np.int64); held = np.zeros(m, np.int64)
        ver = np.zeros(m, np.int64); wake = collections.defaultdict(list)
        for i in range(m):
            if pend[i] == 0: wake[int(r[i])].append((i, 0))
        L = 3 * (b - a) + 20000                 # room for up to ~200 % slowdown before 'overflow'
        ring = np.zeros((L + K, self.A5.shape[1])); Ihat = np.zeros(L + K); BG = self.BG[a:b]; SOON = self.soon[a:b]
        sim = pdn5.Sim(design, Im_A=(BG[:min(3000, len(BG))].mean(0) + 0.5 * self.A5[:, :, :].sum() / max(1, self.n) * 0) * 1e-3)
        g = self.kernel(design) if prm['fc_kernel'] == 'design' else self.g; Kg = len(g); H = int(prm['H'])
        dq = max(1, int(sensor_delay))
        vq = collections.deque([np.full(T, pdn5.VNOM)] * dq, maxlen=dq)
        bq = collections.deque([1.0] * dq, maxlen=dq)
        sq = collections.deque([0.0] * dq, maxlen=dq)
        has_bank = bool(design.get('bank')); src = 0; slip = 0.0; ready = []; done = 0; t = 0
        vmin_all = []; arm_left = 0; arm_ewma = 0.0; E_ref = 0.0; E_dum = 0.0; Ilp = None; recov = 0; nhold = 0; nreo = 0; ndum = 0; narm = 0; mu = 0.0
        I_ref = BG.sum(1).mean() + self.A5[sel].sum() / max(1, (b - a))
        orc = bool(prm.get('oracle', 0))
        kv = (lambda i: self.A5[sel[i]].sum(0)) if orc else (lambda i: self.hp[cls[i]])
        # --- Phase 5 state ---
        timing = prm['G_mV'] is not None; rng = np.random.RandomState(int(prm['seed']))
        if timing and prm['rep_mode'] not in ('razor', 'state'): raise ValueError("timing errors need rep_mode 'razor' or 'state'")
        uname = np.array([TM.FU.get(int(c), 'ALU') for c in fu]); lat = np.array([TM.LAT[u] for u in uname])
        VF = {u: TM.v_fail(TM.R[u], prm['G_mV']) for u in TM.R} if timing else {}
        VU = {u: TM.v_fail(TM.R[u], prm['G_mV'], prm['W']) for u in TM.R} if timing else {}
        chk = collections.defaultdict(list); vworst = []; stretched = []; pend_str = set()
        rbuf = []; inrb = np.zeros(m, bool); nfail = np.zeros(m, np.int64); safe = np.zeros(m, bool); next_rel = 0
        nerr = 0; nundet = 0; nsquash = 0; nforced = 0; nstretch = 0; defer = []; err_ewma = 0.0; cls_err = collections.Counter()
        s_short = s_long = 0.0; debt = 0.0; evs = 0.0; evs_sum = 0.0
        LEAK = pdn5.LEAK_T.astype(np.float32)[None, :]
        def s12():                                               # predicted current of the next two greedy selects
            rd = sorted(ready, key=lambda i: age[i]); b0, _ = replay3.port_fill(rd, fu, iq); rem = [i for i in rd if i not in b0]
            ok = lambda x: [i for (i, v) in (wake[x] if x in wake else []) if v == ver[i] and pend[i] == 0 and tnew[i] < 0 and not inrb[i]]
            c1 = sorted(rem + ok(t + 1), key=lambda i: age[i]); s1, _ = replay3.port_fill(c1, fu, iq)
            c2 = sorted([i for i in c1 if i not in s1] + ok(t + 2), key=lambda i: age[i]); s2, _ = replay3.port_fill(c2, fu, iq)
            S = np.zeros(H + 1)
            for i in s1: S[1:] += fc7.hp[cls[i]][:H]
            for i in s2: S[2:] += fc7.hp[cls[i]][:H - 1]
            return S
        def forecast(extra):
            if FC7:
                if fpc[0] is None: fpc[1] = s12(); fpc[0] = fc7.forecast(Ihat[t:t + H + 1], np.zeros(10), fpc[1])
                return fpc[0] + sum((fc7.E[cls[i_]] for (i_, k0) in extra), np.zeros(H))
            fut = Ihat[t:t + H + 1].copy()
            for (i_, k0) in extra: fut[k0:k0 + K][:H + 1 - k0] += kv(i_)[:H + 1 - k0]
            lo = max(0, t - Kg)
            if MEAS:                                              # measured past, background held at its last value
                bg = (Imeas[t - 1] - Ihat[t - 1]) if t > 0 else 0.0
                seq = np.concatenate([np.zeros(Kg - (t - lo)), Imeas[lo:t] - mu_m, fut + bg - mu_m])
                return np.array([g[::-1].dot(seq[j + 1:j + 1 + Kg]) for j in range(H)]) + obs_corr
            past = Ihat[lo:t] - mu
            seq = np.concatenate([np.zeros(Kg - (t - lo)), past, fut - mu])
            return np.array([g[::-1].dot(seq[j + 1:j + 1 + Kg]) for j in range(H)]) + obs_corr
        def rise(extra):                                         # predicted issue-current rise over edge_w cycles (mA)
            w_ = int(prm['edge_w']); fut = Ihat[t:t + 3].copy()
            for (i_, k0) in extra: fut[k0:k0 + K][:3 - k0] += kv(i_)[:3 - k0]
            prev = Ihat[max(0, t - w_):t]
            return float(fut.max() - (prev.min() if len(prev) else 0.0))
        def viol(f, extra, th_eff):
            if f.min() < -th_eff: return True
            return prm['edge_dI'] is not None and rise(extra) > prm['edge_dI']
        def squash(i):
            nonlocal done
            S = set(); st = [i] + ([x for x in order[rank[i] + 1:rank[i] + 513] if tnew[x] >= 0] if prm['flush'] else [])   # younger in-flight uops
            while st:
                x = st.pop()
                if x in S: continue
                S.add(x); st.extend(c for c in cons[x] if tnew[c] >= 0)
            for x in S: tnew[x] = -1; done -= 1
            aff = set(S) | set(c for x in S for c in cons[x] if tnew[c] < 0)
            for x in aff:
                pend[x] = sum(1 for p_ in prods[x] if tnew[p_] < 0); ver[x] += 1
                if x in ready: ready.remove(x)
            for x in aff:
                if x != i and pend[x] == 0: wake[t + (int(prm['rep_pen']) if prm['flush'] else 1)].append((x, ver[x]))
            return len(S)
        str_left = 0; acc = 0.0; t_extra = 0.0; vref_rows = []; rng2 = np.random.RandomState(int(prm['seed']) + 7919)
        ROBN = int(prm['rob']); order = np.argsort(age, kind='stable'); rank = np.empty(m, np.int64); rank[order] = np.arange(m)
        commit = 0; rob_wait = []; nrobstall = 0; freeze = 0; navoid = 0; lock = np.zeros(3, np.int64) - 1; nlock = 0
        MEAS = prm['fc_mode'] == 'meas'; Imeas = np.zeros(L + K); mu_m = 0.0; flog = []
        FC7 = prm['fc_mode'] == 'fc7'; fpc = [None, None]; Ib_last = 0.0
        if FC7:
            import forecast7 as F7
            fc7 = F7.Forecaster7(F7.wfile(prm['fc7_design'], bool(prm['fc7_lag0']), prm['fc7_tag']))
            c7d = dict(F7.cache_sigs(self.w)); c7arr = [c7d[nm] for nm in fc7.sig]
        OBS = float(prm['fc_obs']); fnow = []; vbar = pdn5.VNOM; obs_corr = np.zeros(H); DEC = 0.7 ** np.arange(H)
        while done < m or t < (b - a):
            if t >= L - 2: raise RuntimeError('overflow')
            e_now = 0
            # A: timing check for uops whose execute window ended last cycle
            if timing and (t - 1) in chk:
                for i in chk.pop(t - 1):
                    if tnew[i] < 0: continue
                    t0 = tnew[i] + 2; t1 = t - 1; ws = vworst[t0:t1 + 1]
                    if not len(ws): continue
                    u = uname[i]; wsu = [w_ for k2, w_ in enumerate(ws) if not stretched[t0 + k2]]
                    if safe[i] or not wsu: continue
                    vmn = min(wsu)
                    if vmn < VF[u] and rng.rand() < TM.P[u]:
                        if vmn < VU[u]: nundet += 1; continue                 # beyond the detection window: silent
                        nerr += 1; e_now += 1; cls_err[cls[i]] += 1; nfail[i] += 1; nsquash += squash(i)
                        if prm['flush']: freeze = max(freeze, t + int(prm['rep_pen']))
                        if prm['err_arm'] and has_bank and arm_ewma < prm['arm_duty']: arm_left = int(prm['arm_len'])
                        if prm['rep_mode'] == 'razor': wake[t + int(prm['rep_pen'])].append((i, ver[i]))
                        else: rbuf.append([i, t]); inrb[i] = True
            err_ewma += (e_now - err_ewma) / 512.0
            if prm['evs_gain']: evs = float(np.clip(evs + prm['evs_gain'] * (err_ewma * 1000.0 - prm['evs_target']), 0.0, 80.0))
            evs_sum += evs
            if ROBN:                                              # in-order commit, then admit waiting uops into the ROB window
                while commit < m:
                    i_ = order[commit]
                    if tnew[i_] >= 0 and tnew[i_] + 2 + lat[i_] <= t: commit += 1
                    else: break
                if rob_wait:
                    ok_ = [i for i in rob_wait if rank[i] < commit + ROBN]
                    if ok_:
                        rob_wait = [i for i in rob_wait if rank[i] >= commit + ROBN]
                        ready.extend(i for i in sorted(ok_, key=lambda x: rank[x]) if pend[i] == 0 and tnew[i] < 0 and not inrb[i] and i not in ready)
            if t in wake:
                for (i, v) in wake.pop(t):
                    if v == ver[i] and pend[i] == 0 and tnew[i] < 0 and not inrb[i] and i not in ready:
                        if ROBN and rank[i] >= commit + ROBN:
                            if i not in rob_wait: rob_wait.append(i); nrobstall += 1
                        else: ready.append(i)
            pick = []
            vs = vq[0]; margin = (vs.min() - (pdn5.VNOM - prm['vlow'] * 1e-3)) * 1e3; bank_frac = bq[0]
            if OBS and t > dq:                                    # observer: sensed (age dq) vs now-cast made for that cycle
                vbar += (vs.min() - vbar) / 256.0
                obs_corr = OBS * ((vs.min() - vbar) * 1e3 - fnow[t - dq]) * DEC * 0.7 ** dq
            if FC7:                                               # start of cycle t: measured t-1, sensor t-2, cache flags at t
                tg = t - int(round(slip)); sp = max(0, tg) if t == 0 else max(src, min(src + 1, tg)); sp = min(sp, len(BG) - 1)
                fc7.begin(Imeas[t - 1] if t else None, Ib_last if t else None, vs.min() * 1e3 if t >= 2 else None,
                          [arr_[min(a + sp, len(arr_) - 1)] for arr_ in c7arr])
                fpc[0] = fpc[1] = None
            # A3: circuit-state-aware replay release (one per rep_space cycles)
            if rbuf and t >= next_rel:
                for ent in sorted(rbuf, key=lambda e_: e_[1]):
                    i, td = ent; d_ = t - td
                    if d_ < prm['rep_pen'] or pend[i] > 0 or tnew[i] >= 0 or i in ready: continue   # producer re-squashed: wait
                    go = d_ >= prm['rep_Dmax'] or nfail[i] >= 2
                    if nfail[i] >= 2: safe[i] = True; nforced += 1
                    if not go:
                        f_ = forecast([(i, 0)]); u = uname[i]
                        go = (f_.min() > -prm['rep_th'] and vs.min() > VF[u] + prm['rep_vm'] * 1e-3
                              and (not has_bank or bank_frac >= prm['rep_bmin']) and not SOON[min(src, len(SOON) - 1)])
                    if go:
                        ready.append(i); rbuf.remove(ent); inrb[i] = False; defer.append(d_); next_rel = t + int(prm['rep_space']); break
            if ready and t >= freeze:
                ready.sort(key=lambda i: age[i]); base, _ = replay3.port_fill(ready, fu, iq); pick = base
                if timing and prm['avoid_vm'] is not None and pick:     # predictive error avoidance
                    fa = forecast([(i, 0) for i in pick]); vnow = vs.min(); keep = []
                    for i in pick:
                        j1 = min(H - 1, 1 + int(lat[i])); vp = vnow + (fa[2:j1 + 1].min() - fa[0]) * 1e-3
                        if vp < VF[uname[i]] + prm['avoid_vm'] * 1e-3 and held[i] < prm['maxd']: navoid += 1; held[i] += 1
                        else: keep.append(i)
                    pick = base = keep
                if prm['lock_k'] and pick:                             # regional heavy-op lockout
                    keep = []
                    for i in pick:
                        hv = self.charge[cls[i]] >= prm['lock_q']; rg = int(iq[i])
                        if hv and held[i] < prm['maxd'] and any(lock[r_] > t for r_ in range(3) if r_ != rg):
                            nlock += 1
                            if pol == 'A': held[i] += 1
                            continue
                        keep.append(i)
                        if hv:
                            for r_ in range(3):
                                if r_ != rg: lock[r_] = max(lock[r_], t + int(prm['lock_k']))
                    pick = base = keep
                if pol == 'TOK' and base:
                    extra = 0.0
                    if prm['fb_err']: extra += prm['fb_err'] * err_ewma * 1000.0
                    if prm['win_dI'] is not None: extra += prm['win_k'] * max(0.0, s_short - s_long - prm['win_dI'])
                    if prm['debt_k']: extra += prm['debt_k'] * max(0.0, sq[0])
                    if prm['look_k'] and SOON[min(src, len(SOON) - 1)]: extra += prm['look_k']
                    if prm['tone_thr'] is not None:
                        tn = abs(TONE.dot(Ihat[max(0, t - 24):t + 4][-28:]) if t >= 24 else 0.0) / 28.0
                        extra += prm['tone_k'] * max(0.0, tn - prm['tone_thr'])
                    if prm['ready_k'] and has_bank: extra += prm['ready_k'] * max(0.0, prm['ready_min'] - bank_frac)
                    th_eff = float(np.clip(prm['th'] + prm['kv'] * np.clip(margin, -20, 20) + prm['kb'] * (bank_frac if has_bank else 0.0) - extra, 2.0, 40.0))
                    fb = forecast([(i, 0) for i in base])
                    if fb.min() < -prm['th_arm'] and has_bank and arm_ewma < prm['arm_duty']: arm_left = int(prm['arm_len'])
                    if prm['look_arm'] and has_bank and SOON[min(src, len(SOON) - 1)] and arm_ewma < prm['arm_duty']: arm_left = int(prm['arm_len'])
                    if viol(fb, [(i, 0) for i in base], th_eff):
                        if prm['reo']:
                            ck = (lambda i: kv(i).sum() + prm['cls_w'] * cls_err[cls[i]]) if prm['cls_w'] else (lambda i: kv(i).sum())
                            light = sorted(ready, key=lambda i: (ck(i), age[i]))
                            forced = [i for i in ready if held[i] >= prm['maxd']]
                            sel2, _ = replay3.port_fill(forced + [i for i in light if i not in forced], fu, iq)
                            if len(sel2) == len(base):
                                f2 = forecast([(i, 0) for i in sel2])
                                if f2.min() > fb.min(): pick, fb = sel2, f2; nreo += 1
                        if viol(fb, [(i, 0) for i in pick], th_eff) and prm['hold']:
                            keep = [i for i in pick if held[i] >= prm['maxd']]
                            for i in sorted([i for i in pick if held[i] < prm['maxd']], key=lambda i: kv(i).sum()):
                                ex = [(j, 0) for j in keep + [i]]
                                if not viol(forecast(ex), ex, th_eff): keep.append(i)
                            if len(keep) < len(pick): nhold += 1; recov = 4
                            pick = keep
                    elif recov > 0 and len(pick) > prm['rel_cap']:
                        pick = sorted(pick, key=lambda i: (-held[i], age[i]))[:max(1, int(prm['rel_cap']))]
                    if recov > 0: recov -= 1
                    if prm['cs_pred'] is not None and pick:
                        fp = forecast([(i, 0) for i in pick])
                        if fp.min() < -prm['cs_pred']: pend_str.add(t + int(np.argmin(fp)))
            if pol == 'TOK' and prm['dummy_th'] is not None:
                fd = forecast([(i, 0) for i in pick])
                if fd.max() > prm['dummy_th']:
                    used = {q: [False] * len(p) for q, p in replay3.PORTS.items()}; replay3.port_fill(pick, fu, iq, used)
                    for _ in range(sum(1 for x in used[0] if not x)):
                        ring[t:t + K] += self.dummyA.T; Ihat[t:t + K] += self.hp[0]; ndum += 1
                        E_dum += float(self.dummyA.dot(self.WB).sum()) * 1e-3 * pdn5.VNOM * pdn5.T_CYC
            for i in pick:
                ready.remove(i); tnew[i] = t; done += 1; ring[t:t + K] += self.A5[sel[i]].T; Ihat[t:t + K] += kv(i)
                d = t - tt[i]; slip += (max(d, 0) - slip) / 32.0
                if timing: chk[t + 1 + int(lat[i])].append(i)
                if safe[i]: pend_str.update(range(t + 2, t + 2 + int(lat[i])))
                for c in cons[i]:
                    delay[c] = max(delay[c], d); pend[c] -= 1
                    if pend[c] == 0: wake[int(max(r[c] + max(delay[c], 0), t + 1))].append((c, ver[c]))
            for i in ready:
                if i not in pick and pol != 'A': held[i] += 1
            mu += (Ihat[t] - mu) / 256.0
            if FC7:
                nc_ = np.bincount(cls[pick], minlength=10).astype(float) if pick else np.zeros(10)
                if prm['fc_log']: flog.append(fc7.forecast(Ihat[t:t + H + 1], nc_, fpc[1] if fpc[1] is not None else s12()))
                fc7.end(Ihat[t], nc_)
            if (prm['fc_log'] or OBS) and not FC7:
                f_ = forecast([]); fnow.append(f_[0] - obs_corr[0])
                if prm['fc_log']: flog.append(f_)
            s_short += (Ihat[t] - s_short) / 30.0; s_long += (Ihat[t] - s_long) / 300.0
            tgt = t - int(round(slip)); src = max(0, tgt) if t == 0 else max(src, min(src + 1, tgt)); src = min(src, len(BG) - 1)
            I_t = BG[src] + ring[t].dot(self.WB)
            vref = pdn5.VNOM + evs * 1e-3
            if prm['R_LL'] > 0:
                tot = I_t.sum(); Ilp = tot if Ilp is None else Ilp + (tot - Ilp) / 300.0
                vref = pdn5.VNOM + prm['vg'] * 1e-3 - prm['R_LL'] * (Ilp - I_ref) * 1e-3
            arm = np.full(T, 1 if arm_left > 0 else 0, np.uint8); narm += int(arm_left > 0); arm_ewma += (float(arm_left > 0) - arm_ewma) / 256.0
            arm_left = max(0, arm_left - 1)
            if t >= warm: E_ref += I_t.sum() * 1e-3 * pdn5.VNOM * pdn5.T_CYC
            str_this = str_left > 0 or t in pend_str
            if str_this:                                          # stretched period: same charge over (1 + s) T
                pend_str.discard(t); str_left = max(0, str_left - 1)
                I_t = pdn5.LEAK_T + np.maximum(I_t - pdn5.LEAK_T, 0) / (1 + prm['cs_s']); acc += prm['cs_s']; t_extra += prm['cs_s']
            vmin, vavg, vres = sim.run(I_t[None, :].astype(np.float32), vref=np.array([vref]), arm=arm[None, :], want_res=has_bank, want_state=True)
            xrows = []
            while str_this and acc >= 1.0:
                vx, _, _ = sim.run(I_t[None, :].astype(np.float32), vref=np.array([vref]), arm=arm[None, :], want_state=True)
                xrows.append(vx[0]); acc -= 1.0
            vq.append(np.round(vmin[0] / (q_mV * 1e-3)) * q_mV * 1e-3)
            if has_bank:
                vr = vres[0][[t_ for t_ in design['bank']]]; bq.append(float(np.clip((vr.mean() - 1.3) / (2.15 - 1.3), 0, 1)))
            iL1 = float(sim.state[0][0]); debt += (I_t.sum() * 1e-3 - iL1) * pdn5.T_CYC * 1e9; debt *= (1 - 1 / 1000.0); sq.append(debt)
            vmin_all.append(vmin[0]); vref_rows.append(vref)
            for xr in xrows: vmin_all.append(xr); vref_rows.append(vref)
            vw = float(min([vmin[0].min()] + [xr.min() for xr in xrows])); vworst.append(vw)
            Ib_last = 1e3 * float(sim.state[0][2]); Imeas[t] = float(I_t.sum()) - Ib_last; mu_m += (Imeas[t] - mu_m) / 256.0
            stretched.append(str_this)
            # B2: detector (latency, noise) -> stretch the following cs_len periods
            if prm['cs_on'] and t >= prm['cs_lat']:
                thr = pdn5.VNOM - (prm['cs_thr'] + prm['cs_noise'] * rng2.randn()) * 1e-3
                if vworst[t - int(prm['cs_lat'])] < thr:
                    if str_left == 0: nstretch += 1
                    str_left = max(str_left, int(prm['cs_len']))
            t += 1
        V = np.array(vmin_all); VR = np.array(vref_rows, np.float32)[:, None]; d = (VR[warm:] - V[warm:]) * 1e3; worst = d.max(1)   # droop vs the rail's own target
        st = sim.stats; E = st[0]
        exec_new = (int(tnew.max()) if m else 0) + t_extra; exec_old = int(tt.max()) if m else 1
        res = dict(p999=float(np.percentile(worst, 99.9)), max=float(worst.max()), p99=float(np.percentile(worst, 99)),
                   tile_spread=float(np.ptp(np.percentile(d, 99.9, axis=0))), slowdown=(exec_new - exec_old) / max(1, exec_old),
                   loss=float((st[2] + st[3] + st[4] + st[5]) / E), E=E, E_cmos=float(E * V[warm:].mean() / pdn5.VNOM),
                   energy_over=float((st[2] + st[3] + st[4] + st[5] + E_dum) / max(E - E_dum, 1e-30)), E_off=float(st[1]), vmean=float(V[warm:].mean()),
                   nhold=nhold, nreo=nreo, ndum=ndum, arm_frac=narm / max(1, t),
                   nerr=nerr, nundet=nundet, nsquash=nsquash, nforced=nforced, nstretch=nstretch, stretch_time=t_extra, evs_mean=evs_sum / max(1, t), nrobstall=nrobstall, navoid=navoid, nlock=nlock, nuops=m,
                   defer_mean=float(np.mean(defer)) if defer else 0.0, defer_max=int(max(defer)) if defer else 0)
        if prm['fc_log']: res['flog'] = np.array(flog); res['vavg_t'] = np.array(vworst)
        if record: res['vworst_pdn'] = V.min(1); res['vworst'] = np.array(vworst); res['stretched'] = np.array(stretched); res['tnew'] = tnew
        return res
