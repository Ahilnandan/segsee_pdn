# Checklist: hotspot-mapped closed-loop PDN + token-aware issue protocol

Rule for every step: a step is DONE only when its gate check passes and the number is written down here.
No sweeps before the Phase-0 sizing table is complete. No conclusion from a policy until the apparatus can show the
effect it is testing (inject a deliberately bad and a deliberately good pattern first).

## Phase 0 — physical inputs (no simulation yet)
- [ ] 0.1 Block list: every BoomTile block with a current source (20 LACPo blocks, I$, D$ data/meta/MSHR, PTW, frontend rest, core glue).
- [ ] 0.2 Block areas from PTPX leakage shares (same library) scaled to the die estimate; flag as estimate; cross-check
          with flop counts implied by cache sizes (16 KB I$ / D$ as flop arrays).
- [ ] 0.3 Floorplan: 8x8 tiles (~250 um), contiguous placement per block, caches as large blocks, IQ/RF/ALU/LSU small and adjacent
          to their consumers; execution ports mapped to distinct tiles (needed for port steering).
- [ ] 0.4 Grid per tile edge: R and L from a stated 45 nm metal stack (sheet R per layer, strap width/pitch) -> table.
- [ ] 0.5 Bumps per tile from bump pitch; package L/R per bump; total must equal the lumped 0.5 nH / 50 mOhm.
- [ ] 0.6 Intrinsic C per tile from each block's switched + non-switching capacitance estimate.
- [ ] 0.7 Current mapping: block current -> tiles by area; issue-coupled current -> tile of the port/unit that executed the uop.
- [ ] 0.8 Sizing table reviewed: currents per tile (mean, p99.9, max step), grid IR at p99.9, bump current density.
      GATE: lumped equivalent of the tiled network reproduces the 4-region impedance curve within 10 %.

## Phase 1 — simultaneous closed-loop simulator
- [ ] 1.1 C++ PDN made steppable per clock cycle with retained state (one call = 12 sub-steps), callable from the replay.
- [ ] 1.2 Elements, each placeable per tile: damped decap; LDO pass device fed by (a) external rail or (b) charge-pump bank;
          role-switchable on-die buck (injector / regulator); gated active decap (deadband); off-chip buck with load line.
- [ ] 1.3 Sensors per tile: delay 1-2 cycles, 5 mV quantisation, offset error; bank charge level per tile.
- [ ] 1.4 GATES: (a) linear part vs exact frequency solve; (b) coarse setting reproduces earlier 4-region numbers;
          (c) charge conservation over a run; (d) open-loop replay reproduces real issue timing (>= 97 % exact);
          (e) closed loop with policy = baseline is bit-identical to open loop.
- [ ] 1.5 Speed check: full qsort trace closed-loop in < 10 min; otherwise port the replay inner loop to C++.

## Phase 2 — token-aware issue protocol (design slowly, one mechanism at a time)
Inputs available to the protocol each cycle: ready uops + ports + predicted per-tile current (1-8 cycles), per-tile voltage
(sensed, delayed), per-tile bank charge, regulator state, DVFS state.
- [ ] 2.1 Predictor: per-tile current forecast from the issue queue; measure accuracy (R^2, event recall/precision) per tile.
- [ ] 2.2 Token definition: per-tile current budget for the next window, refilled from PDN state (voltage margin, bank charge),
          spent by issuing uops whose predicted current lands in that tile.
- [ ] 2.3 Mechanisms (each ablated alone, then combined):
      a. voltage-aware port steering (zero slowdown)
      b. resonance-aware pacing (break ~4-5 cycle rhythms only when the forecast shows build-up)
      c. predictive load line (target V follows predicted current)
      d. PDN-state gating (throttle issue into a tile whose bank/headroom is low; recharge banks in predicted quiet cycles)
      e. joint DVFS (pre-raise rail, regulator role switch, issue reshaping during the transition)
      f. complementary current budgeting (keep the TOTAL current's di/dt under a limit by trading issue between tiles)
- [ ] 2.4 Guard rails: starvation bound (max hold per uop), forward-progress proof, no holds that create release bursts
          (check the worst case, not only p99.9 — holds made rsort/memcpy worse before).
- [ ] 2.5 GATE: apparatus test — synthetic resonant pattern vs smooth pattern at equal mean must separate clearly
          before any policy result is believed.

## Phase 3 — Pareto optimisation at fixed decap budgets
- [x] 3.1 Budgets: 0.2 / 0.35 / 0.5 mm^2 total capacitor area (decap + bank + active).
- [x] 3.2 Objectives: droop (p99.9 and max, worst tile), DVFS settle time, energy overhead, slowdown.
- [x] 3.3 Variables: per-tile placement/sizing of every element, protocol parameters, load-line slope.
- [x] 3.4 Fitness windows: worst events per (bench, policy) re-found in the tiled model; constraints checked there.

## Phase 4 — validation and reporting
- [ ] 4.1 Full-trace validation, all 10 benches; per-tile droop maps before/after.
- [ ] 4.2 Jitter and timing guardband per tile; DVFS step tests.
- [ ] 4.3 Sensitivity: grid strength x0.5/x2, die area 3-6 mm^2, sensor delay 0-3 cycles, predictor degraded.
- [ ] 4.4 Prior-art check of the winning protocol (search snippets only, no patent-page fetches).
- [ ] 4.5 Write-up: Pareto fronts per budget, mechanism ablation table, hotspot maps; mark older results superseded.

## Phase 0 results (2026-10-22)
0.1-0.3 DONE (floorplan.py -> floorplan.json): 22 blocks, areas from PTPX leakage shares (median report), die 4 mm^2 est.,
8x8 tiles of 250 um. I$ 26.5 %, D$ data 25.2 % (equal 16 KB flop arrays -> consistent), BPD+BTB 9.3 %, FP 8.1 %, int RF 5.5 %,
PTW 5.0 %, ALU0 3.4 %, LSU 2.6 %, IQs/rename/ROB ~1-1.5 % each. Layout: I$ cols 0-1 | BPD/BTB/fetch col 2 | rename/ROB/IQ/RF
col 3 | ALU0/ALU1/FP col 4 | LSU/D$ ctrl col 5 | D$ data cols 5-7 | PTW col 7 bottom.

0.4 Grid (stated assumptions): power mesh on the two top metals (FreePDK45-class M9/M10, ~2 um thick Cu, ~0.011 Ohm/sq),
straps 8 um wide at 40 um pitch per net per direction (20 % coverage) -> effective 0.055 Ohm/sq per direction; one tile
pitch = one square; + lower-metal/via share -> R_edge = 0.08 Ohm between adjacent tile centres, L_edge = R_edge x 10 ps = 0.8 pH.
0.5 Package: flip-chip area array, bump pitch 200 um -> 100 bumps on 2 x 2 mm, 40 VDD / 40 GND / 20 I/O -> 0.625 VDD bumps per
tile. Per bump 0.25 nH / 20 mOhm (bump + package via); bumps in parallel 6 pH; package plane + balls + board vias lumped
0.45 nH / 40 mOhm from the package node to the bump field. Total ~0.456 nH / ~40.5 mOhm (vs 0.5 nH / 50 mOhm lumped before).
0.6 Intrinsic C: 3 nF split by tile area share of cells (leakage share), ESR via tau 0.3 ns.
0.7 DONE (tilecurrents.py): per-tile currents for all 10 benches; core glue spread over core tiles. Hotspots (qsort):
BPD counter-table tiles mean 8 mA / p99.9 19.5 / step 13 mA; core column 5-7 / 10-14 / 6-9 mA; D$ data 2 / 10.8 / 7.7 (refills).
0.8 GATE PASSED: tiled network with the old lumped package + stiff grid reproduces the pdn4 4-region impedance within 4.3 %
(peak 843 vs 852 mOhm @ 72 MHz).
1.1 DONE: pdn5.cpp/pdn5.py — sparse (Eigen SparseLU) trapezoidal, steppable with retained state; 19k cycles/s for 64 tiles.
1.4a GATE PASSED: time domain vs exact frequency solve corr 0.99999, max diff 0.16 mV.

### FINDING (2026-10-23): droop is a global mode on this die
Per-tile droop with the realistic flip-chip package and grid (R_edge 0.08 Ohm): all 64 tiles within 0.4-1.1 mV of each
other (p99.9 and max), every bench, every decap level. Even with a 12x weaker grid (1 Ohm/edge) and ALL decap crammed into
the I$ columns, the spread is only 7 mV (p99.9) / 10 mV (max). Tile currents are only 5-20 mA, so grid IR is ~1 mV;
the droop is the package/die resonance seen by the whole die.
=> Consequences for the protocol: per-tile voltage steering and per-tile decap placement are second-order here (<~10 %);
the lever is the TOTAL current waveform (resonance band), PDN state (bank charge, regulator mode) and DVFS coordination.
Complementary current budgeting (f) reduces to shaping the total current. Keep the tiled model (it proves this and stays
correct if die/current scale up), but design the protocol around global current shaping + bank/regulator coordination.
G (charge-pump bank, open-loop policies) validated: GA picked aggressive holds (H4); full traces miss max 80 mV at
0.43-0.68 mm^2 (82-88 mV); feasible only >= 0.92 mm^2 -> no gain over passive 0.90. Open-loop holds + bank is not enough.
1.x DONE workload5.py: block-level per-uop attribution (9 issue-coupled blocks incl. D$ hit path), identity error 0, movable 6-24 %.
1.4d GATE PASSED: closed-loop baseline grants identical to open-loop replay (100 %), vvadd.
1.4e GATE PASSED: closed-loop baseline droop == open-loop pdn5 on the same currents (55.21 / 75.45 mV both).
1.5 speed: closed loop 12k cycles/s (Python issue engine + per-cycle C++ PDN).
2.5 GATE PASSED: same 40 mA mean extra load on the core column: 4-cycle on/off -> 53.4 mV droop, smooth -> 4.1 mV.
2.1 predictor (class kernels): R2 0.17 vs actual issue-coupled current, top-1 % event recall 0.06 / precision 0.06,
although issue-coupled current is 63 % of total tile-current variance -> predictor is the bottleneck; evaluate the
protocol with the realistic predictor AND an oracle predictor (true per-uop current) to bound the value of better prediction.
2.x first protocol evaluation (p2_eval.json; 5 benches, 20k-cycle worst-event windows, 0.3 mm^2 damped decap):
 - charge-pump LDO banks (16 tiles, 0.6 nF): p99.9 -15..-24 mV, max -20..-45 mV, energy +4.3..7 % (mostly amp bias
   0.4 mA x 16 banks) -> few banks suffice (global mode); bank count/bias must be optimised.
 - token protocol, realistic predictor: p99.9 0..-6 mV, max mixed, slowdown <= 1 %; neutral once banks are present.
 - oracle predictor: p99.9 -7..-11 mV, but holds cost 5-18 % slowdown and release bursts worsen max; reorder-only oracle mixed.
 - arm duty limiter fixed bank starvation (rsort max 99 -> 62 mV with banks). Load line (R_LL 0.15) no benefit.

## Corrections and Phase 3-4 progress (2026-10-24)
### CORRECTION: forecast kernel sign (closedloop.pdn_impulse)
- Found by code review: the kernel drove the load with the injection sign, so the forecast predicted a voltage RISE for a
  load increase. Every TOK / TOK-oracle result from p2_eval.json and the first GA5 run (archive_buggy_kernel/) used it.
  Those results are INVALID and superseded.
- Fixed (closedloop.py:18). Verified against the sim (simdiag2 / kcheck in archive_buggy_kernel/): pulse response matches the
  kernel in sign and shape; white noise 15 mA: corr 0.978, measured/kernel slope 0.80 (48-tap truncation, unresolved).
- The earlier "corr -0.97" from the review came from the same sign error and is not reproducible on the fixed kernel.
- Note on the sim: loads below the leakage floor (dyn = max(0, I - leak)) are clipped to leakage; tests must start above it.

### Phase 3 DONE (2026-10-23, fixed kernel): GA5 budgets 0.35 / 0.2 / 0.5, pop 36 x 30 gens; fronts 31 / 18 / 36.
### Phase 4.1 full-trace validation DONE except per-tile maps -> results in docs/protocol-results.md.
    p99.9 -26..-33 % holds; max only -9..-23 %; no-bank protocol designs FAIL the worst-case guard rail (worst max
    +12..22 % vs all-decap) -> 2.4 not passed; GA6 must score worst case on more benches / full traces.
### Phase 3 (history: restarted on the fixed kernel)
- GA5 restarted (run_ga5.sh): budgets 0.35, 0.2, 0.5 sequentially, 5 workers. Status in ga5_runner.log.
- 3.1-3.4 pending GA completion.

### Phase 4 (in progress)
- 4.1 val5.py written: front designs (up to 4 reps: min p99.9, min max, min energy, knee) + all-decap baseline, all 10 full
  traces, closed loop, spawn pool 6, one bench cached per worker. Queued by run_val5.sh after ALL_DONE -> val5_runner.log.
- 4.2 DVFS (dvfs5.py, libpdn5b.so, role-split on-die buck, 1.10 -> 1.00 -> 1.10 V under qsort load, 8 core tiles, Ipk 50 mA):
    off-chip loop only : enter +-25 mV at 4.61 us (down), 4.70 us (up); steady error -7.8 / +7.5 mV
    role-2 regulator   : 4.30 us (down), 3.39 us (up); steady error -8.6 / +2.2 mV
    role-1 injector    : does NOT move the rail (injects only when below reference; off-chip in current-share mode)
  Finding: on-die regulator speeds the DVFS transition ~7 % (down) / ~28 % (up); both are limited by the assumed 20 mV/us
  off-chip slew (DEFAULT slew in pdn5.py), which needs a sensitivity check (4.3).

## Phase 5 — MAIN CONTRIBUTION: circuit-state-aware replay (decided 2026-10-23; starts after val5 finishes)
Claim: timing errors caused by supply droop are detected on the datapath and re-executed at a time and in a way chosen
from the LIVE CIRCUIT STATE (forecast voltage, bank charge, buck current / charge debt, ring phase, clock-stretch status),
and the error flags themselves feed back into the scheduler. Everything else in this phase supports this.
Prior art (baselines, not novelty): Razor / resilient-circuit detect-and-replay (immediate, state-blind), error-rate
voltage scaling (slow loop to the regulator), pipeline damping, reactive throttling, charge-budget gating, proactive
boost, adaptive clock stretching. Novelty to defend: the replay decision is driven by PDN state, and errors drive the
issue scheduler (fast loop), not only the regulator.
New code goes in new files (pdn6 / closedloop6 / ga6); do not edit files the running val5 imports.

### 5.A Circuit-state-aware replay (core)
- [x] A0 Prior-art check (2026-10-23, search snippets only, 3 queries; no claims read). Disclosed:
         - Razor I/II, Bubble Razor, Bowman (Intel) EDS + instruction replay; replay cost several cycles.
         - Error rate -> supply voltage / clock (slow loop): Razor DVS, US7480838 (ROB holds results until verified,
           controller lowers f / raises V, refetch), US9727342 (clock from per-stage error counts), Bowman adaptive clock.
         - US7480810 "Voltage droop dynamic recovery": droop SENSOR -> pipeline flush -> reissue with inserted STALLS to
           delay. CLOSEST: discloses delayed reissue after a droop (fixed stall delay, sensor-triggered, not error-triggered).
         - US8327118: scheduler suppresses issue to an error-detecting unit on consecutive cycles (short-path errors).
         Not found: replay timing/placement chosen from PDN state (bank charge, charge debt, ring phase, forecast);
         replay-storm avoidance by current shaping; per-cycle choice replay vs hold vs bank arm vs stretch; error flags
         driving the ISSUE scheduler's forecast / pattern learning / bank arming.
         => Claim must be narrowed to the state-driven decision + placement + scheduler feedback; "delay replay after a
         droop" alone is anticipated by US7480810. Full claims review still needed before filing.
- [x] A1 FIRST VERSION DONE (timing6.py; assumptions in its header: alpha-power VT 0.40 / ALPHA 1.3 -> 0.96 % delay per
         10 mV; unit path ratios R, sensitising probabilities P, latencies; window W). Result, budget 0.35, 5 GA windows,
         expected errors per 1M uops:
           guardband (mV)        80    60    50    40    30    20
           all-decap (max 110)   0.6   73   345   836  1955  3727
           GA min-p99.9 (max 94) 0     0.6  9.0   209  1707  4451
         Undetected = 0 everywhere, even with a 10 % window (droop never passes the detection floor here).
         => With replay, the guardband can drop from the max droop (~95-110 mV) to ~40-50 mV at a few hundred errors per
         1M uops; banks cut the error rate ~4-40x at 40-50 mV. Rates scale linearly with the assumed P (sweep needed);
         windows are worst-event windows (pessimistic). Replay cost not yet simulated (A2).
- [ ] A1b Remaining: sweep R / P / W / delay law; full traces; per-unit breakdown.
- [ ] A1 (original spec) Timing-error model: per uop, voltage during its execution vs its path delay (alpha-power delay law, stated
         per-unit path delays); detection window (errors beyond it are UNDETECTED -> the floor that banks / stretch /
         guardband must hold); datapath only (ALU, MUL, FPU, AGU); control state and SRAM are out of scope.
- [ ] A2 Replay mechanism in the issue engine: flagged uop squashed before commit and returned to the ready list;
         dependents re-woken; stated replay cost.
- [ ] A3 State-aware replay decision (the novelty), using the circuit state:
      - WHEN: re-issue only when the forecast voltage has recovered, bank charge is sufficient, charge debt is low,
        no refill burst is due (lookahead: D$ miss ~6 cy, I$ miss ~14 cy)
      - WHERE: place replays into current dips / ring peaks so they smooth the current; spread them (no replay storms)
      - WHICH: per cycle choose the cheapest of issue-and-replay vs hold vs arm bank vs accept a clock stretch
      - forward progress: repeated failure -> hold until recovered or stretch
- [ ] A4 Error feedback: error flags recalibrate the forecast threshold (absorbs temperature, aging, process) and learn
         which instruction patterns cause errors; also trigger bank arming.
- [ ] A5 Gates: replay never makes droop worse than no replay; immediate replay reproduces Razor-style behaviour;
         synthetic droop beyond the window is reported undetected; forward progress holds on a worst-case pattern.
- [ ] A6 Main result: performance at equal reliability and energy, at reduced guardband, vs
         (1) guardband only, (2) clock stretching, (3) Razor-style immediate replay, (4) Razor + error-rate voltage scaling,
         (5) our scheduler without replay, (6) ours with state-aware replay. At 333 MHz and GHz-scaled.

### 5.B Supporting work (needed by 5.A)
- [x] B0 DONE (2026-10-23): pdn6.cpp/.py + closedloop6.py (new files; pdn5/closedloop untouched for the running val5).
         Fixed: role-1 edge energy in J (5 pJ x (Vin/1.8)^2); off-chip role-1 integrator clamped 0..5 A; sensor_delay is the
         true age (min 1 cycle); energy_over = (device losses + dummy energy) / useful load energy; np.ptp. State export
         per cycle: off-chip buck A, on-die buck A, LDO-bank A.
         GATE PASSED: towers + qsort 20k windows, all-decap and banked design, policy A and TOK: p99.9 / max / p99
         identical to closedloop+pdn5; all-decap energy_over exactly 0 (was -0.11..-0.25 %); banked designs 1.7-3.6 %
         (old formula 1.5-3.3 %, it credited the higher mean voltage). Sensor-delay fix: no change on these windows.
- [ ] B1 Simulator state exports: off-chip + on-die buck current per cycle (charge debt), bank charge, ring phase;
         per-band droop readout (77 MHz / 3.6 MHz / buck loop, sums to total).
- [ ] B2 Adaptive clock model: detector threshold/offset/noise, 1-3 cycle reaction, stretch amount (stretch shifts time).
- [ ] B3 Jitter + timing model: clock insertion delay x delay law -> period / cycle-to-cycle jitter, net timing loss.
- [ ] B4 Scheduler v2 base: state observer (model + delayed sensor); edge-rate limit over 1-3 cycles (worst droops are
         69-89 % fast edges; rhythm only in vvadd/rsort/dhrystone); minor tone penalty; 3.6 MHz window budget;
         back-to-back bank readiness; actions cheapest first (reorder, redistribute, bank, dummy, hold, stretch).
- [ ] B5 Separate-rail bank option with full loss accounting (pump, pass device, rail, bias); compare vs pump-fed.
- [ ] B6 Current redistribution: frontend current vs fetch activity model, then trim other units ahead of heavy uops.
- [ ] B7 Gates: forecast vs sim (pulse, white noise); burst train suppressed / smooth untouched; back-to-back dt 5-1000 cy;
         redistribution lowers total di/dt; stretch triggers correctly; slow ramp -> ~0 period jitter.
- [ ] B8 Experiments: ablation one mechanism at a time; head-to-head vs prior art at equal area; scheduling vs damping
         at 77 MHz (R_d > 0.6 Ohm allowed); sensitivity (predictor, sensor delay, grid, die size, clock speed);
         which units cause the fast edges.
- [ ] B9 GA6: genes for replay policy + scheduler weights, bank source, dummy_th; widen bank tau, R_pass, R_d after a
         buildability check. Objectives: performance at equal reliability, energy, worst band share, back-to-back droop.
         Then full-trace validation vs GA5 fronts and all-decap.

### Phase 5 implementation status (2026-10-23)
Engine: lacpo/policy6.py (all mechanisms default OFF), analysis6.py (B1 bands, B3 jitter), timing6.py (A1), pdn6 (B5),
eval6.py (A6 / B8 on 10 full traces, results eval6_<set>.json).
- [x] A2 replay mechanism: detected error squashes the uop + every dependent already issued (transitive), pend recount,
      wake versions invalidate stale wakes; razor mode re-issues after rep_pen = 5 cycles.
- [x] A3 state-aware release: replay buffer, released one per rep_space cycles when forecast dip > -rep_th, sensed V >
      unit failure V + 5 mV, bank charge >= rep_bmin, no refill due; forward progress: release at rep_Dmax = 24, second
      failure -> protected by a stretched execute cycle.
- [x] A4 error feedback: error-rate tightens the token threshold (fb_err), arms banks (err_arm), per-class error counts
      in the reorder key (cls_w).
- [x] A5 gates PASSED (towers, qsort 20k windows): defaults == closedloop6 exactly; timing on with no possible errors ==
      no-timing run; every razor/state run completes (forward progress). "Replay never worse": pending full traces.
- [x] B1 bands + B3 jitter: synthetic gates PASSED (77 MHz -> band 1, 3.6 MHz -> band 2, DC/ramp -> band 3; slow ramp ->
      0 period jitter; 20 mV 77 MHz ripple -> ~11 ps period jitter).
- [x] B2 clock stretch: FIRST MODEL (pulse masking, leakage-only cycles) WRONG - the current drop rang the PDN (p99.9
      57.8 -> 83.1 mV). Replaced by phase stretch: cs_len periods lengthened by (1 + cs_s), same charge spread; now
      neutral-to-helpful for droop (0-1.8 % time). Predictive stretch (cs_pred) implemented.
- [x] B4 terms implemented: edge-rate limit, 3.6 MHz window budget, charge debt (off-chip buck current export), refill
      lookahead (dc_mshr / ic_miss, causal), 77 MHz tone penalty, bank readiness. Weights untuned (first guesses).
- [x] B5 separate-rail bank (pdn6 src=1, Vx 1.8 V, Rx 5 Ohm, charging loss counted). First look (20k windows):
      better droop AND lower loss than pump-fed (qsort 37.1/51.1 vs 42.1/55.3 mV, 1.28 vs 1.66 %). Separate rail's own
      PDN and bump cost NOT modelled yet.
- [x] Baseline Razor + error-rate voltage scaling (evs_gain / evs_target) implemented.
- [ ] B6 current redistribution via fetch throttling: NOT implemented (needs the frontend current model); only
      issue-side redistribution exists (edge limit + refill lookahead).
- [ ] A6 / B8 full-trace runs: running (eval6.py).
- [ ] B9 GA6.
- [ ] A6 main result: NOT SUPPORTED YET (see protocol-results.md "A6b + apparatus check"). Decision needed on the claim.
- [ ] ROB gate: rob=64 baseline not bit-identical (7-9 stalls); fix completion model. Re-check all hold/replay slowdowns with ROB.
- [x] ROB gate PASSED at rob=80 (64 entries + slack for non-unique trace ages): baseline bit-identical on towers, qsort, mm.
- [x] Forecast accuracy measured (policy6 fc_log; R^2 vs actual worst-tile deviation, 0 / 2 / 4 cycles ahead):
      original 'pred' forecast R^2 0.00-0.14 at all horizons (every TOK decision so far used a near-useless forecast);
      'meas' (measured past current, net of bank current) + per-design kernel + observer 0.5: decap 0.49-0.53 / 0.17-0.25 /
      0.08-0.10, banks 0.17-0.25 / 0.07-0.14 / 0.01-0.03. Future current is the limit. BUG FIXED: kernel was built for the
      default network, not the run's design (fc_kernel='design').
- [ ] Agents running (2026-10-23): forecast7.py (better forecaster), calib6.py (B4 term calibration), frontend7.py (B6
      frontend model), policy6 code review.
- [x] Code review of policy6 (agent, 2026-10-23) -> FIXED: (1) state replay could issue a uop twice / before its producer
      when a replay-buffer uop's producer was squashed (inrb flags, guarded release, wakes skip buffered uops); (2) observer
      had no effect in 'pred' mode; (3) avoidance window off by one (now execute cycles t+2..t+1+LAT); (4) predictive
      stretch one cycle late; (5) a stretched worst cycle exempted the whole window (now: any unstretched cycle below VF
      fails); (6) safe replays stretch every execute cycle; (7) extra stretch PDN steps logged in time order, not folded
      into the main row; (8) droop measured vs the rail's own reference (EVS fairness) + E_cmos (I ~ V) reported;
      (10) kernel cache key on full arrays; stretch noise on its own RNG stream.
      GATES PASSED again: defaults == closedloop6 exactly, rob=80 == defaults (towers, qsort; A and TOK). Stress (P=0.6,
      G=5 mV): state replay, selective and flush -> every uop issued, 0 dependency violations (old code crashed).
      Flush + 20-cycle penalty at P=0.6 still overflows: a real Razor collapse, not a bug (no forward-progress guarantee
      for razor flush).
      Not yet fixed (noted): meas-mode forecasts are ~3.6x larger, so GA thresholds (th_arm) tuned on 'pred' fire
      differently -> retune per forecast mode; pending checks at loop end not evaluated (2-3 per run).
- [x] B6 frontend model (agent; frontend7.py + frontend7_model.npz): frontend = I$ + BPD + BTB + fetch, 25 tiles, 65-79 mA
      mean, 49.9 mA idle. Held-out fit R^2 median 0.75 (0.23-0.88). Fetch gating removes ~14 mA in the first cycle,
      17-20 mA/cycle after (peak 19 mA, 0.26 nC for 4 cycles) = ~10-15 % of issue-coupled p99 current; free only when the
      fetch buffer is full (34 % of burst cycles; 4-10 % on fetch-bound qsort/dhrystone/multiply/median), else ~1 cycle
      delay per gated cycle. Interface FrontendThrottle(bench).gate(t, g). => a modest lever.
- [x] Learned forecaster (agent; forecast7.py + 8 weight files) integrated as policy6 fc_mode='fc7' (default off; gates
      still pass). Closed-loop scores on held-out benches reproduce the agent's offline numbers, e.g. 2 cycles ahead R^2:
      decap qsort 0.21 -> 0.37 / mm 0.17 -> 0.32 / towers 0.25 -> 0.57; banks qsort 0.14 -> 0.26 / mm 0.07 -> 0.21 /
      towers 0.07 -> 0.46; worst-1 % recall roughly doubles. Key inputs: cache-pipeline flags (background current) and a
      learned per-class issue kernel; old class kernels explain little. Forecast scale is shrunk 0.5-0.9x -> thresholds
      need retuning (GA6).
- [ ] A6c running: fixed engine, rob=80, fc7 (held-out weights; qsort/mm/towers are the honest benches), G 30/20/10,
      selective vs flush recovery, state release, predictive avoidance, EVS.
- [x] B4 calibration (agent; calib6.py, calib6_eval.json; banks, TOK, rob=80, 10 full traces): the window / debt / tone /
      readiness terms DID fire (10-25 % of decisions with B8 weights) but almost never flipped a decision: th0 sits at its
      18.3 mV ceiling on 86 % of decisions (margin clipped), the old 'pred' forecast's 1st percentile is only -14.9 mV vs
      ~46 mV real p99.9 droop. Calibrated to fire on ~2 % of decisions: hold=0 -> mean p99.9 change <= 0.1 mV; hold=1 ->
      each term -0.05..-0.55 mV mean, worst-bench p99.9 WORSE by 0.3-1.2 mV. => window, debt, tone: DROP; readiness:
      marginal. all4 + edge 15 mA + holds: mean p99.9 -3.5 mV but 6.2 % worst slowdown and an rsort RUNAWAY (>20k cycles
      behind) -> guard rail 2.4 FAILS for edge+holds; holds need a duration / throughput cap.
      Note: single-bench max droop shifts +-10-20 mV from any decision change (chaotic event alignment) -> differences in
      max below ~20 mV between policies on one bench are not meaningful; use p99.9 and multi-bench means.

## Phase 6 — multi-core (2026-10-26): pdn7.py (2x2 copies of the 8x8 tile on one shared PDN), mc_open.py
- [x] Sizing: 16 mm2 die, 16x16 tiles (250 um), 160 VDD bumps (0.625 / tile, as before), segmented package plane (one
      0.45 nH / 40 mOhm branch per core + lateral links L_lat 0.1 nH), off-die parts scaled for 4 cores (board L, R / 4;
      C1, C_pk x 4; buck L1, R1 / 4), intrinsic C 4 x 3 nF. GATE PASSED: per-core-equivalent |Z(f)| identical to the
      1-core model at 0.1-200 MHz.
- [x] Locality: load on core 0 only -> far core sees 95 % of core 0's |Z| at 70 MHz (L_lat 0.1 nH); even with NO lateral
      package coupling and a 6x weaker grid it is 81 %. The continuous on-die grid (8 parallel edges ~0.06 Ohm) couples
      cores more than the package branches (0.45 nH ~0.2 Ohm at 70 MHz) isolate them -> first droop is chip-wide on a
      shared grid; spatial steering stays second-order. Local droop needs separate per-core grids / domains.
- [x] Burst alignment across cores (open loop, 60k cycles, per-core decap 0.35): 1 core qsort 50.8 / 75.8 mV;
      4x qsort lockstep 51.2 / 76.2 (total-current step p99.9 608 mA); 4x qsort random offsets 22.0 / 31.3 (222 mA);
      4 different programs 25.9 / 32.4. => alignment alone changes droop 2.3x. Cross-core coordination (stagger bursts)
      is the large lever for the scheduling scheme. Banks per core help little here (random offsets 27.9 vs 22.0: less
      decap, banks rarely reach their 19 mV offset).
- [ ] Closed-loop multi-core coordinator (agent): 4 issue engines in lockstep on one pdn7 PDN; scenarios lockstep /
      barrier-like jitter / random; coordinator = global current-rise budget with per-core holds.
- [x] CORRECTION (grid check, 2026-10-26): the claim "the continuous grid dominates the coupling" was WRONG. Locality at
      70 MHz (far-core / near-core |Z|): baseline 0.95; VDD+VSS (grid + bumps x2) 0.93; grid 10x / 25x weaker 0.91 / 0.90;
      separate per-core grids with shared package 0.91; separate grids AND no lateral package link 0.14. Cores couple
      through BOTH grid and package plane; grid values barely matter. Full isolation (FIVR-per-core-like domains) makes
      droop local but raises each core's own |Z| ~3x (546 vs 184 mOhm): sharing lets cores borrow neighbours' C.
      => shared supply: total current matters (cross-core staggering is the lever); isolated domains: per-core
      scheduling / banks matter; spatial steering is not the lever in either case.
- [x] PDK grid (Nangate45 tech LEF + OpenROAD reference grid_strategy-M1-M4-M7, downloaded): M1 0.38, M2-3 0.25, M4-6
      0.21, M7-8 0.075, M9-10 0.03 Ohm/sq (old assumption M9/M10 0.011: 2.7x optimistic). Per 250 um tile edge: PDK flip-chip
      grid (reference M1/M4/M7 + M9/M10 8 um @ 40 um) 0.133 H / 0.149 V Ohm (VDD), ~0.27/0.30 with VSS; PDK reference grid
      alone (no M9/M10, peripheral-pad style) 1.24 H / 24.5 V (x2 with VSS). pdn6/pdn7: R_edge_h / R_edge_v added (default =
      old R_edge, old results unchanged). Results (qsort 60k, 1 core): flip-chip PDK grids change droop <= 2 mV, tile spread
      <= 1.3 mV, 4-core locality 0.92-0.94 -> conclusions hold. Reference grid: decap design 63.3/87.0 mV, but the 4-bank
      design COLLAPSES (50.3/58.2 -> 88.9/129.9): banks cannot feed the die through 49 Ohm edges -> on weak grids banks must be
      many and spread. Team's actual Innovus PDN not available locally.
- [x] Multi-core coordinator (agent; mc_engine.py, mc_eval.py, mc_eval_decap/banks.json; windows 70k-130k; gate: 1 core
      == policy6 'A' exactly, with and without rob=80): global current-rise budget over 3 cycles, per-uop maxd 4, per-core
      hold-rate cap 25 %. Budget 30 mA (decap): lockstep 4x qsort p99.9 48.3 -> 22.9 mV at 0.19 % mean / 0.25 % worst
      slowdown (random-offset reference 21.3); barrier jitter 0-20 cy 29.6 -> 22.5 at 0.17 %; 4 different benches 27.6 ->
      25.2 at 1.6 % / 3.5 % worst; banks design similar (lockstep 49.0 -> 28.9). => works for synchronised parallel
      workloads; little for mixed ones. mc_open.py window moved to 70k (load-phase contamination); conclusions unchanged
      (lockstep 48.8 vs random 21.3 mV).
- [x] Single-core steering with the learned forecaster (agent; steer7.py / steer7.json; 10 full traces, rob=80, fc7 held-out
      weights; thresholds re-set to fc7 scale: act on worst ~3-4 % of forecast dips, kv=0). Held-out (qsort, mm, towers)
      mean p99.9: decap A 52.5 -> fc7 holds 44.5 mV (-8.0, -15 %) at 0.19 % mean / 0.23 % worst slowdown; banks 45.8 ->
      42.9 (-2.9) at 0.11 %. In-sample benches agree (decap -8.4..-9.1, up to 5.5 % worst slowdown; banks -2.9..-3.2).
      Old forecast with holds: decap -3.6, banks -0.6. Reorder-only stays ~0 (banks +0.3..+1.6 worse). => steering WORKS
      when the forecast is good, via holds; reorder alone is useless.
