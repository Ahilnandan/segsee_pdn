# Power-aware scheduling, charge banks and selective replay for BOOM: final results (2026-10-26)

Written for: the project team. Simulator: LACPo-based per-tile currents of the MediumBoom tile (10 riscv-tests, full
traces) -> tiled PDN (8x8 per core; 2x2-core shared-PDN variant) -> closed-loop issue engine with real ready sets,
ports, dependencies, 80-entry ROB, timing-error detection and selective replay. All numbers are full traces unless noted.

## 1. Result in designer terms (single core; reference = all-decap 0.35 mm^2, no policy, guardband 110.4 mV)
| configuration | guardband | performance at same V | power at same f | slowdown (mean / worst) |
|---|---|---|---|---|
| all-decap, no policy | 110.4 mV | 0 | 0 | 0 |
| charge-pump banks, no replay | 85.7 mV | +2.9 % | -0.5 % | 0 |
| power-aware scheduling (forecast holds), no replay | ~100-102 mV | +0.2..+0.9 % | -1.5..+2 % | 0.4-0.7 % / 1.6-4.5 % |
| **all-decap + selective replay** | **10 mV** | **+11.4 %** | **-17.4 %** | 0.02 % / 0.17 % |
| banks + selective replay | 10 mV | +11.3 % | -14.4 % | 0.12 % / 1.0 % |
| scheduling + selective replay | 10-30 mV | +8.5..+10.8 % | -11..-17 % | 0.4-0.7 % / 1.7-4.5 % |
No undetected errors in any configuration (detection window 25 % of the period). Power view: V^2 dynamic power at the
saved guardband, bank losses added back. Performance view: alpha-power delay law (VT 0.40 V, alpha 1.3).

## 2. Multi-core (2x2 cores, one shared PDN; 60k-cycle windows; mc_full.json)
- Droop is chip-wide on a shared supply (far core sees 90-95 % of a burst, also with PDK-derived grids); only fully
  separate grid AND package rail makes it local, at ~3x higher per-core impedance.
- Burst alignment sets droop: 4x qsort in lockstep 48.3 mV p99.9 vs random offsets 21.3 mV.
- Cross-core coordinator (power-aware scheduling at chip scale: hold a few uops on some cores so the total current rise
  stays within 30 mA over 3 cycles): lockstep 48.3 -> 22.9 mV, errors at 30 mV 184 -> 1 per 1M uops, at 0.19 % mean
  slowdown; minimum guardband 20 -> 10 mV. Barrier jitter 29.6 -> 22.5 mV. Mixed programs: -1..-3 mV, 1.6 % / 3.5 %.
- On the shared die the single-core bank design is worse than plain decap at 10-20 mV guardband (~2 % loss, 3-4x the
  errors); banks would need re-sizing for the shared die.

## 3. What each idea contributed
| idea | verdict | evidence |
|---|---|---|
| selective timing-error replay | **decisive** | guardband 110 -> 10 mV, ~0 slowdown, 0 undetected |
| charge-pump LDO banks | helps droop, not the replay regime | p99.9 -27..-33 % at equal area; +2.9 % perf without replay; adds 1.5-4 % loss |
| separate-rail banks | better banks | 18-40 % less loss than pump-fed, similar droop |
| power-aware scheduling, single core | works on droop, little on guardband | forecast holds p99.9 -15 % (decap) at 0.2-0.7 % slowdown; halves errors under replay but its slowdown outweighs that |
| power-aware scheduling, multi-core | **works for synchronised workloads** | lockstep / barrier droop halved, guardband 20 -> 10 mV at 0.2 % |
| learned voltage forecaster | enabler | R^2 2 cycles ahead ~0.3-0.6 vs ~0 before; made holds work |
| regional heavy-op lockout rule | no effect | heavy ops in different regions almost never collide; droop set by total bursts |
| state-aware replay timing / predictive avoidance | no gain | = Razor; avoidance costs more slowdown than it saves |
| reorder-only, 4 scheduler-v2 terms, refill-lookahead arming | no gain / harmful | see protocol-results.md |

## 4. Caveats (read before quoting numbers)
- Timing-error model is assumed (unit path ratios, sensitising probabilities, window, delay law); no per-path timing
  reports were available. Error rates scale with the assumed probabilities.
- Max droop on one bench is chaotic (+-10-20 mV from any decision change); no-replay guardbands use worst max over 10
  benches and are fragile. p99.9 comparisons are robust.
- Control logic (rename, ROB, issue) is assumed covered by margin; replay covers datapath only.
- Workloads: 10 small riscv-tests (+ a miss-driven stressmark that matches the worst benchmark p99.9, not its max);
  multi-core uses replicated single-thread traces, not real barrier-synchronised code.
- PDN: estimated die area, PDK-derived but not layout-extracted grid, lumped package; LACPo cache model ~40 % low on
  the largest steps; the separate 1.8 V rail is ideal.
- Not done: parallel FIVR in the current model; head-to-head against published schemes; GA re-optimisation (GA6).

Details: protocol-checklist.md (all gates and findings), protocol-results.md (all experiment tables).
