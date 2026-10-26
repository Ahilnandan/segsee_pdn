# Protocol results: GA5 Pareto designs on full traces (2026-10-23)

Closed loop (dependency-aware issue replay + tiled 8x8 PDN), all 10 benchmarks, full traces (28k-436k cycles each).
Designs: the all-decap baseline plus up to 4 GA5 Pareto picks per capacitor budget (lowest p99.9, lowest max, lowest
energy, knee), found by NSGA-II on the worst 20k-cycle window of 5 benches. Forecast kernel sign fixed before the GA run.
Droop = worst tile below 1.1 V, mV. "Worst" = worst bench; "mean" = mean over benches. Energy is shown relative to the
all-decap run on the same bench (cancels the known bias of the pdn5 energy formula). Slowdown is mean / worst bench.

## Budget 0.2 mm^2
| design | worst p99.9 | worst max | mean p99.9 | mean max | energy vs base | slowdown |
|---|---|---|---|---|---|---|
| all-decap | 78.8 | 128.9 | 63.0 | 100.7 | 0 | 0 |
| lowest p99.9 (16 banks, holds) | **57.8** | 116.7 | 49.1 | 83.8 | +9.4 % | 0 |
| lowest max (16 banks, holds) | 62.3 | 124.3 | 51.7 | 83.3 | +9.4 % | 0 |
| lowest energy (no banks, holds) | 71.8 | **156.7** | 53.0 | 100.1 | 0 | 0.02 / 0.16 % |

## Budget 0.35 mm^2
| design | worst p99.9 | worst max | mean p99.9 | mean max | energy vs base | slowdown |
|---|---|---|---|---|---|---|
| all-decap | 66.7 | 110.4 | 53.9 | 87.8 | 0 | 0 |
| lowest p99.9 (2 banks) | **44.9** | 94.4 | 41.6 | 73.8 | +3.2 % | 0 |
| lowest max (4 banks) | 50.8 | **84.5** | 46.6 | 66.0 | +4.1 % | 0 |
| lowest energy (no banks, holds) | 64.6 | 124.1 | 48.7 | 82.4 | 0 | 0 |

## Budget 0.5 mm^2
| design | worst p99.9 | worst max | mean p99.9 | mean max | energy vs base | slowdown |
|---|---|---|---|---|---|---|
| all-decap | 56.6 | 100.5 | 46.6 | 76.7 | 0 | 0 |
| lowest p99.9 (4 banks) | 41.8 | 92.4 | 36.0 | 66.8 | +4.5 % | 0 |
| lowest max (4 banks) | **40.3** | 92.9 | 37.2 | 65.3 | +6.6 % | 0 |
| lowest energy (no banks) | 55.7 | 99.2 | 45.9 | 75.3 | 0 | 0 |
| knee (2 banks) | 42.0 | **80.5** | 37.2 | 63.2 | +2.8 % | 0 |

## What holds and what does not
1. **p99.9 gains hold on full traces:** -27 % (0.2), -33 % (0.35), -26..-29 % (0.5) vs all-decap, at +3..9 % energy and
   no slowdown. All good p99.9 designs use banks.
2. **Max gains shrink on full traces.** Windows predicted -38 %; full traces give -9 % (0.2), -23 % (0.35), -20 % (0.5,
   knee). The GA's "lowest max" pick is not the full-trace max winner at 0.2 or 0.5: worst-case droop on benches outside
   the GA windows (towers, spmv) was not seen during the search.
3. **Protocol without banks fails the worst-case guard rail.** The cheap no-bank designs make worst max WORSE than
   all-decap: 156.7 vs 128.9 mV (0.2, rsort / spmv) and 124.1 vs 110.4 mV (0.35, rsort), although their mean droop
   improves. Holds and reordering create release bursts on traces they were not tuned on (checklist 2.4).
4. **Area saving:** 0.35 mm^2 with banks (lowest max: 50.8 / 84.5 mV; lowest p99.9: 44.9 / 94.4) beats 0.5 mm^2 all-decap
   (56.6 / 100.5) on both metrics -> ~30 % less capacitor area. The window result "0.2 with banks beats 0.35 all-decap"
   holds only for p99.9 (57.8 vs 66.7); its worst max is 6 % higher (116.7 vs 110.4).

## Caveats
- Energy uses the pdn5 formula corrected by subtracting the baseline per bench; pdn6 removes the bias (Phase 5 B0).
- Realistic predictor (R^2 0.17); cache current extremes ~40 % low; die area estimated; benchmarks, not stressmarks.
- GA fitness used only 5 benches' worst windows; GA6 should score more benches or full traces for the max objective.

# Phase 5 first full-trace results (2026-10-23; policy6 engine, eval6_A6.json / eval6_B8.json)
## A6: replay at reduced guardband (0.35 mm^2; worst over 10 benches; errors per 1M uops; no undetected errors anywhere)
| design | G | config | worst max | errors/1M | slowdown |
|---|---|---|---|---|---|
| banks | 40 | Razor | 104.6 | 391 | 0 |
| banks | 40 | Razor + error-rate voltage scaling | 87.6 | 72 | 0 (rail +4.6 mV) |
| banks | 40 | clock stretch + Razor | 81.4 | 251 | 2.1 % |
| banks | 40 | state-aware replay | 92.4 | 370 | 0 |
| banks | 40 | OURS (with refill-lookahead arming) | 122.4 | 391 | 0 |
| banks | 30 | Razor | 96.1 | 1430 | 0 |
| banks | 30 | state-aware replay | 109.8 | 1408 | 0 |
| banks | 30 | OURS | 120.9 | 1349 | 0 |
Findings:
1. Any replay lets the guardband drop from ~85-110 mV (max droop) to 30 mV at ~0 slowdown: errors are rare (<= 0.14 %
   of uops) and a 5-cycle replay is cheap. This is the Razor result (prior art), not ours.
2. State-aware replay: errors -2..-5 % vs Razor, droop unchanged within noise. At these rates placement barely matters.
3. OURS bundle made worst max droop WORSE (122 vs 105 mV): refill-lookahead bank arming blocks bank recharge and starves
   the banks (B8: that term alone 84.5 -> 103 mV). Removed in OURS2.
4. Error feedback did not lower the error rate on the banks design (-5 % on all-decap).
5. Razor + voltage scaling looks best on errors and droop, but its 5-13 mV higher rail costs ~1-2.4 % dynamic power that
   energy_over does not count, and droop is measured vs nominal, not the raised target. Must be fixed before comparing.
## B8: scheduler-v2 ablation (banks design, droop only)
- edge-rate limit 15 mA: worst max 84.5 -> 80.5 mV, energy 4.10 -> 3.94 %: the only term that helped.
- 3.6 MHz window, charge debt, 77 MHz tone, bank readiness: results identical to TOK base -> never triggered with the
  first-guess weights (calibration needed).
- refill lookahead with arming: harmful (max 84.5 -> 103). All terms together: 113.
- separate-rail banks: loss 4.10 -> 3.04 %, mean p99.9 -4 %, max ~same.
## Status of the main claim
Not supported yet: circuit-state-aware replay shows no measurable gain over Razor at the error rates of this model.
Next test (A6b, running): guardband 20 / 10 mV and a full-flush replay penalty (20 cycles), where Razor should degrade.

## A6b + apparatus check (2026-10-23)
A6b (10 full traces, G 20 / 10 mV, replay penalty 5 / 20): Razor and state-aware replay both ~0 % slowdown even at
G = 10 mV (1.8 % of uops failing); state-aware errors -1..-4 %; droop unchanged within noise.
APPARATUS FINDING: the issue engine had no ROB capacity / in-order commit, so delayed or replayed uops never stalled
younger work. Added rob=64 + in-order commit (policy6, off by default). Gate: rob=0 == closedloop6 exactly; rob=64
baseline differs slightly (7-9 ROB stalls per run; completion model issue+1+LAT is a little off) -> not ticked.
With the ROB, SELECTIVE replay (failed uop + dependents) is still nearly free (<= 0.1 % at G=10, penalty 20): 64 entries
hide the delay. Every earlier "0 % slowdown" for holds/replay was measured without the ROB and must be re-checked.
FLUSH recovery (Razor / resilient-core style: squash all younger uops, freeze issue 20 cycles), qsort 60k window, banks:
| G | config | errors | slowdown | max |
|---|---|---|---|---|
| 30 | selective replay | 31 | 0.00 % | 59.7 |
| 30 | flush | 20 | 0.12 % | 62.4 |
| 30 | flush + predictive avoidance 5 mV | 20 | 0.06 % | 53.8 |
| 10 | selective replay | 420 | 0.00 % | 59.4 |
| 10 | flush | 440 | 8.16 % | 66.6 |
| 10 | flush + avoidance 5 mV | 395 | 7.89 % | 57.8 |
| 10 | flush + state-aware release | 430 | 7.30 % | 68.4 |
- Recovery type decides the cost: selective replay ~free, flush 8 % at G = 10 mV.
- Predictive avoidance barely prevents errors (-10 % at G=10 while holding 12k uop-cycles): the forecast cannot see the
  droops that cause errors (predictor R^2 0.17). State-aware release with flush: -11 % relative slowdown.
=> In this model the big win is detect + SELECTIVE replay in the OoO core (guardband ~85-110 -> 10-30 mV at ~0 %
   slowdown). Circuit-state awareness adds little because errors are rare, cheap to replay selectively, and not
   predictable with the current forecast. Need a prior-art check on selective replay of timing errors in OoO cores
   before claiming that part.

## A6c: fixed engine (review fixes), ROB 80, learned forecaster fc7 (held-out weights), 10 full traces (2026-10-24)
Held-out benches (qsort, mm, towers), banks design; slowdown = mean (collapse = fell > 20k cycles behind, i.e. >= ~6 %):
| G | recovery | mean slowdown | errors / 1M uops | collapses |
|---|---|---|---|---|
| 30 | selective replay | 0 % | 556 | 0 |
| 30 | flush (Razor) | 0.48 % | 498 | 0 |
| 30 | flush + state-aware release | 0.39 % | 506 | 0 |
| 30 | flush + predictive avoidance 5 mV | 2.63 % | 305 | 0 |
| 30 | flush + error-rate voltage scaling | 0 % | 100 (rail +9 mV) | 0 |
| 10 | selective replay | 0.34 % (worst of 10 benches 1.75 %) | 29330 | 0 |
| 10 | flush (Razor) | 12.3 % | 16261 | 2 of 3 |
| 10 | flush + predictive avoidance | 19.7 % | 12234 | 2 of 3 |
| 10 | flush + voltage scaling | 0 % | 123 (rail +28 mV) | 0 |
Findings:
1. Selective replay wins: guardband 85-110 mV -> 10 mV at <= 0.6 % mean slowdown (held-out), no undetected errors.
2. State-aware release = Razor within noise (errors and slowdown).
3. Predictive avoidance (fc7) prevents 20-50 % of errors but holds 10-35 % of uop-cycles, so flush slowdown RISES
   (0.5 -> 2.6 % at 30 mV); only win: it prevented the single rsort collapse at 30 mV (all 10 benches).
4. Voltage scaling never slows down but raises the rail 9-28 mV (~2-5 % power), giving back most of the guardband gain.
5. Flush at 10-20 mV collapses on rsort / mm / qsort whatever the variant.
=> Circuit-state-aware replay NOT supported in this model (forecast R^2 ~0.3 at 2 cycles is not enough to avoid errors
   cheaply; selective replay already makes errors nearly free). Strong measured result: banks + selective timing-error
   replay in the OoO core. Collapsed configs not re-run with the enlarged buffer yet (exact slowdowns unknown, >= ~6 %).
