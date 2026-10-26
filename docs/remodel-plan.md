# PDN remodel on LACPo + cache-model currents (MediumBoom, 45 nm, 1.1 V, 333 MHz)

Written 2026-10-10 before any sweep. Everything earlier (curmodel, PaperMatchBoomV3, 1 GHz) is superseded.

## Load (from `traces/lacpo_cy100/bench/<wl>_out/tile_P_t.csv`)
| quantity | value | note |
|---|---|---|
| clock | 3 ns (333 MHz) | LACPo training clock |
| VDD | 1.1 V | LACPo/PTPX corner |
| tile mean | 150-235 mW -> 0.14-0.21 A | 10 riscv-tests |
| tile p99.9 | ~420 mW -> 0.38 A (truth likely ~0.47 A) | model tail ~17-25 % low |
| per-cycle step p99.9 | 170-300 mW -> 0.16-0.28 A per 3 ns | model ~40 % low on the extreme tail |
| leakage | 63 mW (35 %) | PTPX median |

## Die and network (real units; assumptions marked)
| item | value | basis |
|---|---|---|
| tile area | 4 mm^2 (sweep 3-6) | ASSUMED: 63 mW leakage / ~50 nW per Nangate45 cell -> ~1.3 M cells x ~2 um^2 / 0.7 util; flop-array caches |
| intrinsic die C | ~3 nF | switched C = 0.115 W x 3 ns / 1.1^2 = 285 pF/cycle; non-switching ~10x |
| explicit MOS decap density | 12 nF/mm^2 (thin oxide) | 45 nm MOS decap 10-15 fF/um^2 |
| decap at 1.8 V (active decap reservoir) | 5 nF/mm^2 (thick oxide) | ASSUMED |
| package to die | L 0.5 nH, R 50 mOhm (sweep 0.2-1 nH) | small wire-bond/flip-chip test package |
| package decap | 1 uF, ESL 100 pH, ESR 20 mOhm | |
| board + VR | off-chip buck 3.3 -> 1.1 V, 2 MHz, loop ~100 kHz | the "one big off-chip path" |
| droop budget | 55 mV (5 %); also report 33 mV (3 %) | |
| die-pkg resonance | 1/(2 pi sqrt(0.5 nH x ~8 nF)) ~ 80 MHz = 0.24 cycles/cycle | inside the current's 0.15-0.5 band (~21 % of variance) |

Small on-die regulators (alternate parallel path), per region: 1.8 V input rail via its own package path, buck with
air-core on-die inductor 1-5 nH, Ron 0.5-4 Ohm, fsw <= 200 MHz, triggered by sensed droop (comparator delay) and/or
the issue-stream predictor (feed-forward lead). Active decap: reservoir C at 1.8 V switched onto the rail on trigger
(Meng/Saleh-style charge dump), recharged after. Regions (area share ~ leakage share): R0 frontend+I$ 0.39,
R1 integer core 0.22, R2 LSU+D$+PTW 0.33, R3 FP 0.06; grid R between regions, Lg = Rg x 10 ps.

## Phases
1. Issue trace from VCD: per-IQ slot request/grant, fu, rob, pdst/prs + rtypes -> uop instances (ready r, grant t), real deps.
2. Per-class current kernels: ridge-fit tile/region current on per-class issue counts (lags 0-12), held-out-benchmark R^2;
   residual = unschedulable background (time-stretched on stalls).
3. Replay (MediumBoom port map, age-order baseline must reproduce real grants): reorder, holds, dummy ops, combos;
   policies driven by the issue-stream current predictor (no PDN dependence -> traces cached); predictor triggers exported.
4. PDN solver in C++ (trapezoidal, dt = 0.25 ns, 4 regions, switch-level regulators and active decap); cross-check vs Pulsim.
5. GA (NSGA-II): genes = per-region decap area + damping, active decap size/threshold, regulator count/L/Ron/trigger/lead,
   scheduling policy variant; objectives = on-die PDN area, energy + slowdown; constraint droop <= budget; fit on
   droop-heavy windows of 5 benches, validate full traces on all 10. Report usable-area gain vs passive-only baseline.
6. Jitter: alpha-power clock-buffer delay vs V (Vt 0.45 V, alpha 1.3), insertion delay 0.5 ns (sweep), period and
   cycle-to-cycle jitter + timing guardband per design.
7. DVFS reaction: 1.0 <-> 1.1 V reference steps under load for off-chip-only, on-die regulated, hybrid; settle time to 1 %.
