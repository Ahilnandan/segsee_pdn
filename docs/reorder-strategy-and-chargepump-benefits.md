# Power/droop-aware issue reorder + active decap: strategy and validated benefits

Status as of 2026-09-30. Machine: `PaperMatchBoomV3Config` (chipyard, INT IQ 48/w4,
MEM 24/w2, FP 24/w2, ROB 384). Apparatus and all numbers below are reproducible
from `pdn_work/` — file references throughout.

---

## 1. Goal

Cut the PDN guardband (the passive decap + FIVR margin reserved to survive
worst-case di/dt) by (a) reordering issue so the PDN sees fewer/smaller current
edges, and (b) adding a small switched-capacitor "active decap" bank that
supplies a burst of current exactly when an edge is about to happen. Baselines
we're positioned against: P1 pipeline damping (slew-limit issue, no active
decap), US7937563, US11204766, P13 switched decap, P3 current ramping.

## 2. Reorder strategy

### 2.1 Three policies, tested on real hardware ready-sets

Everything runs on `BACT_CAND` — a per-cycle trace of BOOM's *entire* ready
set for each issue queue (not just what got granted), carrying the real
dependency graph (rob_idx/pdst/prs1-3) and, as of this session, each uop's
power token (§4). Policy **A** (oldest-first) exactly reproduces hardware
issue times on all three workloads (0 uops moved) — this is the correctness
gate every other policy is checked against
([`replay-must-use-real-ports-and-deps`](../../../.claude/projects/-home-ahilnandan-Documents-PMIC/memory/replay-must-use-real-ports-and-deps.md)).

| Policy | Rule |
|---|---|
| **A** | oldest-ready-first (= real BOOM grant) |
| **T** | same issue *count* as A, but pick the subset that minimizes the bank's projected charge need over the next H=6 cycles (pure reorder, work-conserving) |
| **TD** | T, then additionally **hold** a uop back (up to MAXD=6 cycles, then force-issue) if even the best subset still exceeds available tokens |

### 2.2 Result: reordering alone does nothing; holding does everything

On real traces (mm/qsort/dhrystone), **T is statistically indistinguishable
from A** — full-event-service% moves ±1–5 points, droop is sometimes *worse*,
despite the picker genuinely choosing differently on 1,000–9,000 cycles per
run. The explanation (not proven, but consistent with every trace): moving an
expensive op by 1–2 cycles just relocates it inside the same 5-cycle
detection window: the edge still happens, so the bank/PDN sees the same rise.

**TD (hold-when-tokens-run-out)** is what actually works: it caps *when*
things are allowed to issue, not just *which order*, effectively turning the
scheduler into a physical rate limiter on the bank's real recharge rate.

| Workload | TD+bank, TH100 slowdown | droop | ≈ equivalent passive decap | damping-only (no bank), same slowdown budget |
|---|---|---|---|---|
| mm | 1.93% | 106 mV | 11.7 nF | 25.09% slowdown for the same droop |
| dhrystone | 1.03% | 85.7 mV | 10.7 nF | 13.01% |
| qsort | 0.49% | 100.6 mV | 12.6 nF | 10.34% |

(Actual on-die decap = 5 nF; all droop numbers are ±30–40% absolute per the
approximate current model, §6 — ratios, not absolutes, are load-bearing.)
Full sweep: `pdn_work/sweep_table.txt`. Source memory:
[`bank-token-scheduling-result`](../../../.claude/projects/-home-ahilnandan-Documents-PMIC/memory/bank-token-scheduling-result.md).

**Bottom line on reorder**: the scheduling lever is *when to hold*, not *what
order to pick among the ready set*. A pure "pack HIGH-with-HIGH" reordering
policy — which is what most literature on this topic proposes — measurably
does not help on this machine; a token-aware admission gate does, and cuts
the slowdown needed for a given droop target by roughly 5–13× versus plain
slew-limiting.

## 3. Architecture

```
reservoir bank Cb (0.5 nF) ── SW_release (comparator-gated, current-controlled)
        │                              │
        │ recharge, Ir mA,             ▼
        │ current-limited,           Vdd grid (vdie)
        │ self-gated off at            │
        │ full charge                  C_pdn (5 nF, passive)
        ▼                              │
   1.8 V FIVR input rail             gnd

token counter (issue-select stage) mirrors the bank's charge in real time;
gates issue when a predicted rise would exceed available tokens.
```

Recharge times (`pdn_work/recharge.py`, 250 pC bank): 10/25/50 ns at 25/10/5 mA
off the 1.8 V rail. A 2:1 SC doubler off Vdd is possible but a *hard*-started
one self-droops 37.4 mV (comparable to half the droop it's meant to prevent);
soft-started (20 mA cap, 20 ns ramp) it self-droops only 1.3 mV.

## 4. Power-token granularity (Chisel + RTL)

**Decision** (explicit ask this session): amplitude/duration levels derived
from the already-validated per-class current shapes in `curmodel.py`, not a
new guess — an 8-bit token per uop, `amp[7:4]|dur[3:0]`: amp = peak mA
quantized in 5 mA steps (0–14), level 15 reserved for the 330 mA CSR/fence
spike class; dur = shape length in cycles − 1. Generated programmatically
from the Python table into a Chisel ROM (`power-tokens.scala`) so hardware
and the Python replay can't drift apart.

Implemented in the real chipyard checkout, `PaperMatchBoomV3Config`:
- `BACT_CAND` extended with the 2 token fields (purely additive; a teammate's
  parser on the same shared checkout reads only the first 6 fields and is
  unaffected — see [`outoforder-team-repo-conflict`](../../../.claude/projects/-home-ahilnandan-Documents-PMIC/memory/outoforder-team-repo-conflict.md)).
- A real hardware token-bucket register in `core.scala` (`BACT_TOK`,
  `+boom_tok_qb`/`+boom_tok_ir` runtime-tunable), same-cycle no-lookahead —
  it can only *observe* when demand exceeds supply, which is the honest
  hardware-realizable half of the python replay's H=6-cycle-lookahead
  admission control.

Verified against curmodel.py exactly (CSR→15/5, ALU→4/3, shift→8/3, branch→9/3
on every sampled instance of a real `hello.riscv` run) and a full clean
**PASSED** dhrystone run (1.22M cycles). Details:
[`hardware-power-tokens`](../../../.claude/projects/-home-ahilnandan-Documents-PMIC/memory/hardware-power-tokens.md).

## 5. The charge pump: is it actually implemented?

**Before this session: no.** Only a token-bucket *number* existed (both in
the Python replay and the RTL accounting register above) — it assumes
instant, lossless delivery with no capacitor voltage sag. That's an
admission-control model, not a circuit.

**Now: yes**, in Pulsim (native C++23 MNA solver, cross-validated against the
project's closed-form exact-ZOH PDN solver to 1.3–2.3% peak-droop agreement
on real traces — [`pulsim-pdn-cross-validation`](../../../.claude/projects/-home-ahilnandan-Documents-PMIC/memory/pulsim-pdn-cross-validation.md)).
Circuit: 0.5 nF reservoir precharged 0.25 V above nominal Vdd, released
through a comparator-gated switch (fires only on real droop, never backfeeds),
recharged through a current-limited path that self-shuts-off at full charge.
(First attempt used a plain diode for release — wrong: an unregulated diode
between a reservoir held above the rail and the rail has nothing holding it
up, it just bleeds down instantly and becomes an always-on current source,
not a triggered reservoir. Fixed by gating release on an actual droop
comparator, matching how a real circuit would do it.)

### 5.1 Per execution unit, not lumped

A single pooled current trace hides which physical block needs the
reservoir. Decomposed by the *real* `fu_code` each instruction issued on
(BOOM's own opcode→hardware-unit routing — not a guessed grouping),
verified to sum back to the lumped total with **zero error**
(`pdn_work/perunit.py`). dhrystone, 20,000 cycles, same bank point as §2.2:

| unit | mean mA | peak mA | events | equivalent-decap benefit |
|---|---|---|---|---|
| ALU | 127.4 | 243.9 | 1006 | +1.6 nF |
| **CSR** | 4.2 | **333.9** | 13 | **+12.3 nF** |
| MUL | 5.0 | 58.9 | 0 | none — never crosses the trigger |
| DIV | 5.4 | 23.9 | 0 | none — never crosses the trigger |
| MEM | 106.3 | 199.9 | 1035 | +6.0 nF |

CSR needs the *most* equivalent decap despite the lowest mean current and the
rarest events (13 in 20k cycles) — a single 330 mA serializing-drain spike
dominates. MUL/DIV need nothing at all on this workload; a pump sized for
them here would be wasted silicon. Spot-checked the dominant unit (ALU) with
the real circuit against the idealized token-bucket math: 71.1 mV vs 74.6 mV
— within 5%, so the earlier token-bucket benefit claims are not optimistic
fiction. Details: [`chargepump-circuit-and-perunit-result`](../../../.claude/projects/-home-ahilnandan-Documents-PMIC/memory/chargepump-circuit-and-perunit-result.md).

### 5.2 Reading "benefit"

`benefit (nF)` = how much *extra passive* decap would be needed to match what
the active circuit achieved. Positive is good; it means the active circuit
is doing more than its own silicon footprint would as plain decap. As a
multiplier over the 0.5 nF physically spent: **ALU 3.2×, MEM 12×, CSR 24.6×**.
Zero (MUL/DIV here) isn't a failure — it means the trigger never fires, so
neither the benefit nor the switching-power cost is paid.

## 6. Cross-check against published silicon: Meng & Saleh, JSSC 2009

*"An Improved Active Decoupling Capacitor for 'Hot-Spot' Supply Noise
Reduction in ASIC Designs"* — a 90 nm 1 V-core test chip, latch-based
comparator-triggered active decap, fabricated and measured (not simulation).
Different mechanism from §5 (theirs **redistributes** existing charge by
switching two decaps from parallel to series for a momentary local boost, no
external recharge; ours **adds** new charge from an external current-limited
rail) — not apples-to-apples on the circuit, but a useful literature anchor
for design-space bounds.

**Independently re-derived every numeric example in their gain analysis**
(not taken on faith):

| quantity | their eq. | their stated value | independently recomputed |
|---|---|---|---|
| ideal gain, n=2, k=0.1 | (4) | 2.75 | **2.750** ✓ |
| k at which G=1, n=2 | (4) | 1/3 | **0.3333** ✓ |
| loss factor b = n·f(R_on)·g(C_decap) | (6) | ≈1.7 | **1.710** ✓ |
| lossy gain | (5) | 2.0 | **2.025** ✓ |
| final settled voltage factor a | (7) | 1.3 | **1.310** ✓ |

All five check out exactly — no arithmetic errors in the paper.

**Update — resolved by looking at the actual figure, not the OCR'd text.**
An earlier pass of this document reported being unable to reconstruct Fig. 1(c)'s
topology with confidence. Looking at the figure image directly found the fix:
Fig. 1(b) explicitly relabels the top rail **"2V_DD"** (not "V_DD") in the
series state — the tell that one outer plate must *disconnect and float*
during boost, not stay pinned to the rail (my text-only reading had assumed
both outer plates stay hard-wired, which cannot produce a boost — see below).
With that correction, an independent circuit-level derivation from Kirchhoff's
laws confirms the claim, not just their algebra:

Two isolated floating nodes at the series-mode transition: `OUT` (cap_L's top
plate alone) and `MID` (cap_L's bottom + cap_R's top, merged by the
cross-coupled switch). Charge is conserved on each independently through the
transition:
```
Q_OUT frozen at C*V_DD (cap_L's charge the instant before isolation)
  => V_OUT = V_MID + V_DD                                    (i)
Q_MID conserved at 0 (cap_L's -C*V_DD bottom-plate charge + cap_R's
  +C*V_DD top-plate charge, from being independently parallel-charged)
  => V_MID = V_OUT / 2                                        (ii)
```
Solving (i) and (ii) together: **V_OUT = 2·V_DD** exactly, matching Fig. 1(b)'s
label. Energy check: before the merge, two caps each at V_DD store `C·V_DD²`
total; after, cap_L sees `V_DD` and cap_R sees `V_DD` (same magnitude, energy
unchanged) — lossless in this ideal `R_on→0` limit, exactly consistent with
their eq. (6) (`b = n·f(R_on)·g(C_decap) → n` as `R_on→0`). Reproducible:
`pdn_work/docs/` derivation script inline above.

### 6.1 Their circuit, run on our workload

Built the validated topology in Pulsim (`pdn_work/pulsim_mengsaleh.py`):
2×0.25 nF decap (0.5 nF total — same silicon budget as our own bank, for a
fair comparison), R_on=15Ω (their quoted 10–20Ω design range), comparator
threshold 50 mV (their stated design sensitivity). One necessary
approximation: the paper's off-state delivery path is really subthreshold
MOSFET leakage, fit from their own unpublished HSPICE curves (Fig. 2) — we
approximate it as one fixed linear resistor and checked the result is
qualitatively stable across a 10× range of that resistance (15–150 Ω) before
trusting it.

Same 6000-cycle windows as §5's real-circuit spot-checks, per execution unit:

| unit | no-help droop | **our bank** | **Meng & Saleh circuit** | eq nF: none / ours / M&S |
|---|---|---|---|---|
| ALU | 95.72 mV | 71.11 mV | 84.35 mV | 9.0 / 18.2 / 11.4 |
| CSR | 68.31 mV | 45.75 mV | **69.56 mV (worse)** | 19.6 / 35.4 / **19.0** |
| MEM | 99.92 mV | 84.65 mV | 95.16 mV | 8.5 / 11.3 / 9.1 |

As a multiplier over the identical 0.5 nF physically spent: **our bank 18.4×
(ALU), 31.6× (CSR), 5.6× (MEM); their circuit 4.8× (ALU), −1.2× (CSR, net
harm), 1.2× (MEM).** Their technique is real and helps ALU and MEM — but on
CSR it makes the worst-case droop slightly *worse* than doing nothing, and
that result held (got mildly more negative, not less) across the whole 15–150
Ω sensitivity sweep, so it isn't parameter-tuning noise.

**Why CSR specifically**: looked at the actual waveform, not just the peak
number. At the true worst CSR event (a single 330 mA serializing spike), the
circuit works exactly as intended — droop at that instant drops from 68.3 mV
to 8.4 mV, an 8× improvement. But right after, the circuit *overshoots*: -17,
-37, -41, -30 mV (voltage rising above nominal) as the caps' abrupt
recombination rings against the lightly-damped PDN (ζ≈0.125 — the same
resonance flagged in [`pdn-solver-numerical-damping`](../../../.claude/projects/-home-ahilnandan-Documents-PMIC/memory/pdn-solver-numerical-damping.md)).
That ringing creates a *new*, unrelated worst-case droop event elsewhere in
the run that didn't exist before, and the whole-run peak-droop metric reports
that new event, not the (genuinely fixed) original one. This is exactly the
risk the paper itself designs around — latch-based comparators with
hysteresis specifically "to prevent the comparator from switching back to
the standby state" on small variations (Section II.C) — on a rare, sharp,
single-shot spike like CSR's, on a real irregular workload (not their clean
periodic 500 MHz lab test), a cruder hard-switched trigger can lose that
fight. Directly reinforces §5.1: this technique needs per-event, possibly
per-unit tuning, not one fixed threshold applied everywhere.

**Reading this against §5's mechanism**: our bank *adds* new charge from an
external recharge rail, so it's not charge-limited by a one-time
redistribution and keeps helping across repeated events. Their circuit
*redistributes* existing charge with a fixed budget — real benefit, but
bounded, and vulnerable to self-induced ringing on sharp single-shot events
that their own controlled lab test (periodic clock, single buffer load)
wouldn't have exposed. Neither result invalidates the other; they're
different tools for different current profiles.

**Still not independently re-verified**: their specific 90 nm f(R_on)/g(C_decap)
curves (Fig. 2 — their own HSPICE result, needs their transistor models) and
the silicon-measured 10–20% noise reduction / 200 MHz–1 GHz sweep (real chip
measurements, not something to re-simulate from a paper). But the core
circuit mechanism — switching two independently-charged decaps into series
delivers a real voltage boost via ordinary charge conservation, no exotic
physics — is now independently confirmed, both algebraically and at the
circuit-topology level.

**Design lessons that do transfer to our per-unit work**:
- **Sizing has a sweet spot, not "bigger is better."** Their Fig. 12/13: active
  beats passive only for capacitor area between ~0.001 mm² and ~0.6 mm²,
  with the biggest advantage at 0.01–0.1 mm². Too small: not enough switched
  charge to matter. Too large: the fixed-speed switching circuitry can't
  slew the oversized cap fast enough. Our per-unit sizing (§5.1) should be
  checked against an analogous curve per unit rather than assumed monotonic.
- **Placement dominates.** Active decap effectiveness fell from 133→74 mV
  average noise as `R_dist` (distance-proxy resistance) to the hot spot
  dropped from 10Ω→0.1Ω — a far bigger swing than the *passive* decap showed
  (134→124 mV) over the same range. Directly supports §5.1's finding: a
  pump has to sit physically close to the unit it's serving (CSR's own
  local rail, not a shared/lumped one) or most of the benefit is lost.
- **10–20% is the real, measured, silicon bound** for a comparator-triggered
  active decap of the SAME mechanism family (droop-triggered switch) as our
  release path. Our reported 3–25× *equivalent-decap* multipliers (§5.2) are
  a different metric (nF-equivalent, not %-noise-reduction) and a different
  charge mechanism (added charge vs. redistributed charge), so they are not
  directly comparable numbers — but the qualitative shape (real, worthwhile,
  strongly placement- and sizing-dependent, with a genuine sweet spot) matches
  published silicon, which is a meaningful sanity check on the whole approach.

## 7. Caveats (carried forward, still true)

- Current model is approximate (±30–40% absolute); ratios are what's
  load-bearing, not the mV numbers.
- ROB-capacity gate is optimistic (frees at issue, not commit); IQ capacity
  and memory deps aren't modeled; register-file aliasing adds some
  conservative false dependencies.
- TD is a naive hold policy (holds the costliest op, not criticality-aware)
  — reported costs are an upper bound for this policy class.
- Bank/pump timing in the Pulsim circuit is real (comparator delay via
  step_observer lag, real capacitor sag, real switch R_on) but the die-side
  PDN is still a single lumped series RLC per unit, not a floorplan-aware
  distributed model — consistent with Meng & Saleh's point that placement
  matters a lot, which we can't yet quantify for our own layout.
- Per-unit results are on dhrystone only; MUL/DIV's "zero benefit" finding
  is workload-specific and should be re-checked on a MUL/DIV-heavy kernel.

## 8. Reproduce

```
cd pdn_work
python3 sweep.py mm qsort dhrystone          # §2 reorder sweep
python3 perunit.py dhrystone                  # §5.1 per-unit decomposition + conservation check
python3 chargepump_perunit.py dhrystone 20000 # §5 full per-unit charge-pump table + real-circuit spot-check
python3 pulsim_validate.py dhrystone 3000     # §5 PDN solver cross-validation (Pulsim vs closed-form)
```
RTL (chipyard, `PaperMatchBoomV3Config`): `sims/verilator/simulator-chipyard.harness-PaperMatchBoomV3Config-debug
+permissive +max-cycles=N +verbose +permissive-off <elf>` for `BACT_CAND`/`BACT_TOK`.
