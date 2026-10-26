# Dummy instructions, a novelty search, and real power-supply chains

Status 2026-10-03. Machine `PaperMatchBoomV3Config`; workloads mm / qsort / dhrystone (full traces, HTIF load phase cut).
Everything below runs on the validated replay (`replay2.run`, regression-checked: every new option
disabled reproduces `sweep.jsonl` exactly) and on Pulsim circuit models validated against closed-form
results before use. Current model is still the approximate `curmodel.py` (±30–40% absolute): compare
policies with each other, don't quote absolute mV as silicon.

## 1. Answers in one place

| Question | Answer |
|---|---|
| Do dummy instructions help? | **Yes, but they cost energy.** They cut droop and slowdown, at +3% to +155% dynamic energy depending on the threshold. Fill threshold THF=200 mA is the cheap point: −8..−16% worst droop for +3..4% energy. "Ramp" dummies make things worse. |
| Combined (dummies + token hold + bank)? | Best overall. On the real FIVR chain it lets the FIVR switch at **~59–77 MHz instead of 140 MHz** at the same droop (combined200, +3–4% energy). |
| Off-chip high-current buck vs off-chip multiphase → FIVR? | The FIVR holds DC exactly (0.900 V, no load-line sag). Its switching ripple makes *fast* droop worse on raw current. Once current is shaped, both land in the same 55–90 mV range. |
| Something novel? | Strongest candidate is an **edge-harvesting bank**, but at circuit level it only gives a modest p99.9 gain. The firmer result is a requirement: **the token counter must track the sensed reservoir charge.** An open-loop token model is ~2× optimistic on slowdown. Issue→regulator feed-forward turned out to be prior art. |

## 2. Dummy instructions

Implementation (`replay2.run(dummy=...)`): a dummy is an ALU-class op (shape 10/20/20/10 mA) placed on an
integer port the real picks left free. There is no performance cost, because only free ports are used.

| mode | rule | result |
|---|---|---|
| `fill` | falling edge: keep I[t+1] ≥ (max of last 5 cycles) − THF | works; energy set by THF |
| `ramp` | while the token gate holds real ops, ramp dummies ≤1/cycle | **fails**: the dummies' own current triggers events, so slowdown goes *up* (mm 1.9%→4.7%) |
| `bridge` | fill only if real work is pending or waking within 4 cycles | saves energy only on qsort; mm/dhrystone almost always have work pending |

Fill-threshold Pareto, token hold + bank 250 pC/25 mA, TH=100 mA (cells: slowdown / worst droop / extra dynamic energy):

| | no dummies | THF 100 | THF 150 | **THF 200** |
|---|---|---|---|---|
| qsort | 0.49% / 100.6 mV | 0.05% / 91.1 / +86% | 0.24% / 72.8 / +17% | 0.38% / 84.1 / **+4%** |
| dhrystone | 1.03% / 85.7 | 0.00% / 71.2 / +56% | 0.03% / 67.6 / +29% | 1.02% / 75.6 / **+3%** |
| mm | 1.93% / 106.0 | 0.55% / 111.4 / +97% | 1.06% / 101.0 / +29% | 1.78% / 97.8 / **+4%** |

Without any holds (policy A + bank), fill alone cuts qsort droop 165→99 mV at 0% slowdown, but at +87% energy.

Dummies help for two reasons:
- **They smooth falling edges.** This is the slowdown benefit.
- **They raise the recent-minimum baseline** that later rising edges are measured against. So fewer events need the bank or a hold.

Prior art: Grochowski et al., HPCA 2002, "phantom firing". What's new here is how dummies interact with the bank and the token gate, and the energy Pareto above.

## 3. Novelty search

**Prior art found** (WebSearch only; patent pages were not opened):

| Idea | Status |
|---|---|
| Issue-time detection of a high-power instruction → pre-charge → voltage boost | Prior art: Intel US11204766 |
| Activity-based issue throttling | Prior art: AMD US7937563 |
| Pipeline/issue → voltage-regulator load feed-forward | Prior art: US7421604, US10243456, US9083229, US7957847; load-current feed-forward literature (Duke) |
| Active decaps with reverse charge recovery, boosted second decap, smart switched decaps | Prior art: US11811303 / US12176803, US10972083, US8467261 |
| Phantom firing; resonance-aware throttling | Prior art: Grochowski HPCA'02, Joseph HPCA'03, Powell & Vijaykumar |

**Candidate I couldn't find in prior art: the edge-harvesting bank** (`replay2.run(recyc=η)`, `pulsim_harvest.py`).
- On a falling edge the scheduler already knows about, the bank sinks grid current. That smooths the edge the same way a dummy does.
- Instead of burning the current, it pumps it into the reservoir (η=0.8–0.85). The token gate then spends it on the next rising edge.
- A search isn't a clearance, so this is a candidate, not a claim.

How it held up, step by step:

| Level | Result |
|---|---|
| Replay, idealized bank | Better slowdown **and** droop than matched controls (same bank, same recharge) in 26/27 cells, at ≤0.8% wasted energy. E.g. 500 pC/10 mA: mm 3.45%/102.9→2.63%/92.0 mV; qsort 1.17/96.3→0.64/85.3; dhrystone 1.46/88.4→0.89/72.9. With no external recharge rail at all, the bank lives on harvested charge. |
| Real bank circuit, fixed schedule | Worst droop lower in **6/6** cells (qsort 182→147 mV, mm 166→149, dhrystone 132→118 at Ir=10). |
| Real circuit + token hold, open-loop tokens | Mixed, and at Ir=10 it **reverses on all three workloads**. |
| Closed-loop co-simulation (`cosim.py`) | The scheduler is `replay2.run` itself, turned into a generator stepped from inside the circuit simulator. Tokens = sensed reservoir charge. Regression-exact against the open-loop circuit. Worst droop ≈ unchanged; **p99.9 droop 2–9% better** (qsort 71.6→64.9 mV, mm 73.3→69.8); slowdown unchanged. A real but modest effect. |

**The firmer finding: token-counter fidelity.** At circuit level, an open-loop token model overestimates the deliverable charge:
- The 0.5 nF reservoir holds ~135 pC usable, not the 250 pC assumed.
- It only releases after the comparator trips.

Results with sensed tokens instead, in the real circuit:

| Workload | Open-loop tokens | Sensed tokens |
|---|---|---|
| qsort | 0.49% / 110 mV | 2.03% / 85 mV |
| dhrystone | 1.03% / 104 mV | 1.95% / 78 mV |
| mm | 1.93% / 125 mV | 4.52% / 100 mV |

So worst droop is 20–25% lower and slowdown roughly 2× higher. **Every earlier replay slowdown figure (including `bank-token-scheduling-result`) should be read as ~2× optimistic.** The two cheapest ways to win that slowdown back, both measured in closed loop:
- **A 1 nF bank:** mm droop 100→79 mV at similar slowdown.
- **Fill dummies at THF=200:** slowdown 2.03→1.71% (qsort), 4.52→4.06% (mm).

## 4. Real power-supply chains (`pdnchain.py`, Pulsim, switching models)

| Case | Path | Control |
|---|---|---|
| **A: off-chip high-current buck, direct** | 4-phase synchronous buck, 12 V→0.9 V, 500 kHz/phase → bulk cap → package (507 pH / 79.6 mΩ) → 5 nF die → load | average-current-mode, 100 kHz inner / 30 kHz outer loop |
| **B: off-chip multiphase → FIVR** | 4-phase buck 12 V→1.8 V, 500 kHz → package → 10 nF FIVR input decap → on-die 4-phase interleaved FIVR 1.8→0.9 V, fsw swept 35–200 MHz (L = 2 nH·140/fsw, so ripple stays constant) → 20 mΩ grid → 5 nF die | off-chip VR as above; FIVR voltage-mode with nominal-duty feed-forward |

Off-chip components are a consistent per-core share of a real 100 A VR: L, DCR and ESR ×250, C ÷250. That keeps f₀, Q and ripple ratio real.

The full traces (245k–600k cycles) go through every converter with no windowing. Droop is measured per cycle on the sub-cycle minimum:
- *droop* = below 0.9 V absolute.
- *fast* = below a 5 µs moving average, which separates fast droop from slow regulation drift.

Controller lessons found by bisection (each was a real bug before it was fixed):

| Controller | What went wrong | Fix |
|---|---|---|
| FIVR, plain voltage-mode PI | Oscillated at ~55 MHz against the output LC (f₀ 71 MHz, Q≈14) | — |
| FIVR, current-mode with instantaneous current sampling | Unstable at D≈0.5 (subharmonic) | — |
| FIVR, average-current-mode | Stable, but with only 5 nF it turns the FIVR into a current source: 570 mV droop on a 200 mA step | — |
| FIVR, final | — | Voltage-mode with nominal-duty feed-forward, small Kp, slow integrator. Keeps the stage's natural LC stiffness. |
| Off-chip VR | Inconsistent scaling gave Q≈29 and a slow instability; a 1 kHz integral loop drifted 160 mV | Consistent per-core-share scaling + average-current-mode |

Results (droop mV; FIVR at 140 MHz, 5 nF). The reference is raw current, FIVR 140 MHz, 5 nF. iso-fsw = lowest FIVR frequency that matches the reference; iso-C = least die decap that matches the topology's own raw 5 nF droop.

| | slowdown | extra energy | off-chip (fast) | FIVR 140 | **iso-fsw** | iso-C off-chip / FIVR |
|---|---|---|---|---|---|---|
| mm raw | 0 | 0 | 124.7 | 155.7 | 140 MHz | 5 / 5 nF |
| mm hold+bank | 1.93% | 0 | 83.5 | 84.2 | 103 MHz | 2.5 / 3.5 |
| mm combined200 | 1.78% | 3.8% | 85.0 | 83.4 | **59 MHz** | 2.5 / 3.5 |
| qsort raw | 0 | 0 | 142.9 | 142.2 | 140 MHz | 5 / 5 |
| qsort hold+bank | 0.49% | 0 | 72.0 | 87.9 | 101 MHz | 2.5 / 3.5 |
| qsort combined200 | 0.38% | 3.3% | 71.2 | 88.8 | **77 MHz** | 2.5 / 3.5 |
| dhrystone raw | 0 | 0 | 87.8 | 122.8 | 140 MHz | 5 / 5 |
| dhrystone hold+bank | 1.03% | 0 | 55.6 | 77.5 | 79 MHz | 2.5 / 3.5 |
| dhrystone combined200 | 1.02% | 3.0% | 54.5 | 64.9 | **72 MHz** | 2.5 / 3.5 |

The extra-energy column here is relative to total current including the floor.

What the table shows:
- **Project deliverables at iso-droop:**
  - **FIVR switching frequency:** 140 → 79–103 MHz with hold+bank, and → 59–77 MHz with the combined policy.
  - **On-die decap:** 5 → 2.5 nF behind the off-chip VR, and 5 → 3.5 nF behind the FIVR.
  - **Overshoot:** dummies also cut FIVR overshoot from 140–177 mV to 75–107 mV.
- **FIVR vs off-chip:**
  - The FIVR regulates DC perfectly. The off-chip rail sags with load and its loop takes ~20 µs to recover.
  - On raw current the FIVR's ripple (~35 mV peak-to-peak) makes fast droop worse.
  - With shaped current the two are comparable.
- **Issue→FIVR feed-forward is inconsistent:**
  - It helps mm under every policy (raw 155.7→145.0, hold+bank 84.2→75.1).
  - It hurts dhrystone raw (122.8→142.9).
  - It does nothing at 35 MHz, because a 3-cycle prediction lead is useless when each phase holds its duty for 28 ns.
  - It is prior art anyway.
- **Same open-loop caveat as §3:** these chain traces come from the open-loop token replay, so the slowdowns listed are ~2× optimistic. The co-sim shows droop holds up or improves when tokens are sensed.

## 5. Caveats
- The current model is approximate. The FIVR and VR are per-core-share models, not a floorplanned multi-core rail.
- Bank and harvest efficiency (η) are assumed, not designed. The harvest path is behavioral: a boost transfer, not a switch-level converter.
- Worst-case droop is a single-event statistic, so it's noisy. p99.9 is reported where it matters.
- Prior-art checks are from search snippets only.

## 6. Reproduce (`pdn_work/`)
```
python3 sweep_dummy.py; python3 sweep_bridge.py; python3 sweep_thf.py      # dummies
python3 sweep_recyc.py; python3 sweep_recyc2.py; python3 sweep_recyc_ctrl.py  # edge-harvesting (replay)
python3 harvest_circuit_exp.py                                              # harvesting, real bank circuit
python3 cosim.py; python3 cosim2.py                                         # closed-loop scheduler<->circuit
python3 chain_exp.py; python3 analyze_chain.py chain_exp.jsonl droop|fast    # power-supply chains
```
