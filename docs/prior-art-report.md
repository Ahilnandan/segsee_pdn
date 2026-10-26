# Prior-art report: what is covered and what is still open

Prepared 2026-10-04 for the SEGSEE power project: issue-level current shaping plus an active decap bank plus FIVR.

**Method.** For each patent below I read the claims and description in full. The papers I read as PDFs. This is an engineering prior-art map, not a legal freedom-to-operate opinion. Claims are paraphrased from the documents themselves; quotes are marked.

## 1. Sources read

### Patents

| # | Patent | Assignee, priority | Core claim (from the document) |
|---|---|---|---|
| P1 | **US7937563B2** | AMD, 2008 | Scheduler keeps an **activity weight table** (per-instruction estimated node-capacitance switching). It compares recent vs old moving-average activity against thresholds and **throttles issue width** (issues fewer instructions) when exceeded. No reordering, no decap or VR interaction, no falling-edge or dummy activity. |
| P2 | **US11204766B2** | Intel, 2019 | Decoder/scheduler identifies **high-power instructions** (e.g. matrix multiply) and asserts a pre-charge signal to the VR. After a programmable start delay it asserts **Vboost**: the digital FIVR's PWM setpoint is switched up, with a stop delay afterwards. Issue is **never held or reordered**, no capacitor storage, no overshoot handling. |
| P3 | **US11625087B2** (and US11860707B2, same family) | IBM (now Blue Heron), 2021 | Per-core **power proxy → charge value → charge accumulators**, plus a **global "charge tank" counter that models VRM delivery** (a recharge value subtracted every cycle). When charge requests exceed delivery capacity it **throttles fetch rate / in-flight instructions** (64-level FTX) and delays grants. The tank is a **modeled counter, not a physical capacitor**; it gets no VR state feedback. |
| P4 | **US11966786B2** | IBM, 2021 | Proxy-based energy targets per period with a hardware **PI loop** setting the throttle index. Energy and power management, not di/dt. |
| P5 | **US7607028B2** | IBM, 2006 | **Senses actual voltage per execution unit**; throttles scheduling to that unit below a threshold. It **schedules "dummy executions" above a higher voltage threshold** and claims reduced decoupling (~500 nF → 50 nF). It does not control or interact with the VR. |
| P6 | **US11960338B2** | Microsoft, 2021 | Hierarchical **activity smoothener** (tile/cluster/node) enforcing a di/dt limit **on both rising and falling edges**. Credits plus **"fake work"** on falling edges, floor and ceiling on current. No voltage sensing, no VR or decap interaction. |
| P7 | **US11586265B2**, **US11275644B2**, **US11029742B2** | IBM, 2017–2021 | Predict droop from early pipeline events, stall events or critical-path-monitor timing margin, then throttle. No decap, VR loop or dummy activity. |
| P8 | **US7421604B1** | NVIDIA, 2005 | GPU pipeline-activity → **feed-forward into the VR feedback error signal**. Notes it allows smaller or higher-ESR output caps. No issue holds; switching frequency not addressed. |
| P9 | **US10243456B2** | NXP, 2017 | LDO gate boost triggered by **instruction-execution monitoring** (ALU or memory ops) that predicts a load increase. |
| P10 | **US7957847B2** | WD/HGST, 2005 | VR feed-forward from **deterministic load gate signals** (HDD controller), with weighted sum. No processor-side load modulation; no claim on reducing VR bandwidth, capacitance or frequency. |
| P11 | **US10972083B2** | IBM, 2019 | **Second decap pre-charged to N×Vdd** by a booster or charge pump. It connects to the rail on **detected droop *or* a predicted high-power operation**, notified by the functional block. Recharged from the booster; **overshoot not addressed**; no scheduling. |
| P12 | **US8467261B2** | IBM, 2010 | **Smart switched decap** (series/parallel) **triggered by logic's own high-current event signal** (bitline precharge, wake-up), recharged through a weak resistor. Unidirectional. |
| P13 | **US11811303B2** | Apple, 2021 | PIN-diode **stored-charge reverse recovery** decap. It charges during low activity (from the supply) and sources a boost current on voltage drop. **Not charged from overshoot.** |
| P14 | **US7839129B2** | Univ. Michigan, 2007 | Capacitor charge injection on undershoot; **overshoot handled by a separate load transistor that dissipates.** |
| P15 | **US12424846B2** | NVIDIA, 2023 | **Dummy load circuits** (digitally sized) that activate on overshoot above a threshold. **Dissipative.** No scheduling interaction. |
| P16 | **US7893671B2** | TI, 2007 | LDO with current source plus **current sink** for overshoot. Dissipative. |

### Papers

| # | Paper | Covers |
|---|---|---|
| A1 | **Grochowski, Ayers, Tiwari, HPCA 2002** (Intel) | On-die current estimate convolved with the PDN impulse response → **two-threshold feedback**. Below the lower bound: clock-gate or stop issue. Above the upper bound: **force idle units on**. "Cut … variation in half with less than 2% performance loss." Future work names an IIR filter for the ~1 MHz VR LC ring. |
| A2 | **Joseph, Brooks, Martonosi, HPCA 2003** | Control-theory threshold controller: voltage-low → disable units; voltage-high → **"phantom firings"**. Analyses sensor delay. Positioned as letting packaging run at **200% of target impedance**. |
| A3 | **Powell & Vijaykumar, ISCA 2003, "Pipeline damping"** | **Issue-stage current allocations with 4-bit per-op integers**, and a history register of future-cycle allocations (≈ our `need()` projection). **Upward damping = delay issue; downward damping = issue extraneous integer ALU ops**: 33% worst-case di/dt cut at 7% slowdown. |
| A4 | **Powell & Vijaykumar, ISCA 2004, "Resonance tuning"** | Senses chip current. Reacts only to *repeated* variations in the resonance band and shifts current-variation frequency out of the band (stall issue widths). |
| A5 | **Powell & Vijaykumar, 2005, "Balancing resource utilization…"** | Spatial balancing across replicated ALUs/IQ/RF to cut *thermal* power density. |
| A6 | **Zou et al., MICRO 2018 (and TCAD 2020), voltage-stacked GPUs** | Distributed **charge-recycling IVR plus architecture control under control theory**: dynamic issue width scaling, **fake instruction injection**, dynamic current compensation, proportional state feedback. **Architecture handles low/mid-frequency noise, the IVR handles high frequency**, which cuts IVR area 88%. The problem is inter-layer imbalance in a voltage stack. |
| A7 | **Lefurgy et al., MICRO 2011** (IBM POWER7) | Critical-path monitors drive clock frequency within cycles, plus a throttling loop; >50% fewer worst-case droop events. |
| A8 | **Meng & Saleh, JSSC 2009** | Comparator-triggered series/parallel active decap; 10–20% measured; sizing sweet spot; placement-sensitive. |
| A9 | **Gu, Eom, Kim, VLSI 2007**; **Gu & Kim, TVLSI 2009** (= P14 lineage) | Switched decap for resonance damping (11× effective decap); op-amp active decap. **Overshoot handled by a dissipative artificial-load transistor.** |

## 2. Coverage matrix: our mechanisms vs prior art

| # | What we built and measured | Covered by | Status |
|---|---|---|---|
| E1 | Per-uop 4-bit current tokens at issue (amp/dur ROM, `BACT_CAND`) | A3 (4-bit current allocations), P1 (weight table) | **Covered** |
| E2 | Hold issue when the projected charge need exceeds the budget (TD) | A3 (upward damping), P3 (charge tank + VRM recharge model → throttle), P1 | **Covered** |
| E3 | Dummy fill on falling edges using free ALU ports | A3 (extraneous ALU ops), A1, A2, P5, P6, A6 (fake-instruction injection) | **Covered** |
| E4 | Pre-charged active bank, triggered on droop or predicted from issue | P11 (N×Vdd decap, droop or prediction), P12 (logic-event-triggered switched decap), P2 (issue-time boost, VR setpoint) | **Covered** |
| E5 | **Scheduler tokens = the *sensed* charge of a physical active reservoir** (closed loop scheduler↔bank) | P3 uses a *modeled* tank of package capacitance with a constant VRM rate; no source gates issue on the measured state of an active charge reservoir | **Open (narrow)** |
| E6 | **Bidirectional bank that *stores* overshoot into its reservoir** (edge harvesting) and spends it on the next droop | P14/A9, P15, P16 dissipate overshoot; P13 charges from the supply in quiet periods; P11 recharges from a booster | **Open** |
| E7 | CPU regulates its own supply: hold on deficit, phantom/dummy on surplus | A1, A2, P5, A6 | **Covered** |
| E8 | Architecture control lets the regulator be cheaper | A6 (IVR area −88%), P5 (decap ÷10), A2 (200% impedance), P8 (smaller output caps) | **Covered in principle** |
| E8a | **The specific split: FIVR keeps only a slow DC integrator; issue does the fast (ns–100 MHz) regulation; the FIVR's switching frequency is cut 2–4×** | A6 uses the *opposite* split (architecture low/mid frequency, IVR high frequency) and targets area, not fsw; none removes the fast compensator of a single-rail buck FIVR | **Open (narrow)** |
| E9 | **Issue decisions driven by the FIVR's sensed inductor-current surplus** (VR state → scheduler) | P8–P10 and P2 run processor → VR only; A1/A2/P5 sense voltage or estimate current; A4 senses chip current | **Open** |
| E10 | Pure reordering at equal issue count | Not claimed anywhere. Our measured result is that it does *not* help, which is a negative result, not a claim. | n/a (publishable insight) |
| E11 | Per-execution-unit throttling, bank placement | P5 (per-unit sensing/throttle), P11 (per functional block), A8 | **Covered** |
| E12 | Spatial ALU-port steering | A5 (thermal version) | **Covered** (and only ≤3–4% droop benefit for us) |
| E13 | Issue → VR feed-forward | P8, P9, P10, P2 | **Covered** |
| E14 | Spread-spectrum / resonance-frequency shifting | A4 | **Covered** (and it failed for us) |

## 3. Loose ends we can attack

1. **L1 – Storing vs dissipating overshoot (E6).** Every overshoot mechanism found (P14–P16, A9) dissipates the charge. Our bidirectional bank stores it. A head-to-head at equal area and identical trigger would show whether storing gives better droop and better energy than the prior-art dissipative sink. If it does, that is the strongest open claim.
2. **L2 – VR-state signal vs voltage sensing (E9).** The prior art senses voltage (A1, A2, P5), current (A4) or a model (A3, P3). Test the same actuator driven by (a) voltage only (≈ prior art), (b) inductor-current surplus only, (c) both.
3. **L3 – Sensed vs modeled tokens (E5).** IBM P3's tank is a model. Compare a modeled charge-tank gate (P3-style, swept conservativeness) against sensed-reservoir tokens on a droop-vs-slowdown Pareto.
4. **L4 – Fast-loop removal and fsw reduction (E8a).** Two open questions:
   - Does the energy balance favour it? FIVR switching-loss saving vs dummy or harvest energy.
   - Is the combined loop stable? There is no formal analysis yet.

   Quantifying the FIVR loss trade makes the claim concrete.

Results of attacking these are appended in §4.

## 4. Attacking the loose ends: results

Each gap was tested head-to-head against its prior-art counterpart in the closed-loop full-system model:
- the scheduler is `replay2.run` turned into a generator, stepped inside the switching FIVR/off-chip chain;
- the bank replaces 0.5 nF of the 5 nF passive decap, keeping total capacitance equal;
- workloads are qsort, dhrystone and mm as full traces;
- the FIVR runs at 140 MHz and 70 MHz.

Code: `gaps_exp.py`; data in `gaps.jsonl` and `gaps_hybrid.jsonl`. Current model ±30–40% absolute, so compare variants against each other.

### L1 – Overshoot handling: store (ours) vs dissipate (prior art: P14–P16, A9)

Same trigger (V_die > V_nom + 10 mV) and same sink magnitude in all three variants.

| | prior-art dissipative sink | storing bank only | **hybrid: store until full, then dissipate** |
|---|---|---|---|
| p99.9 droop | best | 7–36% worse than dissipative | within 1–8% of dissipative |
| overshoot | best | much worse (reservoir fills, then stops absorbing) | ≈ dissipative (±5 mV) |
| wasted energy, % of load energy | 2.1–11.3% | 0.2–0.7% | **1.1–7.6%: 35–60% less than dissipative** |

Example cells:

| Workload / FIVR / policy | Dissipative: p99.9, overshoot, energy | Hybrid: p99.9, overshoot, energy |
|---|---|---|
| mm, 140 MHz, no scheduling | 82.8 / 112.9 mV / 7.19% | 84.7 / 109.7 mV / **2.84%** |
| dhrystone, 140 MHz, holds | 31.2 / 63.4 mV / 2.06% | 32.4 / 63.9 mV / **1.06%** |

**Verdict.**
- A storing-only bank is *not* a better overshoot suppressor than prior-art dissipative sinks.
- The hybrid reaches near-equal noise at roughly half the wasted energy. It is the most defensible open element I found: no source read in §1 stores overshoot charge for reuse.

### L2 – Issue signal: voltage only (≈ A1/A2/P5) vs FIVR inductor-current surplus vs both

Same actuator (holds only), same storing bank, slow-integrator FIVR.

| | at 140 MHz | at 70 MHz |
|---|---|---|
| inductor-current term (vs voltage only) | p99.9 6–13% better (qsort 52.2→46.3, mm 46.8→40.5, dhrystone 36.3→32.7), but slowdown about +1 percentage point (1.9→3.1%) | no better; voltage-only is equal or better (qsort 93.4 vs 94.9–107.1) |

**Verdict.** No clean win. Part of the gain is simply more holding. The inductor-current signal (E9) is not established as better than voltage sensing.

### L3 – Tokens: sensed reservoir charge (ours) vs modeled charge tank (≈ IBM P3)

TD hold, FIVR at 140 MHz. The modeled tank's capacity was swept from 250 down to 30 pC to trace a slowdown/droop Pareto.

| | sensed | best modeled tank (Qb = 30) |
|---|---|---|
| qsort | 1.75% slowdown, p99.9 54.7 mV | 1.56%, 60.6 mV |
| dhrystone | 1.79%, 43.7 | 2.09%, 45.6 |
| mm | 4.02%, 49.7 | 4.35%, 51.5 |

**Verdict.**
- A well-tuned modeled tank gets most of the way.
- Sensed tokens are modestly better on the Pareto: 3–10% lower p99.9 at similar or lower slowdown on 2 of 3 workloads, about equal on the third.
- The untuned model (Qb = 250, the value from the open-loop replay) is far worse: p99.9 is 30–45% higher. So the practical value of sensing is that it **needs no calibration**, not a large margin over a carefully calibrated model.

### L4 – Removing the FIVR's fast compensator and lowering fsw (E8a): weakened

**What Intel's FIVR paper (Burton et al., APEC 2014) shows:**
- 140 MHz switching, chosen so the output filter fits on package (air-core inductors, 7 nF/A on-die MIM).
- A programmable **type-3 compensator with 78–80 MHz unity-gain bandwidth and 40° phase margin**, plus nonlinear duty saturation on large steps.
- ~90% efficiency.

**Why our result does not transfer to it:**
- I tried to build an equivalent high-bandwidth baseline on our model. With realistic interleaved-PWM delay (≥1.8 ns) and our per-core output filter (f₀ ≈ 100 MHz, Q ≈ 10), the only robust single-loop designs have sub-MHz bandwidth.
- A ~110 MHz design that looks fine on the averaged plant (0.9 ns delay) is **unstable in the switching simulation**: 305 mV step droop vs 67 mV.
- So every "normal FIVR" baseline in the earlier reports was the best compensator achievable *on our modeled plant*, not a real 78 MHz FIVR.

**Verdict.** The claimed 2–4× fsw reduction is measured against a weak baseline. It may shrink substantially against a real FIVR. Closing this needs a faithful FIVR plant: 16 phases, double-edge modulation, Intel's filter, a type-3 design.

**Energy balance (not yet favourable).** Lowering fsw 140 → 35 MHz saves at most the switching part of the ~10% FIVR loss, i.e. a few % of delivered power. The dummy/shunt energy that very-low-fsw operation needs is 9–42% of core energy. Only the holds-only design point (zero dummy energy) at 140 MHz is energy-positive.

## 5. Bottom line

| Element | Status after reading the sources and testing |
|---|---|
| Issue tokens, holds, dummy fill, CPU-as-regulator, issue-triggered boost, VR feed-forward, per-unit throttling, spatial steering, resonance shifting | **Prior art** (A1–A6, P1–P13) |
| **Hybrid store-then-dissipate overshoot bank** | **Open and better than prior art on energy**: equal noise at roughly half the wasted energy vs dissipative sinks. Best candidate. |
| Scheduler gated on *sensed* active-reservoir charge | Open; modest benefit. Mainly removes calibration; small margin over a tuned model. |
| Inductor-current surplus as issue signal | Open; benefit not established (mixed vs voltage sensing). |
| FIVR fast-loop removal / fsw reduction | Open in form, but **evidence weak**: our baseline was far below a real 78 MHz FIVR; the energy balance is unfavourable except holds-only at 140 MHz. |
| Pure reordering | Not claimed anywhere; **measured not to help**. Publishable as a negative result with the resonance-band explanation (mm puts 17.8% of current spectral energy in the PDN band; holds cut it to 3.7%, reorder cannot move it). |

**Suggested next steps:**
1. Design the hybrid bank at circuit level: a boost path into the reservoir plus a dissipative fallback, with real switch losses.
2. Build a faithful FIVR plant model with a 78 MHz type-3 loop and re-test the fast-loop-removal claim.
3. Get a formal patent search on the hybrid-bank and sensed-token claims.
