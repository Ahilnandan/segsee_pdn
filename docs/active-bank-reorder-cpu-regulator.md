# Active bank in the real supply, does reordering help, and the CPU as its own regulator

Status 2026-10-03. Every result here is **closed-loop**:
- The scheduler is `replay2.run`, turned into a generator (`cosim.run_gen2`) and stepped once per cycle from inside the switching circuit simulation (`pdnchain.simulate_chain(sched=...)`).
- The scheduler receives two signals each cycle:
  - the *sensed* bank charge;
  - a supply-state signal, SIG = sensed FIVR inductor current − load + (V_die − V_nom)/Z₀.
- Every new option was regression-checked before use:
  - With defaults, the chain model reproduces all 255 stored open-loop results exactly.
  - The closed-loop age-order run reproduces the open-loop run exactly.
  - Each new policy, given no signal, reproduces age order exactly.

Workloads are mm, qsort and dhrystone, as full traces on `PaperMatchBoomV3Config`. The current model is still approximate (±30–40% absolute), so compare policies against each other rather than quoting absolute millivolts.

**Measurement note.** On a switching FIVR, the single worst-cycle droop carries about ±6% "phase noise": rounding the input current to float32 instead of float64 flips one PWM edge and gives 142.2 mV instead of 150.85 mV on the same trace. The p99.9 droop is stable (104.6 vs 104.8 mV). So the tables below lead with p99.9.

## 1. The active bank in the real supply

The bank is a 0.5 nF reservoir precharged 0.25 V above V_DD:
- It releases into the rail through a comparator-gated switch.
- It recharges at a limited 25 mA from the FIVR's 1.8 V input rail.
- In the **bidirectional** variant it also absorbs overshoot back into the reservoir, through an energy-correct boost transfer with η = 0.85.

**Equal-capacitance replacement** means 4.5 nF passive + 0.5 nF active, versus 5 nF passive. All runs: FIVR at 140 MHz, normal issue, no scheduling change.

| p99.9 droop (mV) | 5 nF passive | 4.5 + 0.5 bank | 4.5 + 0.5 bidirectional bank |
|---|---|---|---|
| qsort | 104.8 | 88.8 | 88.4 |
| dhrystone | 101.5 | 80.9 | 66.7 |
| mm | 119.2 | 97.5 | 87.5 |

That is **−16 to −34% p99.9 at zero performance cost**. Overshoot also drops (dhrystone 139→88 mV).

Going to a larger bank (4 + 1 nF) adds nothing over 0.5 nF. This is a sweet spot, consistent with Meng & Saleh's sizing curve.

## 2. Does reordering help?

**Not reliably, across every formulation tried.**

| Reorder policy | Result |
|---|---|
| Bank-need-minimizing reorder (T), real supply + active bank | No change (qsort 150.8→150.5 mV) |
| Supply-following reorder (SF): pick heavy ops when the supply has surplus, light ops when in deficit | Mixed. Helps some cells (dhrystone with bank: p99.9 75.2→62.3) but hurts others. Costs 0.2–2.2% slowdown. |
| Spread-spectrum issue: random heavy-first/light-first, to decorrelate loops from the PDN resonance | Fails. Resonance-band energy unchanged (mm 17.8%→17.9% of current spectral energy). Worst-case droop is flat or worse. Up to 2.1% slowdown. |
| Spatial ALU-port steering (pure spatial freedom, zero performance cost) | ≤3–4% droop reduction (local droop is small next to the global resonance). Thermal-balancing prior art: Powell & Vijaykumar 2005. |

**Why reordering fails and holding works.** mm's loop puts 17.8% of its current-spectrum energy right in the PDN resonance band. Holding cuts that to 3.7%; reordering cannot move it. The freedom to swap which ready op goes this cycle only shifts current by a cycle or two, which is too small to change the spectrum.

## 3. The outrageous one: the CPU regulates its own supply

**Idea.** Remove the FIVR's fast compensator: either fully open-loop (fixed duty), or keep only a cheap slow integrator for DC. Then make **instruction issue the fast regulation loop**:
- When SIG says the supply is in deficit, hold the costliest picks, with a budget that shrinks as the deficit grows.
- When the supply has a surplus, shunt it, either by issuing dummy ALU ops or by absorbing it into the bidirectional bank.

FIVR fast loop removed, slow integrator only, bank replacing 0.5 nF of passive. Each cell is p99.9 droop / slowdown / extra dummy energy:

| | normal FIVR, passive 5 nF, normal issue | CPU-regulator (dummy shunt) | **holds + bidirectional bank** | **CPU-regulator + bidirectional bank** |
|---|---|---|---|---|
| qsort 140 MHz | 104.8 mV | 72.8 / 2.3% / +42% | **46.3 / 3.1% / 0** | **36.9 / 2.3% / +29%** |
| qsort 70 MHz | 214.6 | 69.5 / 1.8% / +39% | 94.9 / 2.7% / 0 | 54.8 / 2.0% / +28% |
| qsort 35 MHz | 338.2 | 94.2 / 1.4% / +41% | 188 / 2.2% / 0 | 83.2 / 1.8% / +27% |
| dhrystone 140 | 101.5 | 61.1 / 2.8% / +19% | **32.7 / 2.9% / 0** | **30.5 / 2.9% / +11%** |
| dhrystone 70 | 133.9 | 65.1 / 2.3% / +17% | 59.8 / 2.9% / 0 | 46.1 / 2.8% / +10% |
| dhrystone 35 | 208.7 | 68.3 / 2.0% / +15% | 102 / 2.4% / 0 | 60.2 / 2.8% / +9% |
| mm 140 | 119.2 | 69.3 / 5.2% / +42% | **40.5 / 6.5% / 0** | **34.8 / 5.2% / +29%** |
| mm 70 | 190.4 | 63.2 / 4.5% / +39% | 71.9 / 6.0% / 0 | 52.8 / 4.7% / +28% |
| mm 35 | 298.7 | 86.9 / 4.0% / +38% | 139 / 5.6% / 0 | 83.9 / 4.4% / +26% |

### What this shows
- **The CPU regulator plus the bidirectional bank at 140 MHz cuts p99.9 droop by 65–71%** and worst-case droop by 50–57% (e.g. qsort 150.8→66.6 mV, mm 155.7→67.1 mV).
- **At 35 MHz it still beats a conventional FIVR at 140 MHz** on all three workloads. That is a 4× lower switching frequency with a dumber converter.
- **DC is exactly 0.9000 V** when the FIVR keeps a slow integrator. Fully open-loop, the rail sits 6–11 mV low.

### Costs
- **Holds-only (no dummies), with the bidirectional bank:** zero dummy energy, −56 to −68% p99.9 at 140 MHz, 2.9–6.5% slowdown. It degrades at 70 MHz and fails at 35 MHz: a slow FIVR overshoots on falling edges, and holds can only fix deficits.
- **The surplus shunt costs energy:** +9% to +29% with the bidirectional bank, versus +15% to +42% when burned in dummies. Lowering FIVR switching frequency saves switching loss, but probably not that much. So the energy-efficient design point is **holds + bidirectional bank at 140 MHz**. The very-low-fsw point trades energy for converter simplicity.
- **Slowdowns are closed-loop:** sensed tokens, not the optimistic open-loop model.
- **The surplus signal needs the inductor-current term.** A voltage-only trigger is clearly worse (qsort 140 MHz: 93.5 vs 72.8 mV).

### Prior-art check (search only, not a clearance)
- **Already exists:**
  - Dummy loads triggered by an overshoot threshold, e.g. on-die current sink US12424846.
  - Droop-triggered throttling, e.g. US7607028, US7599808.
  - Issue throttling (US7937563) and issue-triggered capacitor boost (US11204766).
- **Not found:** removing the converter's fast compensator and making instruction issue (holds plus a surplus shunt) the primary fast regulation loop, co-designed with a bidirectional active bank, so a simpler FIVR can run at a fraction of its switching frequency.

## 4. Caveats
- **Behavioral parts:** the bank's harvest/absorb path is a behavioral boost transfer with assumed η, not a switch-level converter.
- **The supply-state signal assumes the scheduler can see the FIVR's average inductor current every cycle.** On-die FIVRs do sense it; the latency is assumed to be one cycle.
- **Per-core-share supply, approximate current model.** Both as before.
- **Not yet checked:** the REG gains `reg_k` and the 60 mA offset were not tuned; one setting was used throughout.
- **Stability:** no formal analysis yet of the combined CPU+FIVR loop. The loop stayed stable on all 3 full traces × 3 frequencies, but that is evidence, not proof.

## 5. Reproduce (`pdn_work/`)
```
python3 fullsys_exp.py      # active bank configs x 6 closed-loop policies, FIVR + off-chip
python3 sweep_dith.py       # spread-spectrum issue
python3 spatial.py          # spatial ALU-port steering
python3 reg_exp.py; python3 reg_exp2.py; python3 reg_exp3.py   # CPU-as-regulator, taming, bidirectional bank
```
