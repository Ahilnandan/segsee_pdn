# Core current: time-frequency analysis, total and per region

Status 2026-10-07. This is the characterisation step before designing any regulator.
- Code: `stft_analysis.py` (bands, tiles, figures), `stft_units.py` (unit coherence, bursts), tone search (inline, saved to `stft_tones.json`).
- Data: `stft_analysis.json`, `stft_units.json`, `stft_tones.json`.
- Figures: `docs/figs/stft_*.png`, `tiles_*.png`, `tile_time_*.png`.

**Inputs:** per-cycle current of the unmodified BOOM core (policy A, `PaperMatchBoomV3Config` traces), split over a 4×4 floorplan by the port and unit each op issued on (`replay2.run(tiles=True)`). The current model is ±30–40% absolute. The impedance weighting |Z(f)| comes from the real-unit derivation in `parallel-buck-derivation.md`: direct path, 1.155 Ω peak at 92 MHz.

## 1. Which frequencies matter

| Workload | Current energy <10 / 10–50 / **50–200** / 200–500 MHz | Droop energy (|Z|²-weighted), same bands |
|---|---|---|
| qsort | 17 / 31 / **46** / 6% | 2 / 7 / **90** / 1% |
| dhrystone | 16 / 35 / **32** / 17% | 2 / 5 / **92** / 1% |
| mm | 27 / 24 / **41** / 8% | 2 / 2 / **95** / 0% |

**90–95% of droop comes from the 50–200 MHz band**, even though only a third to a half of the current energy is there. Anything below ~10 MHz is the board regulator's job and contributes ~2%.

## 2. The droop band is mostly loop tones (dhrystone, mm)

The strongest narrow spectral lines (Welch spectrum, ≥4× above the local floor), ranked by droop weight:

| Workload | Tones' share of current / droop energy | Strongest tones (MHz → loop period in cycles) |
|---|---|---|
| dhrystone | 76% / **93%** | 76.9 (13.0), 92.3 (10.8), 107.7 (9.3), 115.5 (8.7), 100.1 (10.0), 84.7 (11.8); 110–490× above floor |
| mm | 63% / **73%** | 99.4 (10.1, 114× above floor), 103.0, 101.1, 81.8, 90.8 |
| qsort | 1% / 0% | only weak high-frequency lines (250–400 MHz); its droop band is broadband |

Hot loops of about 9–13 cycles repeat at 77–115 MHz, right on the die-C / package-L resonance (~100 MHz, i.e. a 10-cycle ring). **For these workloads, droop is a periodic loop driving the supply at resonance**, the classic di/dt stressmark mechanism, occurring naturally. qsort's droop comes from data-dependent, broadband current.

This also explains, after the fact, why holds worked and reordering did not. A hold changes the loop's effective period and phase, which detunes the forcing. Swapping which ready op goes this cycle leaves the loop period unchanged.

## 3. Steady background vs bursts

With 64-cycle windows, the 50–200 MHz power is mostly steady. The p99 window is only 2.6–3.4× the median, and p99.9 is 2.9–4.8×.

| Workload | Events above p99.9 band power | Duration | Typical gap |
|---|---|---|---|
| qsort | 0.06 per µs | 16–32 ns (1–3 resonance periods) | 3.7 µs |
| dhrystone | 0.07 per µs | 16 ns | 0.43 µs |
| mm | 0.06 per µs | 16–32 ns | 0.8 µs |

So there are two regimes. There's a **continuous background** of resonance-band activity, which sets typical droop. On top of it are **rare, short spikes** (one to three ring periods, roughly one per 15 µs) that set the worst case.

The resonance-band power correlates only weakly with the windowed average current (−0.37, −0.03, +0.29). **Droop activity can't be predicted from average power**; it needs cycle-level information, i.e. the issue-stage tokens.

## 4. Where on the core it comes from (per region)

Share of droop energy by tile (`tiles_*.png`):

| Workload | Dominant regions |
|---|---|
| qsort | ALU0 44%, ALU1 15%, LSU0 13%, D$ 12%, ALU2 6% |
| dhrystone | ALU0 25%, ALU1 22%, ALU2 14%, LSU0 11%, LSU1 7%, D$ 6% |
| mm | LSU0 20%, LSU1 19%, D$ 38%, ALU0 10%, ALU1 5%, FPU0 3% |

- **The hot region depends on the workload:** the ALU cluster (port 0 above all, because BOOM's port selection fills the lowest free port first) for integer code, and the load/store unit plus D-cache for mm. A per-region scheme has to cover both; one fixed hotspot isn't enough.
- **Front end, CSR and L2 interface contribute almost nothing** to the droop band. CSR spikes are big but rare.

## 5. Do regions move together? (global vs local)

**Coherence gain** in the 50–200 MHz band = energy of the summed unit currents ÷ sum of the units' own energies. It's 1 when units are independent and equal to the number of units when they're fully in phase. It's computed over real issue units, which removes the fixed floorplan splits that tie LSU to D$ and ALU to the register file by construction.

| Workload | Active units | Coherence gain, 50–200 MHz | Most coherent pairs |
|---|---|---|---|
| qsort | 7 | 2.6 | ALU0–MEM0 0.45, ALU0–ALU1 0.30 |
| dhrystone | 7 | 1.5 | ALU0–ALU2 0.51, ALU0–MEM0 0.46 |
| mm | 9 | 2.1 | MEM0–MEM1 0.61, ALU2–ALU3 0.40 |

Units are **partly in phase**: more than independent, far from fully coherent. They drive a shared **global** resonance together, while each region also has its own fluctuation. The global die-C / package-L mode dominates droop (section 1), so global control is primary. Per-region action helps where one region dominates, such as ALU0 in qsort (44%) or LSU/D$ in mm (77%).

## 6. What this means for design

1. **Target band: 50–200 MHz.** Regulators slower than ~50 MHz (board VR, conventional FIVR loops) can't touch 90%+ of the droop energy. Die-side measures and issue-level control can.
2. **Much of the droop is tonal and periodic** (dhrystone, mm). That opens options:
   - **Detune:** move the supply resonance away from typical 9–13-cycle loop periods (e.g. more die capacitance moves it down; a parallel inductive path moves it up).
   - **Break the periodicity:** small, well-timed holds change the loop period (this is why predictive holds work).
   - **Cancel the tone:** a periodic disturbance is predictable many cycles ahead, so a predictor-driven injector could supply anti-phase current at the loop frequency.
3. **Broadband droop (qsort-like) and rare spikes** need the predictive, event-based response: hold or injector fired on predicted large steps (~0.06 events/µs above p99.9).
4. **Per region:** the ALU cluster (port 0 first) and the LSU/D-cache are the places for local decap, local injectors or local predictor action. The front end, CSR and L2 interface aren't.
5. **Average power is a poor proxy.** Detection has to use cycle-level token information.

## Caveats
- **Floorplan:** the 4×4 floorplan and the fixed current splits (memory op 50% LSU / 50% D$; ALU op 85% ALU / 15% register file) are my assumptions. Unit-level coherence (section 5) avoids the splits, but tile maps (section 4) depend on them.
- **Workloads:** three single-thread microbenchmarks, unmodified issue. Loop tones are a property of these kernels; SPEC-like code will be less tonal, and stressmarks more.
- **Current model:** ±30–40% absolute. Ratios and frequency locations are more reliable than absolute power.
