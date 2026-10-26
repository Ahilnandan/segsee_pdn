# Best fix for core supply droop (real units)

Status 2026-10-07. Builds on `current-stft-analysis.md` (droop = 50–200 MHz, mostly loop tones on the ~100 MHz resonance) and `parallel-buck-derivation.md` (real-unit network).

**Scripts and data:**
- `passive_fix.py` → `passive_fix.json`
- `active_fix.py`, `active_fix2.py` → `active_fix2.json`
- `active_fix3.py`, `injector_sw.py`, `injector_sweep.py`, `injector_gate.py` → `injector_*.json`

**Method:** the passive network is linear, so droop on the real per-cycle current traces is exact by superposition. The injector is simulated at switching level (50 ps steps, hysteretic current control) on that linear network.

## 1. Passive: the damped decap (biggest, cheapest step)

A plain capacitor shifts the resonance. A capacitor with a deliberate series resistor (R ≈ √(L/C)) absorbs the resonance energy, so loop tones can't build it up.

| Option (area at an assumed 12 nF/mm² MIM) | p99.9 droop qsort / dhrystone / mm |
|---|---|
| Baseline: 5 nF, 507 pH package | 108 / 89 / 132 mV |
| Half of the existing 5 nF as a damped branch (0.45 Ω), **zero extra area** | 87 / 76 / 90 |
| +5 nF plain | 80 / 63 / 70 |
| **+5 nF damped (0.2 Ω), 0.42 mm²** | **66 / 55 / 61** |
| +10 nF damped + 2× power bumps (250 pH), 0.83 mm² | 36 / 31 / 34 |

## 2. Active, tier "+5 nF damped"

**Predictive feed-forward injector.** The issue stage knows each op's 4-bit token ~13 cycles ahead. The injector sources or sinks the negative of the 50–200 MHz part of that predicted current, filtered with zero phase inside the lead window.

**Hardware:** a 1.8 V synchronous buck slice with L = 2 nH air-core, **1 Ω switches** (conduction is negligible, so small switches minimise switching energy), hysteretic current control with a ±100 mA window, and a slew limit including the switch IR drop.

| Injector mode | Active time | qsort p99.9 / loss | dhrystone | mm |
|---|---|---|---|---|
| Always tracking | 100% | 30.5 mV / 12.5 mW (6.9%) | 33.7 / 7.8 mW (4.3%) | 34.0 / 10.9 mW (6.1%) |
| Gated: only when the predictor sees the top 10% of dips | 55–71% | 43.0 / 9.0 mW | 35.7 / 6.5 mW | 35.5 / 9.4 mW |
| Gated: top 2% | ~20% | 53.7 / 3.6 mW (2.0%) | **37.3 / 2.7 mW (1.5%)** | 48.6 / 2.9 mW (1.6%) |

Loss is a percentage of the 180 mW core.

**Design notes:**
- **Original switch sizing was wrong.** With 0.1 Ω switches and a ±20 mA window the injector chatters at 4–6 GHz and burns 70–108% of the core's power. Switch size and window are the critical design choices.
- **Gating is a dial between energy and droop.** Tonal dhrystone keeps most of its gain at 20% activity. Broadband qsort needs more.
- **Rough energy balance:** a ~25–35 mV smaller p99.9 margin lets VDD drop by about that much, roughly −6 to −8% core power, against 1.5–7% spent in the injector. So it's net-positive at moderate gating, and better still if the margin is used for frequency instead.

**Predictive holds.** On a damped network, holds within a 1% slowdown budget add little: qsort −5 to −16%, mm −2 to −7%. Dhrystone can't stay under 1% at any target, because its droop is a continuous tone. Holds are therefore optional, mainly for broadband code. The hold predictor's kernel has to be high-passed at 10 MHz; otherwise it chases the slow IR droop through the 80 mΩ package and throttles needlessly.

## 3. Recommended configuration

1. **Damped decap**, +5 nF with a 0.2 Ω series resistor (or the zero-area conversion if area is tight).
2. **Predictor-driven feed-forward injector** (2 nH, 1 Ω switches, ±100 mA hysteresis), with gating as the energy/droop dial.
3. **Optional light holds** (target ≈ 0.8–0.9 × p99.9) for broadband workloads.

**Result:** p99.9 drops from 108 / 89 / 132 mV to about **30 / 34 / 34 mV** with the injector always on (−62 to −74%, ~4–7% power, 0% slowdown), or about 54 / 37 / 49 mV with the top-2% gating (~1.5–2% power).

## Caveats (verify before claiming)
- **Not yet in Pulsim:** the injector is a switching-level model on a linear network. The next step is closed-loop Pulsim with real switch models, including the injector's effect on the 1.8 V input rail, which is assumed stiff here.
- **Assumed switch energy:** 3.25 pJ per transition at 1 Ω, scaled from the first-principles 65 pJ per pair at 0.1 Ω.
- **Assumed capacitor density:** 12 nF/mm² MIM; deep-trench would be ~10× denser.
- **Unverified package model:** the 507 pH / 79.6 mΩ package and the board VR loop approximation.
- **Sinking is not credited:** the injector also sinks current, and the energy it returns to the 1.8 V rail isn't counted.
- **Narrow workload set:** three microbenchmarks, unmodified issue.
