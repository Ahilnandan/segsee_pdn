# Parallel mini-buck path: frequency-domain derivation (real units)

Status 2026-10-07. Code: `pdn_work/pbuck_freq.py`. Numbers: `pbuck_freq.json`. Figures: `docs/figs/`.
Everything here is sized for **one BOOM core at its real current**. There is no per-core-share scaling of on-die parts.

## 0. System and parameters

Two supplies feed the same die node:

```
Supply 1, direct:  VR1 (board, 3.3 V -> 0.9 V) -> board trace -> package decap -> package-to-die -> DIE (C_die) -> core current i(t)
Supply 2, buck:    VR2 (board, 1.8 V, stiff)   -> on-die switches -> L_mb (air-core, package trace) --------^
```

| Block | Value | Basis |
|---|---|---|
| Core load | 0.20 A average, 0.43 A peak, 0.9 V (180 mW) | BOOM replay traces (current model ±30–40%) |
| Die capacitance | 5 nF, 10 mΩ ESR | ≈ Intel's 7 nF/A at 0.43 A peak, rounded up |
| Package-to-die | 507 pH, 79.6 mΩ | Project's die-side model, kept for comparability. **Provenance unverified.** |
| Package decap | 1 µF, 100 pH ESL, 20 mΩ ESR | Typical package ceramic |
| Board trace | 2 nH, 3 mΩ | Small single-chip board |
| VR1 | 3.3→0.9 V, 2 MHz, 2.2 µH (40 mΩ), 47 µF (3 mΩ, 0.5 nH), loop crossover 100 kHz | Ordinary ~1 A point-of-load buck |
| VR2 | 1.8 V, treated as stiff | Input-rail impedance is a later check |
| Mini-buck | N phases, **L = 1–5 nH air-core** (20 mΩ DCR), switches 0.1–0.2 Ω, 0.6 ns loop delay | Buildable package-trace inductors (Intel: 59 inductors on 10 rails in 20×8 mm) |
| Switch energy | 65 pJ per phase per cycle at Ron = 0.1 Ω, scaling ∝ 1/Ron | Gate drive of a ~10 mm-wide 22 nm switch pair (~1 kΩ·µm, ~1 fF/µm, 1.8 V). **Assumption; sweep it.** |

## 1. Direct path impedance

VR1 is a closed-loop source: Z_VR1 = Z_ol / (1 + T), with Z_ol = (sL₁ + DCR₁) ∥ (1/sC₁ + ESR₁ + sESL₁) and the loop gain approximated as T = ω_c/s.

The die sees

  Z_direct = [ (Z_VR1 + sL_b + R_b) ∥ Z_pkcap ] + sL_p + R_p,   Z_die-node = Z_direct ∥ (1/sC_die + R_die).

**Result** (`impedance.png`):
- **Peak:** 1.155 Ω at 92 MHz.
- **Mid-frequency floor:** ~90–100 mΩ from 100 kHz to 10 MHz.
- **Smaller peaks:** at ~15 kHz (VR1 loop / bulk capacitor) and ~3 MHz (package decap / board trace).

**Hand check:** package L and die C give f₀ = 1/(2π√(L_pC)) = 100 MHz, Z₀ = √(L_p/C) = 318 mΩ, Q = Z₀/R = 3.55, so the peak ≈ Q·Z₀ = 1.13 Ω. This matches the model; the board decap shifts it to 92 MHz.

**What the core draws** (`load_spectrum.png`, dhrystone): 69% of the current's spectral power is above 10 MHz, and **26% lies in 50–200 MHz, right on that 1.1 Ω peak**. The droop problem is a die-capacitance / package-inductance resonance.

## 2. The mini-buck branch as an admittance

N interleaved phases share one duty command d. Per phase, the switch node averages to V_in2·d, so

  i_b = N (V_in2·d − v) / (sL + R_ph),   R_ph = Ron + DCR.

With feedback d = d₀ − G_c(s)·e^(−sT_d)·v (small signal):

  **Y_b(s) = N (1 + V_in2·G_c(s)·e^(−sT_d)) / (sL + R_ph)**                                   (3)

  Z_total = 1 / ( 1/Z_die-C + 1/Z_direct + Y_b ).

Two regimes:
- **Fixed duty (G_c = 0): a passive inductive path.** N phases of L in parallel give L/N in parallel with the package L. With 4 × 2 nH that's 0.5 nH ∥ 0.507 nH = 252 pH, so the resonance moves to 142 MHz (Z₀ = 224 mΩ). That moves the peak *away from* the 50–200 MHz load band, which is why fixed duty already helps.
- **Closed loop:** the controller's loop gain is T(s) = V_in2·G_c·e^(−sT_d)·N/(sL + R_ph)·(Z_die-C ∥ Z_direct). The type-3 compensator was searched for the lowest impedance peak subject to every crossover ≤ fsw/5, phase margin ≥ 45° and gain margin ≥ 6 dB.

**DC caveat:** two voltage sources (VR1 and a fixed-duty buck) in parallel fight at DC. A slow current-sharing loop must drive the buck's average current to zero, as in our case D. That makes Y_b → 0 at DC, so the low-frequency part of the "fixed duty" curve below ~10 kHz isn't physical.

## 3. Droop on the real trace (linear superposition, LTI)

v(t) = −IFFT{ Z_total(f) · I(f) }, with I(f) the spectrum of the dhrystone current. Policy A:

| Configuration | p99.9 droop | max |
|---|---|---|
| Direct path only | 89.1 mV | 134 mV |
| + 1-phase 2 nH buck, fixed duty | 89.2 | 130 |
| + 4-phase 2 nH, fixed duty | 75.6 | 137 |
| + 4-phase, closed loop, fsw 140 MHz (crossover ≤ 28 MHz) | 75.7 | 137 (loop too slow to matter) |
| + 1-phase, closed loop, fsw 500 MHz (crossover ≤ 100 MHz) | 69.7 | 99 |
| **+ 4-phase, closed loop, fsw 500 MHz** | **57.8** | **92** |

A continuously running mini-buck therefore needs **~500 MHz switching** before its loop helps, because the problem sits at 100 MHz.

## 4. Ripple (large-signal, periodic)

Per-phase triangle ripple: ΔI = V_in2·D(1−D) / (L·fsw). The Fourier amplitude of harmonic k is

  a_k = ΔI · 2|sin(πkD)| / (π²k²D(1−D)).

Interleaving N phases multiplies harmonic k by Σ_j e^(−i2πkj/N), which is N if N divides k and 0 otherwise. Die ripple = Σ_{N|k} N·a_k·|Z_die-C ∥ Z_direct|(k·fsw).

At 2 nH and 160 MHz the phase ripple is 1.40 A peak-to-peak, **7× the load current**.

| | D = 0.48 | 0.50 | 0.52 |
|---|---|---|---|
| N = 1 | 697 mV pk-pk | 687 | 697 |
| N = 2 | 24 mV | **0** (exact cancellation at D = 1/2) | 24 |
| N = 4 | 11 mV | **0** | 11 |

A single phase is impossible. Two or more phases cancel the output ripple, but **not the circulating ripple current**, which is what sets the loss below.

## 5. Loss of a continuously running mini-buck

Per phase: I_rms² = (I_dc/N)² + ΔI²/12. Loss:

  P = N [ I_rms² (Ron + DCR) + fsw · E_sw(Ron) ].

Minimised over fsw (≤ 500 MHz) and Ron (≤ 0.2 Ω, so it can deliver 0.5 A with ≤ 100 mV drop), **at zero DC current**, i.e. pure overhead:

| L per phase | N = 1 | N = 2 | N = 4 |
|---|---|---|---|
| 1 nH | 31 mW (17% of core) | 62 mW (35%) | 124 mW (69%) |
| 2 nH | 19 mW (10%) | 38 mW (21%) | 75 mW (42%) |
| 5 nH | 10 mW (6%) | 20 mW (11%) | 41 mW (23%) |

**A continuously switching mini-buck sized for this core costs 10–70% of the core's power, before it delivers anything.** The configuration that cut droop most (4 phases at 500 MHz) is the expensive end. That's the quantitative reason small domains use LDOs or switched-capacitor converters, and why Intel shares one big FIVR across a domain.

## 6. The alternative: on-demand injector

Keep both switches off between events. When a current step is coming, turn the high side on:
- the inductor current ramps at **di/dt = (V_in2 − V_DD)/L** (0.45 A/ns at 2 nH);
- it holds about one cycle;
- it then freewheels down through the low side, delivering its stored energy to the die.

| L | 0 → 0.25 A ramp | down-ramp |
|---|---|---|
| 1 nH | 0.28 ns | 0.28 ns |
| 2 nH | 0.56 ns | 0.56 ns |
| 5 nH | 1.39 ns | 1.39 ns |

Energy per event ≈ two switch transitions (130 pJ) + conduction (~5 pJ) ≈ **135 pJ**. Power depends on how often it fires (dhrystone, steps per 1 ns cycle):

| Fire on steps larger than | Fraction of cycles | Power (2 nH) |
|---|---|---|
| 60 mA | 12.5% | 16.9 mW (9.4% of core) |
| 100 mA | 6.3% | 8.5 mW (4.7%) |
| **150 mA** | **0.65%** | **0.9 mW (0.5%)** |

Loss scales with the event rate, not with time. A small inductor is an advantage here because it ramps fast. The ramp (~0.5 ns) is shorter than one cycle, and the issue-stage predictor knows about 13 cycles ahead which steps are coming, so the injector can be fired **before** the step instead of reacting to the droop.

## 7. Conclusions so far

1. **The droop problem is the die-capacitance / package-inductance resonance at ~100 MHz,** and a quarter of the core's current spectrum sits on it.
2. **A continuously running mini-buck only helps if it switches at ~500 MHz with a ≤ 100 MHz loop,** and then costs 20–70% of the core's power. Rejected.
3. **A fixed-duty multi-phase mini-buck acts as a parallel inductor** that moves the resonance out of the load band (−15% p99.9). But it still pays the circulating-ripple loss and needs DC current-sharing. Rejected for the same loss reason.
4. **An on-demand, predictor-triggered injector** (1–5 nH, idle between events) costs about 0.5–9% of core power depending on how often it fires, and can be timed ahead of the step. **This is the candidate to build next.**

## 8. Next derivation steps (not done yet)

1. **Injector as a controlled current source:** pulse-shape design (amplitude, hold, release), the effect of firing early or late, and false-trigger cost (overshoot if fired into a surplus).
2. **Input-rail check:** VR2 / 1.8 V rail impedance and the on-die input capacitance needed so injector pulses don't disturb that rail.
3. **Time-domain verification in Pulsim** with these real component values: ideal switches plus the E_sw model, then real MOSFET models.
4. **Per-tile placement** on the 4×4 grid, with one injector per hot tile (ALU, LSU, FPU).
5. **Comparison** against digital-LDO and switched-capacitor slices at equal area.

**Caveats:**
- The 65 pJ switch energy is a first-principles estimate. It needs a device-level number (sweep 20–150 pJ).
- The 507 pH / 79.6 mΩ package model has unverified provenance.
- VR1's loop is approximated as an integrator.
- The droop numbers are linear (small-signal) estimates. The injector is nonlinear and has to be checked in the time domain.
