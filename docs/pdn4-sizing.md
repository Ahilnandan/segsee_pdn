# pdn4 component sizing (real units, MediumBoom tile, 45 nm, 1.1 V, 333 MHz) — written before the sweeps

Load (LACPo + cache model): mean 0.14-0.21 A, p99.9 0.38-0.47 A, worst per-cycle step ~0.16-0.28 A. Leakage 57 mA.
Die 4 mm^2 (estimate; no layout data available), intrinsic C 3 nF, MOS decap 12 nF/mm^2, thick-oxide 5 nF/mm^2.

| block | sizing | check |
|---|---|---|
| off-chip buck | 3.3 -> 1.1 V, 2 MHz, L1 2.2 uH (40 mOhm), C1 47 uF (3 mOhm), current-mode PI, fc 100 kHz, DVFS slew 20 mV/us | ripple I = 2.2 V x (1/3) / (2.2 uH x 2 MHz) = 0.17 A pk-pk; ripple V = 0.17/(8 x 2 MHz x 47 uF) = 0.2 mV; conduction 0.2^2 x 0.04 = 1.6 mW (~1 %) |
| board/package | L_b 2 nH / 3 mOhm, package decap 1 uF (100 pH, 20 mOhm), bumps 0.5 nH / 50 mOhm | 1st droop 1/(2 pi sqrt(0.5 nH x C_die)) = 65-90 MHz; 2nd ~0.6 MHz |
| on-die regulator (regulator role) | N = 2-8 phases, package chip inductors 20-200 nH (DCR 10 mOhm + 2 mOhm/nH), fsw 20-200 MHz, Ron 0.5-4 Ohm/phase, fc 2-30 MHz | 100 nH, 100 MHz, 1.8 V: ripple/phase = 0.7 x 0.61/(100 nH x 100 MHz) = 43 mA; output ripple (4 ph, K = 0.26) / (8 x 4 x 100 MHz x 2 nF) ~ 2 mV; switching 2 x 5 pJ x 100 MHz x 4 ph = 4 mW/region (dominant cost). 2 nH would give 2.1 A/phase ripple -> rejected |
| on-die injector (injector role) | 1 phase, package 5-50 nH, rail 1.3-1.8 V, Ron 0.3-2 Ohm, Ipk 50-300 mA, threshold 20-60 mV | slope (Vin - 1.1)/L: 1.3 V / 20 nH = 10 mA/ns (acts at ~5-30 MHz); max current (Vin - 1.1)/(Ron + DCR) must exceed Ipk |
| parallel LDO (droop clamp) | separate 1.2-1.3 V rail through its own bumps (bump share taken from VDD), rail package decap 0.47 uF, pass R 0.1-1 Ohm, gm 1-20 A/V, tau 0.5-10 ns, Iq = 0.5 mA x (1 ns / tau) | max current (V_rail - V_die)/R_pass: 0.1 V / 0.3 Ohm = 0.33 A; loss only while clamping |
| series LDO (per-region header) | rail = V + 50-150 mV, R_pass 0.1-0.5 Ohm, gm 5-20 A/V, rail-side decap 0-0.2 mm^2 | loss ~ headroom / V = 5-14 % always; enables instant DVFS down and up-within-headroom |
| active decap | Ca 0.1-1 nF per region at 5 nF/mm^2, G 2-8, amplifier tau 0.5-2 ns, Pq 0.45 mA x G/4 x speed | |
| decap + damping | 0-0.5 mm^2 per region, series R 0-1 Ohm | optimum R ~ sqrt(L/C) ~ 0.2-0.3 Ohm |

Architectures compared: A off-chip sets V + on-die injector; B off-chip current helper + on-die regulator (role swap);
C on-die series FIVR only; D series LDO; E A + parallel LDO; plus predictive boost and the issue-queue scheduling policies.
