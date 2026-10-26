# PDN co-design on LACPo currents — pdn4 results (2026-10-22)

MediumBoom tile (2019 BOOM, chipyard 717c4658), Nangate-45-class, 1.1 V, 333 MHz. Currents: LACPo (326/326 nets) +
calibrated cache model, issue-level replay with 12 scheduling policies. Network and sizing: `pdn4-sizing.md`.
Code: `pdn_work/lacpo/` (pdn4.cpp/.py, ga4.py, ga4_validate.py, dvfs4.py, windows4.py). Data: ga4_*.json, ga4_validate.json, dvfs4.json.

Targets: droop <= 55 mV at the 99.9th percentile and <= 80 mV worst case, on ALL 10 full benchmark traces. Die 4 mm^2 (estimate).

## 1. What changed vs the earlier (pdn3) study
- Off-chip buck has its real power stage (3.3 V, 2.2 uH, 47 uF, current-mode PI 100 kHz) instead of an ideal source.
  Passive reference rises from 0.70 to **0.90 mm^2** (worst-case events are harsher with the real stage).
- Injector inductor resistance fixed: package inductors now 10 mOhm + 2 mOhm/nH (was 0.3 Ohm/nH = on-die spiral value).
- New devices: role-switchable off-chip buck (voltage setter / current helper), role-switchable on-die buck
  (injector / multiphase regulator, analytic interleaved ripple), parallel LDO droop clamp on its own rail (bumps taken
  from VDD), series LDO per region, predictive voltage boost from the issue-queue predictor.
- GA constraints now include each benchmark's 3 worst droop events for every scheduling policy.

## 2. Validated designs (all 10 full traces; area = on-die PDN area)
Passive reference: 0.90 mm^2 uniform decap (p99.9 47.7, max 78.6 mV). "Freed" = (0.90 - area) / 4 mm^2.

| arch | area mm^2 | freed die | energy+slowdown overhead | p99.9 / max mV | DVFS up / down | jitter p99.9 / guardband | key parameters |
|---|---|---|---|---|---|---|---|
| F (best area) | **0.327** | **+14.3 %** | 2.95 % | 49.7 / 78.7 | 7.7 / 5.5 us | 14.3 ps / 4.8 % | parallel LDO (1.25 V rail, 22 % of bumps, R_pass 0.34 Ohm, gm 16 A/V, tau 1.2 ns, clamp 21 mV below V), damped decap 0.23 mm^2 (R 0.2-0.7 Ohm), reorder+hold+dummy policy; injector OFF |
| B (fast DVFS) | 0.378 | +13.0 % | 2.57 % | 52.2 / 79.9 | **0.60 / 0.62 us** | 20.6 ps / 5.6 % | off-chip = current helper, on-die 4-phase regulator, 50 nH package inductors, 123 MHz, 4 Ohm/phase switches, fc 29 MHz; hold policy H10 |
| E (cheapest) | 0.394 | +12.6 % | **0.74 %** | 52.1 / 79.2 | 6.5 / 5.5 us | 21.0 ps / 5.6 % | parallel LDO (1.2 V rail, 22 % bumps, R_pass 0.42 Ohm, gm 6.4 A/V, tau 1.7 ns), decap 0.33 mm^2 |
| B | 0.416 | +12.1 % | 0.91 % | 52.6 / 80.0 | **0.53 / 0.57 us** | 19.4 ps / 5.4 % | as B, 200 nH, 35 MHz |
| F | 0.467 | +10.8 % | 0.37 % | 54.2 / 79.4 | 7.5 / 5.5 us | 21.8 ps / 6.4 % | parallel LDO + package injector (20 nH, 1.3 V) |
| B | 0.509 | +9.8 % | 0.74 % | 53.5 / 78.9 | 0.49 / 0.51 us | 17.2 ps / 4.9 % | as B, 200 nH, 25 MHz |
| E | 0.541 | +9.0 % | 0.48 % | 53.9 / 79.1 | 6.1 / 5.5 us | 18.6 ps / 5.7 % | parallel LDO, slower (2.6 ns) |
| A | 0.627 | +6.8 % | 0.05 % | 52.7 / 79.8 | 5.5 / 5.6 us | 15.4 ps / 5.3 % | damped decap + reorder; injector present but ~idle |
| D (fastest DVFS) | 0.681 | +5.5 % | 7.42 % | 46.2 / 78.3 | **0.12 / 0.12 us** | 13.4 ps / 4.9 % | series LDO, rail +71 mV, rail decap 0.23 mm^2, gm 17 A/V |
| F | 0.809 | +2.3 % | 0.06 % | 53.3 / 79.7 | 5.8 / 5.5 us | 14.1 ps / 5.3 % | damped decap + injector |

About half of the GA front designs that met the limits on the fit windows miss the 80 mV worst case by 0.1-11 mV on the
full traces (listed in ga4_validate.json); only validated-OK designs are in the table.

## 3. Findings
1. **Parallel LDO droop clamp is the most efficient lever**: it frees ~12.6 % of the die for 0.74 % energy, by
   delivering current within ~1-2 ns without an inductor. With damped decap and mild scheduling (F) it reaches the
   smallest area, 0.327 mm^2 (+14.3 % usable die) at ~3 %.
2. **Role split (off-chip as current helper, on-die multiphase buck regulating) gives ~10x faster DVFS** (0.5-0.6 us
   both ways vs 5.5 us) while still freeing 10-13 % of the die at 0.7-2.6 %. It needs large package inductors (50-200 nH);
   at 2 nH the per-phase ripple is ~2 A for a 0.1 A load. Note the GA tuned the share loop and the
   on-die loop so the board capacitors are moved by both supplies; the earlier fixed-parameter test gave 4 us.
3. **Series LDO is the only option for ~0.1 us DVFS** but costs 7-13 % energy (headroom) and rail-side decap.
4. **On-die buck as a pure injector is not worth it**: inductor slope (10-30 mA/ns) is too slow for the 65-90 MHz
   first droop; the GA kept it idle or switched it off.
5. **Not chosen by the GA**: active decap (amplifier power too high for the gain) and predictive voltage boost
   (raising V ahead of predicted droop costs more energy than the droop margin it buys). Predictor-driven scheduling
   (reorder/hold/dummy) is used only in the smallest designs.
6. **Aggressive holds can worsen the worst case** (rsort 93 -> 105 mV, memcpy 79 -> 100 mV with H4): released work
   issues as a larger burst. The GA picked gentle holds (H10) or none.
7. **Damping** (0.1-0.7 Ohm in series with decap) appears in almost every design; it flattens the 65-90 MHz resonance.

## 4. DVFS reference runs (dvfs4.json, fixed parameters)
| arch | up | down | loss |
|---|---|---|---|
| A off-chip sets V, 20 mV/us | 5.4 us | 5.5 us | 0 |
| A, 50 mV/us | 8.7 us | 14.1 us | 0 (loop overshoots/rings; controller tuning issue) |
| B role swap only during transitions | 0.98 us | 5.3 us | 5.7 % |
| D series LDO, rail follows | 3.6 us | 0.11 us | 10.3 % |
| D2 series LDO, rail pre-raised 6 us | 0.14 us | 0.11 us | 10.9 % |
| E parallel LDO | 6.2 us | 5.5 us | 3.0 % |

## 5. Caveats
- Currents: cache model's rarest steps ~40 % low -> worst-case droop optimistic; LACPo 3.6 % cycle MAE on unseen code (paper).
- Die area (4 mm^2), decap densities, package L/R, LDO speed (0.5-10 ns explored; ~1-2 ns chosen is aggressive for 45 nm),
  regulator switch-energy model (5 pJ*Ohm at 1.8 V) are assumptions; comparisons are more reliable than absolute numbers.
- On-die regulator in regulator role is an averaged model with analytic ripple; injector is switch-level.
- The off-chip loop is a simple current-mode PI; faster slew made it ring (not re-tuned).
- Prior art: parallel/series LDOs, droop clamps, hybrid on/off-chip regulation and envelope-tracking style role splits
  all exist; the contribution here is the co-design quantification on a learned-power-model workload.
