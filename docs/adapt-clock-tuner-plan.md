# Plan: ADAPT-style adaptive-clock entry/exit tuning (policy6 B2 -> B2v2)

Written 2026-10-27 before any code. Source paper: Radke and Roy (Intel), "ADAPT - Adaptive Droop and Performance
Tuner for Droop Mitigation Schemes", EPEPS 2025, IEEE 11346569. Only the abstract and reference list were readable
(IEEE blocked the download from this machine); the full PDF must be read before 2.3 is finalised. The second paper
(IEEE 10689633) is not in Crossref, OpenAlex or Semantic Scholar and is left out of this plan until the PDF is on disk.

## 0. What ADAPT says (abstract-level)
- Adaptive clock modulation compares supply noise with one or more voltage thresholds and modulates the clock to
  soften dI/dt.
- Entry into and exit from the modulated profile follow specific sequences ("staggered restore of the clock
  frequency") because the PDN and Icc(t) interact; the sequence trades performance against droop.
- An automated tuner searches the settings and returns the best mitigation configuration. Its references are the
  adaptive-clock silicon papers (AMD Steamroller 2014, POWER9 2017, 16 nm auto-calibrating clock 2015, Qualcomm 7 nm
  Hexagon 2020, 80 ns fast-lock clock 2021), Intel's 2021 PDN review, and two Bayesian-optimisation papers, so the
  tuner is most likely BO over a simulated droop/performance score.

## 1. Where the repo stands (numbers from eval6_A6.json / eval6_B8.json, 10 full traces)
B2 today (`lacpo/policy6.py`): one detector threshold `cs_thr` on the worst-tile voltage delayed `cs_lat` = 2 cycles,
fixed stretch `cs_len` = 4 periods by a fixed amount `cs_s` = 0.25, **instant restore** when `str_left` reaches 0.
Stretched cycles are exempt from timing failure.

| design | G | config | worst max (mV) | mean p99.9 | errors / 1M uops | mean / worst slowdown |
|---|---|---|---|---|---|---|
| banks | 40 | razor | 104.6 | 46.1 | 391 | 0 / 0 |
| banks | 40 | stretch (thr 30) + razor | 81.4 | 45.3 | 251 | 2.07 % / 7.36 % |
| banks | 30 | razor | 96.1 | 46.1 | 1430 | 0 / 0 |
| banks | 30 | stretch (thr 20) + razor | 77.2 | 44.5 | 607 | 5.31 % / 12.64 % |
| decap | 40 | razor | 110.4 | 51.9 | 487 | 0 / 0 |
| decap | 40 | stretch (thr 30) + razor | 94.9 | 49.1 | 279 | 1.86 % / 6.94 % |
| banks | none | A (no policy) | 85.7 | 46.4 | - | 0 / 0 |
| banks | none | A + stretch thr 40 | 81.5 | 46.1 | - | 0.71 % / 2.45 % |

Reading: the stretch is the only mechanism in the repo that cuts **worst max** droop by ~20 mV (banks, G 40), but it
costs 2-5 % mean slowdown, which is more than the guardband it saves is worth (alpha-power law: 10 mV ~ 0.96 % delay).
Selective replay holds a 10 mV guardband at 0.02 % slowdown, so the stretch cannot compete in the replay regime.
It can matter where replay is absent, where recovery is a full flush (12 % slowdown at G = 10), and on the shared
multi-core die where the droop is one chip-wide dI/dt event.

Why the current model should be improvable:
- Instant restore is a current step of up to 1 - 1/(1 + cs_s) = 20 % of the dynamic current, in the 77 MHz band
  (0.24 cycles/cycle). The repo already saw this failure mode once: pulse-masking stretch rang the PDN (p99.9 57.8 ->
  83.1 mV) and was replaced by phase stretch. The restore edge was never looked at.
- stretch_time ~ nstretch x cs_len x cs_s (G 40 banks: 3617 x 4 x 0.25 ~ 3617 vs 4238 observed), so most of the
  slowdown is the fixed 4-cycle x 25 % per trigger whatever the droop depth. A shallow level for shallow droops and a
  deep level only for deep droops is exactly the ADAPT multi-threshold.

## 2. Work packages
### 2.1 B2v2 clock-modulation state machine (policy6, new params, defaults OFF)
All new parameters default to the old behaviour so the existing gates (defaults == closedloop6; `cs_on=1, cs_thr` ==
old B2 bit-identical) keep passing.
- `cs_levels`: list of (threshold mV below nominal, stretch s, min dwell cycles). Deeper threshold -> larger s. The
  detector picks the deepest level whose threshold is crossed; a deeper trigger pre-empts a shallower one.
- `cs_exit`: restore ladder, list of (s, cycles), walked after the dwell has elapsed AND the delayed sensed voltage has
  stayed above (threshold + `cs_hyst` mV) for `cs_rec` cycles. Default `[]` = instant restore (today).
- `cs_entry`: optional one-cycle half-step on entry (ADAPT names entry sequencing too; keep it as a knob, expect ~0).
- Per-cycle stretch `s_t` replaces the global `prm['cs_s']` at the two places it is used (current spread, time
  accumulator). The exported `stretch_time` stays the sum of s_t.
- Timing model: a stretched cycle is no longer blanket-exempt. It fails only if vmin < VF - s_t / 0.00096 mV
  (alpha-power law, 0.96 % delay per 10 mV). For s >= 0.1 this equals the old exemption; for the small restore steps
  (s = 1/16 -> ~65 mV headroom) it is the honest version.
- New per-run metrics: `n_restore_droop` (worst-tile dip > G/2 within 8 cycles after a restore step), `stretch_hist`
  (time at each level), `mod_entries`.

### 2.2 Tuner (new `lacpo/tune_clock.py`)
- Search space: thresholds (2 levels, 10-50 mV), s per level (0.0625-0.375), dwell 1-8, exit ladder (1-4 steps,
  1-6 cycles each), hyst 0-10 mV, rec 1-4.
- Score, per candidate: worst max over the 5 GA windows (ga5 convention), mean p99.9, mean slowdown, errors at G
  (when timing is on). Multi-objective, so reuse the ga5 NSGA-II loop first; add a BO variant (GP on the scalarised
  score) only if the GA needs > ~400 evaluations to converge. ADAPT's novelty is the tuner, so the deliverable is the
  Pareto front plus the chosen configs, not a single number.
- Three regimes, each with its own front: (R1) no replay, banks and decap designs; (R2) selective replay at G = 10 /
  20 / 30; (R3) flush recovery at G = 10 / 20 / 30.
- Validate the front picks on all 10 full traces (the ga5 -> val5 pattern). Per-bench max shifts +-10-20 mV from any
  decision change, so decisions use p99.9 and 10-bench worst max only.

### 2.3 Multi-core (mc_engine)
Add the same B2v2 model as a chip-wide clock (one detector on the die's worst tile, one state machine, all four cores
stretched together) and compare against the cross-core hold coordinator on the three mc scenarios (lockstep 4x qsort,
barrier jitter, 4 different programs) at equal mean slowdown. Reference numbers: lockstep p99.9 48.3 -> 22.9 mV at
0.19 % with holds. A chip-wide clock is what the ADAPT silicon references actually build, so this is the most
realistic head-to-head.

### 2.4 PulSim check of the restore edge
Repo already cross-validates the C++ PDN solver against PulSim (`pulsim_validate.py`). Add one case: a 300 mA step
held for 4 cycles then restored (a) instantly, (b) in 4 x 1/16 steps. Confirm the second droop in the C++ solver
matches PulSim within the existing 5 % peak tolerance before trusting the `n_restore_droop` metric.

## 3. Gates and success criteria (repo style: negative results are reported, not hidden)
- G0 defaults unchanged: policy6 with all defaults == closedloop6 bit-identical (existing gate).
- G1 back-compat: `cs_levels=[(cs_thr, 0.25, 4)], cs_exit=[]` == old B2 bit-identical on towers and qsort 20k windows.
- G2 restore edge: synthetic step test shows the instant-restore second droop; the ladder removes >= 50 % of it.
- G3 forward progress: a run never stays modulated > `cs_max_time` cycles continuously (ADAPT-style floor).
- Success (any one, otherwise the plan is a documented negative result):
  - S1 replay-free regime, banks: 10-bench worst max <= 81 mV (today's stretch result) at mean slowdown <= 1.0 %
    (today 2.07 %), or guardband gain >= 10 mV at <= 0.5 % slowdown (the alpha-power break-even is ~0.1 % per mV).
  - S2 flush regime, G = 10-20: mean slowdown below the 12.3 % flush result by >= 3x with no collapse on rsort/mm/qsort.
  - S3 multi-core lockstep: p99.9 <= 25 mV at <= 0.2 % slowdown with the chip-wide clock alone, or a combined
    clock + coordinator point that dominates the coordinator-only point.

## 4. Environment (nothing in 2.x can run on this Windows checkout today)
- Missing, gitignored artifacts: `traces/lacpo_cy100/bench/*` (LACPo per-tile currents), `lacpo/*_cache.npz`,
  `hpos.npy`, `tile_leak_mA.npy`, forecast7 weight files, `libpdn6.so` / `libpdn7.so`. `replay3.B` points at
  `/home/ahilnandan/...`. These come from the teammate's machine or are regenerated with LACPo + the BOOM VCD/PTPX
  runs, which are not in this repo. First action: request the artifact bundle (traces + npy/npz + weights).
- Linux side: WSL2 Ubuntu with its vhdx on D: (C: has 25 GB free, D: 145 GB): `wsl --install -d Ubuntu`, then
  `wsl --manage Ubuntu --move D:\WSL\Ubuntu` (or `--location D:\WSL\Ubuntu` on WSL >= 2.4). Build the solvers with
  `g++ -O3 -shared -fPIC -o libpdn6.so pdn6.cpp` (and pdn7); `pip install pulsim` in the WSL venv (`import pulsim as p`
  is what the repo uses). Make `replay3.B` an env var (`SEGSEE_TRACES`) so the checkout runs anywhere; that is a
  one-line change and belongs in the first PR.
- Compute: eval6 runs 10 benches x ~20 configs; a tuner adds a few hundred window runs. Fine on 6 workers in WSL.

## 5. Order of work and deliverables
1. PR-1 (small, mergeable now): `replay3.B` from env var; `.gitignore` note listing which artifacts a fresh clone
   needs; this plan.
2. PR-2: B2v2 state machine + metrics + gates G0-G3 + PulSim restore-edge case.
3. PR-3: `tune_clock.py`, fronts for R1-R3, full-trace validation table appended to `protocol-results.md`.
4. PR-4: chip-wide clock in `mc_engine` + the head-to-head against the coordinator.
Repository access is READ for this account, so each PR comes from a fork; the plan itself is filed as an issue.

## 6. Implementation & Gate Results (2026-10-27)
Implemented in `lacpo/policy6.py` (`ClockMod`, `sim_clock_openloop`), `lacpo/mc_engine.py` (`run_mc(..., clock=...)`),
`pulsim_validate.py` (`restore_edge_check`), and `lacpo/tune_clock.py` (`run_gates`, `benchmark_b2_vs_b2v2`, `tune_nsga2`),
with native `libpdn5.so` / `libpdn6.so` compiled via `g++ -O3 -shared -static` against Eigen 5.0.1 and `LEAK_T`
auto-derived from `floorplan.json` when `tile_leak_mA.npy` is absent.

### 6.0 Status: no B2 vs B2v2 result yet
An earlier version of this section reported B2 vs B2v2 numbers. They are withdrawn:
- They were produced on a **synthetic fixture**, not on the 10 real traces. The per-uop kernels (`hpos*.npy`) and
  the uop streams were made up in `tune_clock.ensure_fixture_env`, and the table mixed outputs from two fixture runs.
- The `hello_vcd` trace has constant placeholder currents for the cache blocks, so it is not a realistic workload.
- The restore-edge check was redesigned after its first form showed the ladder 28.1 % *worse* than instant restore.
- PulSim was not installed. `pulsim_validate.py restore_edge` ran the Python series-RLC model, not PulSim.

The fixture is now opt-in (`SEGSEE_FIXTURE=1`), writes only to a temp directory, and is labelled `SYNTHETIC_FIXTURE`.
Its numbers are smoke tests, not results.

### 6.1 What is established
- **G0 (`defaults == closedloop6`)** and **G1 (`cs_levels=[(cs_thr, 0.25, 4)], cs_exit=[] == legacy B2`)** are
  code-equivalence checks: with the new parameters off, `policy6` reproduces the old code bit for bit. These hold on any
  input, so the fixture is adequate for them.
- **G3**: `cs_max_time` bounds continuous modulation (code property).
- **G2** (restore-edge suppression) and every B2 vs B2v2 comparison need the real traces and `hpos*.npy` (Issue #1),
  plus a PulSim cross-check, before any number is reported.

