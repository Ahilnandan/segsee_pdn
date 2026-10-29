# Lab handoff (native Ubuntu, this machine's chipyard)

The lab machine is Ubuntu and already has chipyard. Use that checkout. Do not unpack `segsee-runtime.tar.gz` over it, and do not clone a second chipyard. `lab/run_benches_parallel.sh` refuses to run until `CHIPYARD_ROOT` (or `~/chipyard`) contains a MediumBoom simulator.

Branch: `lab-mediumboom-traces` on `https://github.com/AspiringPianist/segsee_pdn`.

## Goal

Regenerate the 10-benchmark traces `segsee_pdn` was evaluated on, then the per-uop kernels `lacpo/hpos.npy` and `lacpo/hpos_blocks.npy`. Those files are not in git. They come from Verilator waveforms of one specific BOOM: chipyard `717c4658`, BOOM `4e9d496d`, `MediumBoomConfig`.

Do not quote the README headline (guardband 110 mV to 10 mV, +11.4% performance, −17.4% power) as confirmed. It assumes a timing model, assumes rename / ROB / issue are covered by margin, and does not charge Razor flop cost. PR #2 (`adapt-b2v2-clock-tuner`) is a draft. Its B2-vs-B2v2 table was produced on a synthetic fixture and is withdrawn. `SEGSEE_FIXTURE=1` is a smoke test. Anything labelled `SYNTHETIC_FIXTURE` is not a result.

## What the lab agent must use

Use the chipyard already on this machine.

1. Set `CHIPYARD_ROOT` to that tree if it is not `~/chipyard`, `~/chipyard-1.0.0`, or `/opt/chipyard`.
2. Set `RISCV` if the toolchain is not `$CHIPYARD_ROOT/riscv-tools-install` or `$CHIPYARD_ROOT/.conda-env/riscv-tools`.
3. If `sims/verilator` has no `simulator-*MediumBoom*` binary, build it with **this** chipyard's make: `cd $CHIPYARD_ROOT/sims/verilator && make CONFIG=MediumBoomConfig debug -j$(nproc)`. A debug build is required because the next step reads a VCD.
4. Still copy `DecouplingCap` (`tools/run_lacpo_flow.py` and `BOOM/power_modeling/`, the sklearn 0.20 pickles). That is not part of chipyard. Unpickling those models needs scikit-learn 0.20 (`PYTHON=/path/to/that/python`).

`cache_events.py` and `issue_trace.py` look up leaves under `boom_tile` (for example `dcache.data_io_read_valid`, `int_issue_unit.issue_slots_0_valid`). The prefix they try first is `TOP.TestHarness.TestHarness.boom_tile.`. If the first VCD's hierarchy differs, change that prefix to this chipyard's path. Do not switch to a different config to avoid doing that.

LACPo's published model was trained on chipyard `717c4658` / BOOM `4e9d496d`. Print `git -C $CHIPYARD_ROOT/generators/boom rev-parse HEAD` in the run log. A different BOOM still runs; say so next to any kernel or droop number.

The laptop tarball `D:\WSL\segsee-runtime.tar.gz` is a fallback simulator only. Do not use it when this machine's chipyard builds.

## Run every benchmark on its own core

One simulation uses one core unless this chipyard's Verilator was configured with `--threads`. `JOBS` processes at once is the parallelism to use. On the build laptop, a full VCD was about 60 cycles/s and about 19 GB for 5,000,000 cycles. A faster CPU shortens the time. It does not shrink the files.

```
cd segsee_pdn
export CHIPYARD_ROOT=/path/to/this/machines/chipyard
export RISCV=/path/to/this/machines/riscv
export SEGSEE_TRACES=$HOME/segsee_traces/lacpo_cy100/bench
export MAX_CYCLES=5000000
export JOBS=$(nproc)
bash lab/run_benches_parallel.sh
```

Shorter proof of the pipeline (not comparable to the paper):

```
MAX_CYCLES=200000 JOBS=$(nproc) bash lab/run_benches_parallel.sh
```

The script skips any benchmark whose `.riscv` file is missing. `memcpy` is one of the fit benchmarks in `lacpo/kernels.py`. If this toolchain has no `memcpy.riscv`, build that ELF with this machine's `RISCV` before fitting kernels.

`*** FAILED *** via trace_count` in `*.sim.log` means the `-m` cap was hit. A non-empty `.vcd` is success.

## Waveforms to kernels

```
export SEGSEE_REPO=$PWD
export DECOUPLINGCAP=/path/to/DecouplingCap
export SEGSEE_TRACES=$HOME/segsee_traces/lacpo_cy100/bench
export PYTHON=/path/to/python-with-sklearn-0.20
export POST_JOBS=3
bash lab/run_trace_post.sh
```

`POST_JOBS=3` is the original limit because feature extraction uses a lot of RAM. Per benchmark this writes `<w>_features.csv`, `<w>_out/composed_P_t.csv`, `<w>_out/tile_P_t.csv`, `<w>_cache.npz`, `<w>_issue.npz`.

Then, with the same `SEGSEE_TRACES` and `PYTHON`:

```
"$PYTHON" lacpo/workload.py
"$PYTHON" lacpo/workload5.py
```

That writes `lacpo/hpos.npy` and `lacpo/hpos_blocks.npy` next to the sources. They are gitignored. Keep them. Do not replace them with the synthetic fixture (`lacpo/tune_clock.py` only builds that fixture when `SEGSEE_FIXTURE=1`).

Only after the real kernels exist, re-run the clock comparison without the fixture:

```
"$PYTHON" lacpo/tune_clock.py gates
python3 pulsim_validate.py
```

Report the number you get, including a worse one. PulSim has been run once, on `restore_edge` only, and that check was redesigned after an earlier form showed the staggered restore 28.1% worse. G0 and G1 only show that the new clock code matches the old code when the new options are off.

## Do not

- Do not clone a fresh chipyard and do not use the laptop runtime tarball when this machine already has one.
- Do not treat `BOOM/synth_PD` in DecouplingCap as a netlist. Those files are ~100-byte path pointers.
- If the first VCD is missing `boom_tile` signals, fix the prefix in `lacpo/cache_events.py` and `lacpo/issue_trace.py`. Do not swap in a different simulator to hide that.

## After the traces exist

1. Clock comparison on the real traces (legacy stretch vs multi-level staggered restore).
2. Timing: Yosys + Nangate45 under `/opt/segsee` is only a stand-in for the missing commercial signoff. Do not replace `lacpo/timing6.py` path ratios with it and call them signoff.
3. Charge Razor flop power and hold padding. Replay covers the datapath. A later RTL hook can reuse `ldspec_miss` in this BOOM's `issue-slot.scala`.
