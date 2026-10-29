# Lab handoff (native Ubuntu)

The lab machine is Ubuntu, not Windows. Do not use WSL there. The simulator was built inside WSL on the laptop. Copy the runtime tarball onto the lab machine and run it natively.

Branch: `lab-mediumboom-traces` on `https://github.com/AspiringPianist/segsee_pdn`.

## Goal

Regenerate the 10-benchmark traces `segsee_pdn` was evaluated on, then the per-uop kernels `lacpo/hpos.npy` and `lacpo/hpos_blocks.npy`. Those files are not in git. They come from Verilator waveforms of one specific BOOM: chipyard `717c4658`, BOOM `4e9d496d`, `MediumBoomConfig`.

Do not quote the README headline (guardband 110 mV to 10 mV, +11.4% performance, −17.4% power) as confirmed. It assumes a timing model, assumes rename / ROB / issue are covered by margin, and does not charge Razor flop cost. PR #2 (`adapt-b2v2-clock-tuner`) is a draft. Its B2-vs-B2v2 table was produced on a synthetic fixture and is withdrawn. `SEGSEE_FIXTURE=1` is a smoke test. Anything labelled `SYNTHETIC_FIXTURE` is not a result.

## On the laptop, before you leave

The pack script writes `D:\WSL\segsee-runtime.tar.gz` from the WSL distro:

```
wsl -d Ubuntu-22.04 -u root -- bash /mnt/c/Users/unnat/Documents/GitHub/segsee_pdn/lab/pack_runtime.sh
```

That archive holds the simulator, `libfesvr`, the riscv-tests ELF files, the Python 3.7 + scikit-learn 0.20.4 env (`lacpo37`), and the PulSim 2.0 venv. It does not hold the chipyard source. `/opt/segsee` inside WSL is about 17 GB if a rebuild is later required. Packing that whole tree is a separate copy, not this tarball.

Also copy the `DecouplingCap` tree (at least `tools/run_lacpo_flow.py` and `BOOM/power_modeling/`, the 21 sklearn 0.20 pickles).

## On the Ubuntu lab machine

```
sudo mkdir -p /opt/segsee
sudo tar -xzf segsee-runtime.tar.gz -C /opt/segsee
git clone -b lab-mediumboom-traces https://github.com/AspiringPianist/segsee_pdn.git
```

`tar -C /opt/segsee` matches the archive layout: the pack script archives paths relative to `/opt/segsee`, so the simulator lands at

`/opt/segsee/chipyard-cy100/sims/verilator/simulator-example-MediumBoomConfig-debug`

Ubuntu 22.04 already has the libraries the binary needs (`libstdc++`, `libgcc`, `libc`). The only extra library is `libfesvr.so`, and it is in the tarball. `LD_LIBRARY_PATH` is set by the run script.

Packages to have: `tar`, `gzip`, `findutils` (`xargs`). No JDK, no sbt, no Verilator install.

## Run every benchmark on its own core

The binary was not compiled with `verilator --threads`. One simulation uses one core. `JOBS` processes at once is the parallelism that works. Verilator 4.016 does have `--threads`, but using it means regenerating the C++ model from the chipyard tree, which is not in the runtime tarball. Waveform writing is the slow part anyway: on the build laptop, about 1,000 cycles/s with no VCD and about 60 cycles/s with `-V`. A 2,000-cycle VCD was 7.5 MB, so 5,000,000 cycles is about 19 GB and many hours per benchmark. A faster CPU shortens the time. It does not shrink the files. Ten jobs at once need disk for about 190 GB and enough RAM for ten copies of the model.

```
cd segsee_pdn
export SEGSEE_TRACES=/opt/segsee/traces/lacpo_cy100/bench
export MAX_CYCLES=5000000
export JOBS=$(nproc)
bash lab/run_benches_parallel.sh
```

Shorter proof of the pipeline (not comparable to the paper):

```
MAX_CYCLES=200000 JOBS=$(nproc) bash lab/run_benches_parallel.sh
```

`memcpy.riscv` is not in this 2019 riscv-tests set. The script skips it. The other nine are: qsort, rsort, mm, dhrystone, spmv, multiply, median, vvadd, towers. `memcpy` is one of the fit benchmarks in `lacpo/kernels.py`. Add that ELF before fitting kernels, or the fit set is short one workload.

`*** FAILED *** via trace_count` in `*.sim.log` means the `-m` cap was hit. A non-empty `.vcd` is success.

## Waveforms to kernels

```
export SEGSEE_REPO=$PWD
export DECOUPLINGCAP=/path/to/DecouplingCap
export SEGSEE_TRACES=/opt/segsee/traces/lacpo_cy100/bench
export POST_JOBS=3
bash lab/run_trace_post.sh
```

`POST_JOBS=3` is the original limit because feature extraction uses a lot of RAM. Per benchmark this writes `<w>_features.csv`, `<w>_out/composed_P_t.csv`, `<w>_out/tile_P_t.csv`, `<w>_cache.npz`, `<w>_issue.npz`.

Then:

```
export SEGSEE_TRACES=/opt/segsee/traces/lacpo_cy100/bench
/opt/segsee/lacpo37/bin/python lacpo/workload.py
/opt/segsee/lacpo37/bin/python lacpo/workload5.py
```

That writes `lacpo/hpos.npy` and `lacpo/hpos_blocks.npy` next to the sources. They are gitignored. Keep them. Do not replace them with the synthetic fixture (`lacpo/tune_clock.py` only builds that fixture when `SEGSEE_FIXTURE=1`).

Only after the real kernels exist, re-run the clock comparison without the fixture:

```
/opt/segsee/lacpo37/bin/python lacpo/tune_clock.py gates
/opt/segsee/venv/bin/python pulsim_validate.py
```

Report the number you get, including a worse one. PulSim has been run once, on `restore_edge` only, and that check was redesigned after an earlier form showed the staggered restore 28.1% worse. G0 and G1 only show that the new clock code matches the old code when the new options are off.

## Do not

- Do not point the parsers at a different BOOM. `cache_events.py` and `issue_trace.py` look up fixed signal names. Chipayard 1.3.0 (`D:\chipyard` on the laptop) does not match.
- Do not treat `BOOM/synth_PD` in DecouplingCap as a netlist. Those files are ~100-byte path pointers.
- Do not rebuild the simulator on the lab machine unless the binary will not start. A rebuild needs the full `/opt/segsee/chipyard-cy100` tree, JDK 8, bison 3.5.4 (3.7+ breaks Verilator 4.016), the local ivy snapshots (`firrtl-interpreter` 1.2-SNAPSHOT is gone from Sonatype), and the `-FI` stripping in `common.mk`. That story is in `D:\WSL\scripts\build_chipyard_cy100.sh` on the laptop drive.

## After the traces exist

1. Clock comparison on the real traces (legacy stretch vs multi-level staggered restore).
2. Timing: Yosys + Nangate45 under `/opt/segsee` is only a stand-in for the missing commercial signoff. Do not replace `lacpo/timing6.py` path ratios with it and call them signoff.
3. Charge Razor flop power and hold padding. Replay covers the datapath. A later RTL hook can reuse `ldspec_miss` in this BOOM's `issue-slot.scala`.
