#!/usr/bin/env bash
# Chipyard 1.0.0 submodules for MediumBoomConfig RTL sim only: generators/tools (recursive), no toolchains,
# no hammer plugins, no FireSim (its tiny midas.targetutils dep gets a local stub if the build needs it).
set -e
# GitHub disabled git:// in 2022: rewrite for this script only (no global git config change)
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=url.https://github.com/.insteadOf GIT_CONFIG_VALUE_0=git://github.com/
cd /home/ahilnandan/Documents/PMIC/chipyard-1.0.0
for name in toolchains/*/*/ ; do git config submodule."${name%/}".update none; done
git config submodule.toolchains/qemu.update none 2>/dev/null || true
git config submodule.sims/firesim.update none
for p in vlsi/hammer-cadence-plugins vlsi/hammer-synopsys-plugins vlsi/hammer-mentor-plugins; do git config submodule.$p.update none; done
git submodule update --init --recursive
# riscv-isa-sim only (for libfesvr) from the toolchain submodules
git config --unset submodule.toolchains/riscv-tools/riscv-isa-sim.update || true
git submodule update --init toolchains/riscv-tools/riscv-isa-sim
echo INIT_DONE
