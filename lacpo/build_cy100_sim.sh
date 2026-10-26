#!/usr/bin/env bash
# Chipyard 1.0.0 MediumBoomConfig Verilator debug (VCD) simulator = LACPo's training design.
set -e
export JAVA_HOME=/home/ahilnandan/Documents/PMIC/cy100-deps/jdk8   # private JDK 8 (system Java 8 is JRE-only, no javac)
export PATH=$JAVA_HOME/bin:/home/ahilnandan/Documents/PMIC/cy100-deps/bin:$PATH
export RISCV=/home/ahilnandan/Documents/PMIC/cy100-deps/riscv-install
export LD_LIBRARY_PATH=$RISCV/lib:${LD_LIBRARY_PATH:-}
cd /home/ahilnandan/Documents/PMIC/chipyard-1.0.0/sims/verilator
java -version 2>&1 | head -1
# 2019 build lists a dead resolver host; override with live mirrors of the same artifacts (this build only)
JA="-Xmx8G -Xss8M -Dsbt.override.build.repos=true -Dsbt.repository.config=/home/ahilnandan/Documents/PMIC/cy100-deps/sbt-repositories"
make CONFIG=MediumBoomConfig debug -j16 JAVA_ARGS="$JA"
cp simulator-example-MediumBoomConfig-debug simulator-example-MediumBoomConfig-debug.boom4e9d496d; ls -la simulator-example-MediumBoomConfig-debug*
echo BUILD_DONE
