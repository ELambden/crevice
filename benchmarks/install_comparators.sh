#!/usr/bin/env bash
# Install the benchmark comparators into <workspace>/tools (nothing is vendored).
#
#   bash benchmarks/install_comparators.sh [WORKSPACE]
#
# WORKSPACE defaults to the first entry of $CREVICE_WORKSPACE_DIR, else ./.crevice.
# Needs: curl, git, make, a C compiler, unzip and conda (miniconda/miniforge) for the
# Fortran compiler, Java and Python 2.7 environments. Each step is skipped when its
# result already exists. Versions and licences (check them before redistributing
# anything; the benchmark only downloads and runs the tools):
#
#   HOLE 2.3.1      github.com/osmart/hole2 tag v2.3.1        Apache-2.0
#   CAVER 3.0.3     caver.cz (caver_3.0.3.zip)                GNU GPL-3.0 (libraries LGPL-3/GPL-3)
#   OpenJDK 17      conda-forge openjdk=17 (to run CAVER)     GPL-2.0 with Classpath exception
#   fpocket 4.1.5   github.com/Discngine/fpocket tag 4.1.5    MIT (mdpocket; run it through a short path: a long argv[0] makes it abort)
#   pyKVFinder      PyPI pyKVFinder (0.9.5 when benchmarked)  GPL-3.0-or-later
#   POVME 3.0.34    github.com/POVME/POVME (Python 2 only)    MIT; run with conda python=2.7
#   mdahole2        PyPI mdahole2 (0.5.0 when benchmarked)    LGPLv2+; runs HOLE per trajectory frame
#
# Then export, for volumes_vs_comparators.py:
#   CREVICE_KVFINDER_PYTHON=$WS/tools/kvfinder-venv/bin/python
#   CREVICE_POVME_PYTHON=$WS/tools/py27/bin/python2.7
#   CREVICE_POVME_SCRIPT=$WS/tools/povme-src/POVME/POVME3.py
set -euo pipefail
WS="${1:-${CREVICE_WORKSPACE_DIR%%:*}}"
WS="${WS:-.crevice}"
T="$WS/tools"
mkdir -p "$T"
CONDA="${CONDA_EXE:-conda}"

# HOLE 2.3.1 (Fortran; gfortran from conda-forge in an isolated prefix)
if [ ! -x "$T/hole2/bin/hole" ]; then
  mkdir -p "$T/hole2"
  curl -L -o "$T/hole2/hole2-v2.3.1.tar.gz" https://github.com/osmart/hole2/archive/refs/tags/v2.3.1.tar.gz
  tar -xzf "$T/hole2/hole2-v2.3.1.tar.gz" -C "$T/hole2"
  "$CONDA" create -y -p "$T/hole2/toolchain" -c conda-forge --override-channels gfortran make
  # Serial make: a parallel make races at link time.
  (cd "$T/hole2/hole2-2.3.1/src" && export PATH="$(cd "$T/hole2/toolchain/bin" && pwd):$PATH" \
     && source ../source.apache && make && make PREFIX="$(cd "$T/hole2" && pwd)" install)
fi

# OpenJDK and CAVER 3.0.3
[ -x "$T/openjdk/bin/java" ] || "$CONDA" create -y -p "$T/openjdk" -c conda-forge openjdk=17
if [ ! -f "$T/caver/caver_3.0.3/caver/caver.jar" ]; then
  mkdir -p "$T/caver"
  curl -L -o "$T/caver/caver_3.0.3.zip" https://www.caver.cz/fil/download/caver30/303/caver_3.0.3.zip
  (cd "$T/caver" && unzip -q -o caver_3.0.3.zip)
fi

# fpocket / mdpocket
if [ ! -x "$T/fpocket/bin/mdpocket" ]; then
  git clone --depth 1 --branch 4.1.5 https://github.com/Discngine/fpocket.git "$T/fpocket-src"
  (cd "$T/fpocket-src" && make)        # serial: a parallel make races on the qhull headers
  mkdir -p "$T/fpocket/bin" && cp "$T/fpocket-src/bin/"* "$T/fpocket/bin/"
fi

# pyKVFinder in its own venv
if [ ! -x "$T/kvfinder-venv/bin/python" ]; then
  python3 -m venv "$T/kvfinder-venv"
  "$T/kvfinder-venv/bin/pip" install pyKVFinder numpy
fi

# mdahole2 (HOLE on every trajectory frame, for timing_trajectory.py) in its own venv
if [ ! -x "$T/mdahole2-venv/bin/python" ]; then
  python3 -m venv "$T/mdahole2-venv"
  "$T/mdahole2-venv/bin/pip" install mdahole2
fi

# POVME 3 (Python 2.7)
[ -d "$T/povme-src" ] || git clone --depth 1 https://github.com/POVME/POVME.git "$T/povme-src"
[ -x "$T/py27/bin/python2.7" ] || "$CONDA" create -y -p "$T/py27" -c conda-forge python=2.7 numpy scipy

echo "comparators installed under $T"
